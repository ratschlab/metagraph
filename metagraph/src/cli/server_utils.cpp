#include <zlib.h>
#include <json/json.h>
#include <server_http.hpp>

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "server_utils.hpp"


namespace mtg {
namespace cli {

using mtg::common::logger;

/**
 * Compress a STL string using zlib with given compression level and return
 * the binary data.
 * Source: https://panthema.net/2007/0328-ZLibString.html
 */
std::string compress_string(const std::string &str,
                            int compressionlevel = Z_BEST_COMPRESSION,
                            bool gzip = false) {
    z_stream zs; // z_stream is zlib's control structure
    memset(&zs, 0, sizeof(zs));

    // zlib's deflate stream by default; a gzip container (window bits + 16) when the
    // client asked for gzip, which is what most HTTP clients send by default
    const int window_bits = 15 + (gzip ? 16 : 0);
    if (deflateInit2(&zs, compressionlevel, Z_DEFLATED, window_bits, 8,
                     Z_DEFAULT_STRATEGY) != Z_OK)
        throw std::runtime_error("deflateInit failed while compressing.");

    zs.next_in = (Bytef *)(str.data());
    zs.avail_in = str.size(); // set the z_stream's input

    int ret;
    char outbuffer[32768];
    std::string outstring;

    // retrieve the compressed bytes blockwise
    do {
        zs.next_out = reinterpret_cast<Bytef *>(outbuffer);
        zs.avail_out = sizeof(outbuffer);

        ret = deflate(&zs, Z_FINISH);

        if (outstring.size() < zs.total_out) {
            // append the block to the output string
            outstring.append(outbuffer, zs.total_out - outstring.size());
        }
    } while (ret == Z_OK);

    deflateEnd(&zs);

    if (ret != Z_STREAM_END) { // an error occurred that was not EOF
        std::ostringstream oss;
        oss << "Exception during zlib compression: (" << ret << ") " << zs.msg;
        throw std::runtime_error(oss.str());
    }

    return outstring;
}

// The content encoding the client accepts, "gzip" preferred over "deflate", or "" when
// it accepts neither (the response is then sent uncompressed).
std::string requested_encoding(const std::shared_ptr<HttpServer::Request> &request) {
    auto encoding_header = request->header.find("Accept-Encoding");
    if (encoding_header == request->header.end())
        return "";
    const std::string &accept = encoding_header->second;
    if (accept.find("gzip") != std::string::npos)
        return "gzip";
    if (accept.find("deflate") != std::string::npos)
        return "deflate";
    return "";
}

bool is_compression_requested(const std::shared_ptr<HttpServer::Request> &request) {
    return !requested_encoding(request).empty();
}

Json::Value parse_json_string(const std::string &msg) {
    Json::Value json;

    Json::CharReaderBuilder rbuilder;
    std::unique_ptr<Json::CharReader> reader { rbuilder.newCharReader() };
    std::string errors;

    if (!reader->parse(msg.data(), msg.data() + msg.size(), &json, &errors))
        throw std::invalid_argument("Bad json received: " + errors);

    return json;
}

std::string json_str_with_error_msg(const std::string &msg) {
    Json::Value root;
    root["error"] = msg;
    return Json::writeString(Json::StreamWriterBuilder(), root);
}

void process_request(std::shared_ptr<HttpServer::Response> &response,
                     const std::shared_ptr<HttpServer::Request> &request,
                     size_t request_id,
                     const std::function<Json::Value(const std::string &)> &process,
                     bool compact) {
    logger->info("[Server] {} request {} from {}", request->path, request_id,
                 request->remote_endpoint().address().to_string());
    Timer timer;
    // Retrieve string:
    std::string content = request->content.string();
    SimpleWeb::CaseInsensitiveMultimap header({ { "Content-Type", "application/json" } });
    SimpleWeb::StatusCode status;
    std::string ret;

    try {
        // Return JSON string
        status = SimpleWeb::StatusCode::success_ok;
        Json::StreamWriterBuilder builder;
        if (compact)
            builder["indentation"] = "";   // the traversal routes: half the bytes of the indented form
        ret = Json::writeString(builder, process(content));
        const std::string encoding = requested_encoding(request);
        if (!encoding.empty()) {
            ret = compress_string(ret, Z_BEST_COMPRESSION, encoding == "gzip");
            header.insert(std::make_pair("Content-Encoding", encoding));
            header.insert(std::make_pair("Content-Length", std::to_string(ret.size())));
        }
    } catch (const CurrentlyInitializingError& e) {
        logger->info("[Server] Got a request during initialization. Asked to come back later");
        status = SimpleWeb::StatusCode::server_error_service_unavailable;
        header.insert(std::make_pair("Retry-After", "60")); // ask to come back in 60 seconds
        ret = json_str_with_error_msg("Server is currently initializing, please come back later.");
    } catch (const std::exception& e) {
        logger->warn("[Server] Error on request {}: {}", request_id, e.what());
        status = SimpleWeb::StatusCode::client_error_bad_request;
        ret = json_str_with_error_msg(e.what());
    } catch (...) {
        logger->warn("[Server] Error on request {}", request_id);
        status = SimpleWeb::StatusCode::server_error_internal_server_error;
        ret = json_str_with_error_msg("Internal server error");
    }
    double processing_time = timer.elapsed();
    response->write(status, ret, header);
    logger->info("[Server] Request {} processing time: {:.3f} sec, response size: {:.1f} KB, "
                 "finished in {:.3f} sec",
                 request_id, processing_time, (double)ret.size() / 1000, timer.elapsed());
}

} // namespace cli
} // namespace mtg
