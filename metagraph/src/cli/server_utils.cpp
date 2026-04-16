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
                            int compressionlevel = Z_BEST_COMPRESSION) {
    z_stream zs; // z_stream is zlib's control structure
    memset(&zs, 0, sizeof(zs));

    if (deflateInit(&zs, compressionlevel) != Z_OK)
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

bool is_compression_requested(const std::shared_ptr<HttpServer::Request> &request) {
    auto encoding_header = request->header.find("Accept-Encoding");
    return encoding_header != request->header.end()
            && encoding_header->second.find("deflate") != std::string::npos;
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
                     const std::function<Json::Value(const std::string &)> &process) {
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
        ret = Json::writeString(Json::StreamWriterBuilder(), process(content));
        if (is_compression_requested(request)) {
            ret = compress_string(ret);
            header.insert(std::make_pair("Content-Encoding", "deflate"));
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
    double transfer_time = timer.elapsed() - processing_time;
    logger->info("[Server] Request {} processing time: {:.3f} sec, response size: {:.1f} KB, "
                 "transfer time: {:.3f} sec, finished in {:.3f} sec", request_id,
                 processing_time, (double)ret.size() / 1000, transfer_time, timer.elapsed());
}

} // namespace cli
} // namespace mtg
