#include <zlib.h>
#include <json/json.h>
#include <server_http.hpp>

#include "common/logger.hpp"
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

void write_response(SimpleWeb::StatusCode status,
                    const std::string &msg,
                    std::shared_ptr<HttpServer::Response> response,
                    bool compress) {
    auto header = SimpleWeb::CaseInsensitiveMultimap(
            { { "Content-Type", "application/json" } });

    std::string msg_send = msg;
    if (compress) {
        msg_send = compress_string(msg);
        header.insert(std::make_pair("Content-Encoding", "deflate"));
        header.insert(std::make_pair("Content-Length", std::to_string(msg_send.size())));
    }

    response->write(status, msg_send, header);
}

Json::Value parse_json_string(const std::string &msg) {
    Json::Value json;

    Json::CharReaderBuilder rbuilder;
    std::unique_ptr<Json::CharReader> reader { rbuilder.newCharReader() };
    std::string errors;

    if (!reader->parse(msg.data(), msg.data() + msg.size(), &json, &errors))
        throw std::domain_error("Bad json received: " + errors);

    return json;
}

std::string json_str_with_error_msg(const std::string &msg) {
    Json::Value root;
    root["error"] = msg;
    return Json::writeString(Json::StreamWriterBuilder(), root);
}

void process_request(std::shared_ptr<HttpServer::Response> &response,
                     const std::shared_ptr<HttpServer::Request> &request,
                     const std::function<std::string(const std::string &)> &process) {
    // Retrieve string:
    std::string content = request->content.string();
    logger->info("[Server] {} request from {}", request->path,
                 request->remote_endpoint().address().to_string());

    try {
        std::string ret = process(content);
        write_response(SimpleWeb::StatusCode::success_ok, ret, response,
                       is_compression_requested(request));
    } catch (const std::exception &e) {
        logger->info("[Server] Error on request\n{}", e.what());
        response->write(SimpleWeb::StatusCode::client_error_bad_request,
                        json_str_with_error_msg(e.what()));
    } catch (...) {
        logger->info("[Server] Error on request");
        response->write(SimpleWeb::StatusCode::server_error_internal_server_error,
                        json_str_with_error_msg("Internal server error"));
    }
}

} // namespace cli
} // namespace mtg
