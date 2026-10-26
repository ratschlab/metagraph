#include "server_checks.hpp"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <limits>
#include <memory>
#include <ostream>
#include <sstream>
#include <streambuf>
#include <zlib.h>

#include "common/unix_tools.hpp"
#include "common/logger.hpp"

namespace mtg::cli {
using common::logger;
namespace {
class DeliveryChecks {
  public:
    DeliveryChecks(const std::string &out, const std::function<void()> &check, size_t interval,
                   double *max_gap_ms)
          : out_(out), check_(check), interval_(interval), next_(out.size() + interval),
            max_gap_ms_(max_gap_ms) {}
    // the bytes |out| may grow by before the next check is due (at least 1)
    size_t room() const { return next_ > out_.size() ? next_ - out_.size() : 1; }
    void tick() {
        if (out_.size() < next_)
            return;
        next_ = out_.size() + interval_;
        gap();
        check_();
        since_.reset();
    }
    void finish() { gap(); }

  private:
    void gap() {
        if (max_gap_ms_)
            *max_gap_ms_ = std::max(*max_gap_ms_, since_.elapsed() * 1000);
    }
    const std::string &out_;
    const std::function<void()> &check_;
    size_t interval_;
    size_t next_;
    double *max_gap_ms_;
    Timer since_;
};

// A streambuf appending to a string under DeliveryChecks: the JSON text of a large response is
// written under the caller's check, byte for byte what an std::ostringstream receives (the same
// writer writes into either). A piece the writer hands over is copied in pieces up to the next
// check, so that a check comes every |interval| bytes however large the piece (otherwise one
// 16 MiB string value would be appended whole, with one check after it)
class CheckedStringBuf : public std::streambuf {
  public:
    CheckedStringBuf(std::string *out, DeliveryChecks *checks) : out_(out), checks_(checks) {}

  protected:
    int_type overflow(int_type c) override {
        if (traits_type::eq_int_type(c, traits_type::eof()))
            return traits_type::not_eof(c);
        out_->push_back(traits_type::to_char_type(c));
        checks_->tick();
        return c;
    }
    std::streamsize xsputn(const char_type *s, std::streamsize n) override {
        for (size_t left = static_cast<size_t>(n); left; ) {
            const size_t piece = std::min(left, checks_->room());
            out_->append(s, piece);
            s += piece;
            left -= piece;
            checks_->tick();
        }
        return n;
    }

  private:
    std::string *out_;
    DeliveryChecks *checks_;
};


} // namespace

std::string compress_string(const std::string &str, int compressionlevel, bool gzip,
                            const std::function<void()> &check, size_t max_piece) {
    z_stream zs; // z_stream is zlib's control structure
    memset(&zs, 0, sizeof(zs));

    // zlib's deflate stream by default; a gzip container (window bits + 16) when the
    // client asked for gzip, which is what most HTTP clients send by default
    const int window_bits = 15 + (gzip ? 16 : 0);
    if (deflateInit2(&zs, compressionlevel, Z_DEFLATED, window_bits, 8,
                     Z_DEFAULT_STRATEGY) != Z_OK)
        throw std::runtime_error("deflateInit failed while compressing.");

    // zlib counts its input in 32 bits (uInt avail_in): the text is handed over in pieces of
    // at most that many bytes, the last one with Z_FINISH. Assigned whole, a text of 4 GiB or
    // more would be cut to its size modulo 2^32, and the server would answer 200 with a
    // well-formed stream of that prefix. A text of one piece — every text below 4 GiB — is
    // compressed exactly as by a single call
    const size_t piece_limit = std::numeric_limits<uInt>::max();
    const size_t piece_max = max_piece ? std::min(max_piece, piece_limit) : piece_limit;
    const char *next = str.data();
    size_t left = str.size();

    int ret;
    char outbuffer[32768];
    std::string outstring;

    // retrieve the compressed bytes blockwise
    do {
        if (check) {
            // a stop between blocks (an attempt at its bound, a client gone): the stream is
            // released before the exception leaves
            try {
                check();
            } catch (...) {
                deflateEnd(&zs);
                throw;
            }
        }
        if (!zs.avail_in && left) {
            // the next piece, once zlib consumed the last one (no input is added after
            // Z_FINISH: that flush is used only once nothing is left)
            const size_t piece = std::min(left, piece_max);
            zs.next_in = reinterpret_cast<Bytef *>(const_cast<char *>(next));
            zs.avail_in = static_cast<uInt>(piece);
            next += piece;
            left -= piece;
        }
        zs.next_out = reinterpret_cast<Bytef *>(outbuffer);
        zs.avail_out = sizeof(outbuffer);

        ret = deflate(&zs, left ? Z_NO_FLUSH : Z_FINISH);

        // append the block to the output string (what this call wrote: total_out is a uLong,
        // 32 bits on some platforms)
        outstring.append(outbuffer, sizeof(outbuffer) - zs.avail_out);
    } while (ret == Z_OK);

    deflateEnd(&zs);

    if (ret != Z_STREAM_END) { // an error occurred that was not EOF
        std::ostringstream oss;
        oss << "Exception during zlib compression: (" << ret << ") " << zs.msg;
        throw std::runtime_error(oss.str());
    }

    return outstring;
}
std::string json_text(const Json::Value &value, bool compact, const std::function<void()> &check,
                      double *max_gap_ms) {
    Json::StreamWriterBuilder builder;
    if (compact)
        builder["indentation"] = "";   // the traversal routes: half the bytes of the indented form
    if (!check)
        return Json::writeString(builder, value);
    std::string out;
    DeliveryChecks checks(out, check, kDeliveryCheckBytes, max_gap_ms);
    CheckedStringBuf buf(&out, &checks);
    std::ostream stream(&buf);
    // an exception of the check reaches the caller rather than setting badbit
    stream.exceptions(std::ios::badbit);
    std::unique_ptr<Json::StreamWriter> writer(builder.newStreamWriter());
    writer->write(value, &stream);
    checks.finish();
    return out;
}

std::string assemble_traverse_response(const Json::Value &envelope,
                                       std::vector<std::string> results,
                                       const std::function<void()> &check,
                                       double *max_gap_ms) {
    Json::Value before(Json::objectValue);
    Json::Value after(Json::objectValue);
    for (const std::string &name : envelope.getMemberNames()) {
        if (name == "results")
            throw std::logic_error("assemble_traverse_response: the envelope holds results");
        (name < "results" ? before : after)[name] = envelope[name];
    }
    std::string head = json_text(before, true, check, max_gap_ms);     // {...}
    // the envelope's members after "results" are written first (the usage, the timing: a few
    // KiB), so that the whole response's size is known before a byte of the results is copied
    const std::string tail = json_text(after, true, check, max_gap_ms);  // {...}
    head.pop_back();
    const bool head_members = head.size() > 1;
    const bool tail_members = tail.size() > 2;
    static const std::string kResults = "\"results\":[";
    // the exact size: `{` and the members before, a comma, "results":[, the texts and their
    // commas, `]`, then a comma and the members after with their `}`, or the `}` alone. Reserved
    // at once, the text never grows by doubling, whose last step held an old and a new buffer
    // of up to its size beside the texts
    size_t size = head.size() + head_members + kResults.size() + 1
                + (tail_members ? tail.size() : 1);
    for (size_t i = 0; i < results.size(); ++i) {
        size += results[i].size() + (i > 0);
    }
    std::string out;
    out.reserve(size);
    out += head;
    std::string().swap(head);
    if (head_members)
        out += ',';
    out += kResults;
    // each text is freed once copied: what was copied and what is left to copy are the
    // response once; texts kept alive beside it until the handler returns would stay through
    // the compression and the transport's copy (2.5x the text held with gzip, 3.5x without)
    if (!check) {
        for (size_t i = 0; i < results.size(); ++i) {
            if (i)
                out += ',';
            out += results[i];
            std::string().swap(results[i]);
        }
    } else {
        // each text copied in pieces up to the next check
        DeliveryChecks checks(out, check, kDeliveryCheckBytes, max_gap_ms);
        for (size_t i = 0; i < results.size(); ++i) {
            if (i)
                out += ',';
            const std::string &text = results[i];
            for (size_t at = 0; at < text.size(); ) {
                const size_t piece = std::min(text.size() - at, checks.room());
                out.append(text, at, piece);
                at += piece;
                checks.tick();
            }
            std::string().swap(results[i]);
        }
        checks.finish();
    }
    out += ']';
    if (tail_members) {
        out += ',';
        out.append(tail, 1, std::string::npos);
    } else {
        out += '}';
    }
    // the size computed is the size written: a mismatch would only cost a reallocation, never
    // a byte, so it is asserted where tests run unoptimised
    assert(out.size() == size);
    return out;
}


std::string traversal_error_json(const std::string &msg) {
    Json::Value root;
    root["error"] = msg;
    return Json::writeString(Json::StreamWriterBuilder(), root);
}

RequestAnswer answer_request(const std::string &content, const std::string &encoding,
                             size_t request_id,
                             const std::function<Json::Value(const std::string &)> &process,
                             bool compact,
                             const ResponseControl *control) {
    RequestAnswer answer;
    std::string &ret = answer.body;
    static const std::function<void()> kNoCheck;
    const std::function<void()> &check = control ? control->check : kNoCheck;

    try {
        // Return JSON string
        answer.status = 200;
        if (control && control->write) {
            ret = control->write(process(content), check);
        } else {
            ret = json_text(process(content), compact, check);
        }
        if (!encoding.empty()) {
            Timer compressing;
            const size_t text_bytes = ret.size();
            // in pieces of kCompressPieceBytes, so that a stop (an attempt at its bound, a
            // client gone) is seen at least every 4 MiB of text: the check per 32 KiB block of
            // output alone runs once per tens of MiB on a text that compresses strongly. The
            // compressed bytes are those of one call over the whole text
            // (GraphletServer.CompressionPiecesOfTheServerKeepTheBytes)
            ret = compress_string(ret, control ? control->compression_level : Z_BEST_COMPRESSION,
                                  encoding == "gzip", check, kCompressPieceBytes);
            if (control && control->on_compressed)
                control->on_compressed(text_bytes, compressing.elapsed());
            answer.header.emplace_back("Content-Encoding", encoding);
            answer.header.emplace_back("Content-Length", std::to_string(ret.size()));
        }
        // once more before the response is handed to the transport: a body shorter than the
        // check's interval is otherwise never checked after it was built
        if (check)
            check();
    } catch (const ClientGone &e) {
        // nobody to answer: nothing is written, and the connection is closed rather than
        // kept for a next request
        logger->info("[Server] Request {}: {}; no response written", request_id, e.what());
        return RequestAnswer();
    } catch (const HttpError &e) {
        logger->warn("[Server] Error on request {} ({}): {}", request_id, e.status(), e.what());
        answer.status = e.status();
        answer.header.clear();
        // the body as the route wrote it, uncompressed like every error
        ret = json_text(e.body(), compact);
    } catch (const CurrentlyInitializingError& e) {
        logger->info("[Server] Got a request during initialization. Asked to come back later");
        answer.status = 503;
        answer.header.emplace_back("Retry-After", "60"); // ask to come back in 60 seconds
        ret = traversal_error_json("Server is currently initializing, please come back later.");
    } catch (const std::exception& e) {
        logger->warn("[Server] Error on request {}: {}", request_id, e.what());
        answer.status = 400;
        ret = traversal_error_json(e.what());
    } catch (...) {
        logger->warn("[Server] Error on request {}", request_id);
        answer.status = 500;
        ret = traversal_error_json("Internal server error");
    }
    // a route that answers nobody who left asks once more, whatever was built: its errors
    // too (the success path's check asked already; an error never passed one). Apart from
    // the deadline, so that a 503 at the deadline still reaches a client that is there
    if (control && control->gone && control->gone()) {
        logger->info("[Server] Request {}: the client is gone or the server stops; the {} "
                     "answer is not written", request_id, answer.status);
        return RequestAnswer();
    }
    return answer;
}

} // namespace mtg::cli
