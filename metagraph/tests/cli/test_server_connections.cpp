#include "gtest/gtest.h"

#include <cerrno>
#include <chrono>
#include <cstring>
#include <cstdint>
#include <future>
#include <memory>
#include <optional>
#include <random>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <arpa/inet.h>
#include <netinet/in.h>
#include <poll.h>
#include <sys/socket.h>
#include <unistd.h>

#include <json/json.h>

#include "cli/server_checks.hpp"
#include "cli/server_utils.hpp"


// H3 (train 3; review 5 #1 and #4; SPEC-labeled-traversal-core.md §10.3): the cap on the
// connections the server holds open (CappedHttpServer: at the cap a further connection waits,
// unaccepted, until one closes) and GET /column_labels' text (the unpaged array byte for byte
// as json_text writes it, the pages, the byte cap and its 400). The suites' names are within
// CI's unit filter (*ServerRequest*, *ServerCapacity*)

namespace {

using namespace mtg::cli;
using namespace std::chrono_literals;

std::string compact(const Json::Value &v) {
    return json_text(v, true);
}

Json::Value array_of(const std::vector<std::string> &labels) {
    Json::Value a(Json::arrayValue);
    for (const std::string &l : labels) {
        a.append(l);
    }
    return a;
}

// the labels the list holds from |from| on
std::vector<std::string> visited(const ColumnLabelList &list, uint64_t from) {
    std::vector<std::string> out;
    list.visit_from(from, [&](const std::string &l) {
        out.push_back(l);
        return true;
    });
    return out;
}

// A raw HTTP/1.1 client connection to 127.0.0.1:|port| (blocking, with timeouts by poll)
class RawClient {
  public:
    explicit RawClient(unsigned short port) {
        sockaddr_in addr {};
        addr.sin_family = AF_INET;
        addr.sin_port = htons(port);
        addr.sin_addr.s_addr = htonl(INADDR_LOOPBACK);
        // a transient failure of the host (no local port free for a moment, an interrupted call)
        // is retried, so that a busy machine does not fail the test; the last errno is kept
        for (int attempt = 0; attempt < 50 && !connected_; ++attempt) {
            if (fd_ >= 0)
                ::close(fd_);
            fd_ = ::socket(AF_INET, SOCK_STREAM, 0);
#ifdef SO_NOSIGPIPE
            int one = 1;
            ::setsockopt(fd_, SOL_SOCKET, SO_NOSIGPIPE, &one, sizeof(one));
#endif
            connected_ = fd_ >= 0
                    && ::connect(fd_, reinterpret_cast<sockaddr *>(&addr), sizeof(addr)) == 0;
            if (connected_)
                break;
            error_ = errno;
            if (error_ != EADDRNOTAVAIL && error_ != EAGAIN && error_ != EINTR
                    && error_ != EMFILE && error_ != ENFILE && error_ != ENOBUFS)
                break;
            std::this_thread::sleep_for(std::chrono::milliseconds(50));
        }
    }
    // why it did not connect (errno), for the test's message
    std::string error() const { return connected_ ? std::string() : std::strerror(error_); }
    ~RawClient() { close(); }
    RawClient(const RawClient &) = delete;
    RawClient& operator=(const RawClient &) = delete;

    bool connected() const { return connected_; }
    void close() {
        if (fd_ >= 0)
            ::close(fd_);
        fd_ = -1;
    }
    // closes with a reset (SO_LINGER 0)
    void reset() {
        linger l { 1, 0 };
        ::setsockopt(fd_, SOL_SOCKET, SO_LINGER, &l, sizeof(l));
        close();
    }
    bool send_raw(const std::string &bytes) {
#ifdef MSG_NOSIGNAL
        const int flags = MSG_NOSIGNAL;
#else
        const int flags = 0;
#endif
        return ::send(fd_, bytes.data(), bytes.size(), flags) == static_cast<ssize_t>(bytes.size());
    }
    // whether the server closed the connection within |wait| (the end of the stream, or a reset),
    // whatever it sent before
    bool closed_by_peer(std::chrono::milliseconds wait) {
        const auto until = std::chrono::steady_clock::now() + wait;
        while (std::chrono::steady_clock::now() < until) {
            pollfd p { fd_, POLLIN, 0 };
            if (::poll(&p, 1, 50) <= 0)
                continue;
            char buf[4096];
            const ssize_t n = ::recv(fd_, buf, sizeof(buf), 0);
            if (n <= 0)
                return true;
        }
        return false;
    }
    bool send_get(const std::string &path) {
        const std::string req = "GET " + path + " HTTP/1.1\r\nHost: x\r\n\r\n";
#ifdef MSG_NOSIGNAL
        const int flags = MSG_NOSIGNAL;
#else
        const int flags = 0;
#endif
        return ::send(fd_, req.data(), req.size(), flags) == static_cast<ssize_t>(req.size());
    }
    // what arrives within |wait| (until the response's head and its Content-Length body)
    std::string read_answer(std::chrono::milliseconds wait) {
        std::string got;
        const auto until = std::chrono::steady_clock::now() + wait;
        while (std::chrono::steady_clock::now() < until) {
            pollfd p { fd_, POLLIN, 0 };
            const auto left = std::chrono::duration_cast<std::chrono::milliseconds>(
                    until - std::chrono::steady_clock::now()).count();
            if (::poll(&p, 1, static_cast<int>(std::max<int64_t>(left, 1))) <= 0)
                continue;
            char buf[4096];
            const ssize_t n = ::recv(fd_, buf, sizeof(buf), 0);
            if (n <= 0)
                break;
            got.append(buf, n);
            const size_t head = got.find("\r\n\r\n");
            if (head != std::string::npos) {
                const size_t cl = got.find("Content-Length: ");
                if (cl != std::string::npos && cl < head) {
                    const size_t len = std::stoul(got.substr(cl + 16));
                    if (got.size() >= head + 4 + len)
                        break;
                }
            }
        }
        return got;
    }

  private:
    int fd_ = -1;
    bool connected_ = false;
    int error_ = 0;
};

// A CappedHttpServer on 127.0.0.1 (an ephemeral port), with GET /ping answering "pong", run on
// |threads| handler threads in its own thread
class TestServer {
  public:
    // |timeout_request|: the header timeout (s); |max_body_mb|: the body cap, answered 413 as
    // the server answers it (refuse_body_too_large), 0: none
    TestServer(size_t max_connections, size_t threads = 2, long timeout_request = 30,
               uint64_t max_body_mb = 0) {
        server_.config.address = "127.0.0.1";
        server_.config.port = 0;
        server_.config.thread_pool_size = threads;
        server_.config.timeout_request = timeout_request;
        if (max_body_mb) {
            server_.config.max_request_streambuf_size = static_cast<size_t>(max_body_mb) << 20;
            server_.on_error = [max_body_mb](std::shared_ptr<HttpServer::Request> request,
                                             const SimpleWeb::error_code &ec) {
                if (ec == SimpleWeb::make_error_code::make_error_code(
                            SimpleWeb::errc::message_size)) {
                    refuse_body_too_large(request, max_body_mb);
                }
            };
        }
        server_.set_max_connections(max_connections);
        server_.resource["^/ping$"]["GET"] = [](std::shared_ptr<HttpServer::Response> response,
                                               std::shared_ptr<HttpServer::Request>) {
            response->write("pong");
        };
        auto port = std::make_shared<std::promise<unsigned short>>();
        auto listening = port->get_future();
        thread_ = std::thread([this, port]() {
            server_.start([port](unsigned short p) { port->set_value(p); });
        });
        port_ = listening.get();
    }
    ~TestServer() {
        // as the server's shutdown: the workers stopped first, the connections closed after
        // they returned (run_server)
        server_.io_service->stop();
        if (thread_.joinable())
            thread_.join();
        server_.stop();
    }
    unsigned short port() const { return port_; }
    CappedHttpServer& server() { return server_; }

  private:
    CappedHttpServer server_;
    std::thread thread_;
    unsigned short port_ = 0;
};

bool answered(const std::string &got) {
    return got.rfind("HTTP/1.1 200", 0) == 0 && got.find("pong") != std::string::npos;
}

// waits up to |wait| for |pred|
template <class Pred>
bool eventually(Pred pred, std::chrono::milliseconds wait = 5000ms) {
    const auto until = std::chrono::steady_clock::now() + wait;
    while (std::chrono::steady_clock::now() < until) {
        if (pred())
            return true;
        std::this_thread::sleep_for(5ms);
    }
    return pred();
}

} // namespace


// At the cap the server stops accepting: a further connection is completed by the kernel and
// its request waits unread (no answer, no refusal) until a connection closes; then it is
// accepted and answered. A request on a connection already open is not delayed
TEST(ServerCapacityConnections, AtTheCapAFurtherConnectionWaitsUntilOneCloses) {
    TestServer t(/* max_connections */ 3, /* threads */ 2);
    RawClient a(t.port()), b(t.port()), c(t.port());
    ASSERT_TRUE(a.connected() && b.connected() && c.connected())
            << a.error() << " " << b.error() << " " << c.error();
    // a, b and c are open (kept alive after their answers)
    for (RawClient *x : { &a, &b, &c }) {
        ASSERT_TRUE(x->send_get("/ping"));
        EXPECT_TRUE(answered(x->read_answer(5000ms)));
    }
    EXPECT_TRUE(eventually([&]() { return t.server().pauses() >= 1; }));
    EXPECT_EQ(3u, t.server().open_connections());
    // the fourth connects (the kernel's listen queue) but is not answered
    RawClient d(t.port());
    ASSERT_TRUE(d.connected()) << d.error();
    ASSERT_TRUE(d.send_get("/ping"));
    EXPECT_EQ("", d.read_answer(400ms));
    // a request on an open connection is answered at once
    ASSERT_TRUE(b.send_get("/ping"));
    EXPECT_TRUE(answered(b.read_answer(5000ms)));
    // one closes: the fourth is accepted and answered
    a.close();
    EXPECT_TRUE(answered(d.read_answer(5000ms)));
    EXPECT_TRUE(eventually([&]() { return t.server().open_connections() == 3; }));
}

// Without a cap (0) the accept loop is the library's: every connection is answered at once
TEST(ServerCapacityConnections, NoCapAcceptsEveryConnection) {
    TestServer t(/* max_connections */ 0);
    std::vector<std::unique_ptr<RawClient>> clients;
    for (int i = 0; i < 20; ++i) {
        clients.push_back(std::make_unique<RawClient>(t.port()));
        ASSERT_TRUE(clients.back()->connected()) << clients.back()->error();
        ASSERT_TRUE(clients.back()->send_get("/ping"));
        EXPECT_TRUE(answered(clients.back()->read_answer(5000ms))) << i;
    }
    EXPECT_EQ(0u, t.server().pauses());
    EXPECT_EQ(0u, t.server().max_connections());
    clients.clear();
    // only the connection awaiting an accept is left
    EXPECT_TRUE(eventually([&]() { return t.server().open_connections() == 1; }));
}

// A server stopped while paused at its cap stops cleanly, and the closes that follow resume
// nothing
TEST(ServerCapacityConnections, StopWhilePausedAtTheCap) {
    auto t = std::make_unique<TestServer>(/* max_connections */ 2, /* threads */ 1);
    RawClient a(t->port()), b(t->port());
    ASSERT_TRUE(a.connected() && b.connected()) << a.error() << " " << b.error();
    for (RawClient *x : { &a, &b }) {
        ASSERT_TRUE(x->send_get("/ping"));
        EXPECT_TRUE(answered(x->read_answer(5000ms)));
    }
    EXPECT_TRUE(eventually([&]() { return t->server().pauses() >= 1; }));
    const auto start = std::chrono::steady_clock::now();
    t.reset();
    // the stop returned at once, and closed both connections (the end of their streams)
    EXPECT_LT(std::chrono::steady_clock::now() - start, 5s);
    EXPECT_TRUE(a.closed_by_peer(2000ms));
    EXPECT_TRUE(b.closed_by_peer(2000ms));
    a.close();
    b.close();
}

// A connection that ends on an error path gives its place back as a normal close does: a reset
// by the client, the header timeout, the 413 of a body over the cap (its drain ends when the
// client closes); after each the count is the pending accept alone, and the cap holds its N
// connections again (the review of H3, N4)
TEST(ServerCapacityConnections, ErrorPathsGiveTheirPlaceBack) {
    // a header timeout of 3 s: well above what the last part takes with its connections idle
    TestServer t(/* max_connections */ 3, /* threads */ 2, /* timeout_request */ 3,
                 /* max_body_mb */ 1);
    auto settled = [&]() {
        return eventually([&]() { return t.server().open_connections() == 1; });
    };
    // a reset after a request was answered on the connection
    {
        RawClient r(t.port());
        ASSERT_TRUE(r.connected()) << r.error();
        ASSERT_TRUE(r.send_get("/ping"));
        EXPECT_TRUE(answered(r.read_answer(5000ms)));
        r.reset();
    }
    EXPECT_TRUE(settled());
    // a header that never ends: the server closes it at the header timeout (3 s)
    {
        RawClient h(t.port());
        ASSERT_TRUE(h.connected()) << h.error();
        ASSERT_TRUE(h.send_raw("GET /ping HTTP/1.1\r\nHost: x\r\n"));
        EXPECT_TRUE(h.closed_by_peer(10000ms));
    }
    EXPECT_TRUE(settled());
    // a body over the cap: the 413, the drain, the client's close
    {
        RawClient d(t.port());
        ASSERT_TRUE(d.connected()) << d.error();
        ASSERT_TRUE(d.send_raw("POST /ping HTTP/1.1\r\nHost: x\r\nContent-Length: 2097152\r\n\r\n"));
        const std::string got = d.read_answer(5000ms);
        EXPECT_EQ(0u, got.rfind("HTTP/1.1 413", 0)) << got.substr(0, 200);
        d.close();
    }
    EXPECT_TRUE(settled());
    // the cap holds its N again: three kept-alive connections answered, the fourth waits
    RawClient a(t.port()), b(t.port()), c(t.port());
    for (RawClient *x : { &a, &b, &c }) {
        ASSERT_TRUE(x->send_get("/ping"));
        EXPECT_TRUE(answered(x->read_answer(5000ms)));
    }
    RawClient w(t.port());
    ASSERT_TRUE(w.send_get("/ping"));
    EXPECT_EQ("", w.read_answer(400ms));
    a.close();
    EXPECT_TRUE(answered(w.read_answer(5000ms)));
}


// The unpaged list is json_text's styled array byte for byte, whatever the labels hold
TEST(ServerRequestColumnLabels, TheArrayIsJsonTextByteForByte) {
    const std::vector<std::vector<std::string>> cases = {
        {},
        { "a" },
        { "ENST00000456328", "taxid:9606", "" },
        { "quote\"d", "back\\slash", "tab\tnew\nline", std::string("nul\0x", 5), "\x7f" "del",
          "\x01\x1f", "caf\xc3\xa9", "\xf0\x9f\xa7\xac", "bad\xff" "utf8", "/a/b.fa.gz" },
    };
    for (const auto &labels : cases) {
        const ColumnLabelList list({ &labels });
        const std::optional<std::string> got = column_labels_text(list, 0);
        ASSERT_TRUE(got);
        EXPECT_EQ(json_text(array_of(labels), false), *got);
        // split over parts (one of them empty), the same list
        if (labels.size() >= 2) {
            const std::vector<std::string> first(labels.begin(), labels.begin() + 1);
            const std::vector<std::string> none;
            const std::vector<std::string> rest(labels.begin() + 1, labels.end());
            const ColumnLabelList split({ &first, &none, &rest });
            EXPECT_EQ(*got, column_labels_text(split, 0).value());
        }
    }
    std::vector<std::string> many;
    for (int i = 0; i < 1000; ++i) {
        many.push_back("label-" + std::to_string(i));
    }
    const ColumnLabelList list({ &many });
    EXPECT_EQ(json_text(array_of(many), false), column_labels_text(list, 0).value());
}

// The byte cap of the unpaged list: the labels' JSON strings (quotes and escapes counted) sum to
// at most the cap, else nothing; the refusal's body and its text (error first)
TEST(ServerRequestColumnLabels, TheUnpagedListStopsAtTheCap) {
    const std::vector<std::string> labels = { "ab", "c\"d", "e" };   // 4 + 6 + 3 = 13 bytes
    const ColumnLabelList list({ &labels });
    EXPECT_TRUE(column_labels_text(list, 13));
    EXPECT_FALSE(column_labels_text(list, 12));
    EXPECT_FALSE(column_labels_text(list, 1));
    EXPECT_TRUE(column_labels_text(list, 0));   // unlimited
    const Json::Value body = column_labels_too_large_json(651, 16777216);
    EXPECT_EQ(R"x({"error":"the column labels (651 labels) take more than 16777216 bytes of )x"
              R"x(JSON strings, more than one answer holds (--max-column-labels-bytes): page )x"
              R"x(through them with GET /column_labels?offset=0&limit=<n>, following )x"
              R"x(next_offset","code":"column_labels_too_large","max_column_labels_bytes":)x"
              R"x(16777216,"total":651})x",
              error_first_json_text(body));
    // the capabilities state exactly the body's members, and the cap and the total
    const Json::Value caps = column_labels_capabilities_json(16777216, 651);
    Json::Value members(Json::arrayValue);
    for (const std::string &m : body.getMemberNames()) {
        members.append(m);
    }
    EXPECT_EQ(compact(members), compact(caps["too_large"]["members"]));
    EXPECT_EQ(400, caps["too_large"]["status"].asInt());
    EXPECT_EQ("column_labels_too_large", caps["too_large"]["code"].asString());
    EXPECT_EQ(R"x({"max_column_labels_bytes":16777216,"paging":["limit","offset"],"too_large":)x"
              R"x({"code":"column_labels_too_large","members":["code","error",)x"
              R"x("max_column_labels_bytes","total"],"status":400},"total":651})x",
              compact(caps));
    const Json::Value unlimited = column_labels_capabilities_json(0, std::nullopt);
    EXPECT_TRUE(unlimited["max_column_labels_bytes"].isNull());
    EXPECT_TRUE(unlimited["total"].isNull());
    EXPECT_TRUE(caps["total"].isUInt64());
}

// The pages: compact, members in byte order (json_text's), at most |limit| labels and as many as
// fit the cap, at least one, next_offset null at the end, an offset past it empty
TEST(ServerRequestColumnLabels, Pages) {
    const std::vector<std::string> p1 = { "l0", "l1", "l2" };
    const std::vector<std::string> p2 = {};
    const std::vector<std::string> p3 = { "l3", "l\xc3\xa9", "l5" };
    const ColumnLabelList list({ &p1, &p2, &p3 });
    ASSERT_EQ(6u, list.total());
    EXPECT_EQ(std::vector<std::string>({ "l2", "l3", "l\xc3\xa9", "l5" }), visited(list, 2));
    EXPECT_EQ(std::vector<std::string>({ "l5" }), visited(list, 5));
    EXPECT_TRUE(visited(list, 6).empty());

    auto page = [&](uint64_t offset, std::optional<uint64_t> limit, uint64_t max_bytes) {
        ColumnLabelsPage p;
        p.paged = true;
        p.offset = offset;
        p.limit = limit;
        return column_labels_page_text(list, p, max_bytes);
    };
    auto expected = [](std::vector<std::string> labels, Json::Value next, uint64_t offset) {
        Json::Value v(Json::objectValue);
        v["labels"] = array_of(labels);
        v["next_offset"] = next;
        v["offset"] = Json::UInt64(offset);
        v["total"] = Json::UInt64(6);
        return json_text(v, true);
    };
    EXPECT_EQ(R"x({"labels":["l1","l2"],"next_offset":3,"offset":1,"total":6})x",
              page(1, 2, 0));
    EXPECT_EQ(expected({ "l0", "l1", "l2", "l3", "l\xc3\xa9", "l5" }, Json::Value(), 0),
              page(0, std::nullopt, 0));
    // the non-ASCII label as json_text escapes it
    EXPECT_EQ(expected({ "l3", "l\xc3\xa9" }, Json::UInt64(5), 3), page(3, 2, 0));
    // the cap: "l0" "l1" are 4 bytes each; 9 bytes hold two
    EXPECT_EQ(expected({ "l0", "l1" }, Json::UInt64(2), 0), page(0, std::nullopt, 9));
    // at least one label, however small the cap
    EXPECT_EQ(expected({ "l0" }, Json::UInt64(1), 0), page(0, 5, 1));
    // the last page and past the end
    EXPECT_EQ(expected({ "l5" }, Json::Value(), 5), page(5, 10, 0));
    EXPECT_EQ(expected({}, Json::Value(), 6), page(6, 10, 0));
    EXPECT_EQ(expected({}, Json::Value(), 1000), page(1000, std::nullopt, 0));
    // walking the pages gives the list
    std::vector<std::string> all;
    uint64_t offset = 0;
    for (int guard = 0; guard < 100; ++guard) {
        const std::string t = page(offset, 1, 0);
        Json::Value v;
        std::unique_ptr<Json::CharReader> reader(Json::CharReaderBuilder().newCharReader());
        ASSERT_TRUE(reader->parse(t.data(), t.data() + t.size(), &v, nullptr));
        for (const Json::Value &l : v["labels"]) {
            all.push_back(l.asString());
        }
        if (v["next_offset"].isNull())
            break;
        offset = v["next_offset"].asUInt64();
    }
    EXPECT_EQ(visited(list, 0), all);
    // an empty list
    const ColumnLabelList empty;
    ColumnLabelsPage p;
    p.paged = true;
    EXPECT_EQ(R"x({"labels":[],"next_offset":null,"offset":0,"total":0})x",
              column_labels_page_text(empty, p, 0));
    EXPECT_EQ("[]", column_labels_text(empty, 0).value());
}

// The parameters: none (unpaged), offset and limit (paged), each refused when unknown,
// repeated or not an integer in its range
// GPT review 6 #13: a label's JSON string is written in place and measured without being built,
// both byte for byte jsoncpp's (Json::writeString, as json_text writes it): escapes, control
// characters, every non-ASCII code point as \uxxxx, and whatever jsoncpp makes of bytes that are
// no UTF-8 (truncated, overlong, surrogates, beyond U+10FFFF, stray continuation bytes)
TEST(ServerRequestColumnLabels, TheJsonStringIsJsoncppsByteForByte) {
    auto jsoncpp = [](const std::string &label) {
        return Json::writeString(Json::StreamWriterBuilder(), Json::Value(label));
    };
    auto check = [&](const std::string &label) {
        std::string out = "x";
        const size_t n = append_json_string(&out, label);
        const std::string want = jsoncpp(label);
        ASSERT_EQ("x" + want, out) << ::testing::PrintToString(label);
        ASSERT_EQ(want.size(), n);
        ASSERT_EQ(want.size(), json_string_length(label)) << ::testing::PrintToString(label);
    };
    for (const std::string &label : {
             std::string(), std::string("a"), std::string("quote\"d back\\slash"),
             std::string("\b\f\n\r\t\x01\x1f\x7f"), std::string("nul\0x", 5),
             std::string("caf\xc3\xa9"), std::string("\xe2\x82\xac"), std::string("\xf0\x9f\xa7\xac"),
             std::string("\xc3"), std::string("\xe2\x82"), std::string("\xf0\x9f\xa7"),
             std::string("\xc0\x80"), std::string("\xe0\x80\x80"), std::string("\xf0\x80\x80\x80"),
             std::string("\xed\xa0\x80"), std::string("\xf4\x90\x80\x80"), std::string("\xf8\x88"),
             std::string("\xff\xfe"), std::string("\x80" "abc"), std::string("a\xc3(b"),
             std::string("\xef\xbf\xbf"), std::string("\xf4\x8f\xbf\xbf") }) {
        check(label);
    }
    // random bytes, and random mixes of the pieces above
    std::mt19937 rng(20261109);
    const std::vector<std::string> pieces = { "a", "\"", "\\", "\n", "\x01", "\x7f", "\xc3\xa9",
                                              "\xe2\x82\xac", "\xf0\x9f\xa7\xac", "\xc3", "\xe2",
                                              "\xf0", "\x80", "\xff", "\xed\xa0\x80", "/" };
    for (int i = 0; i < 20000; ++i) {
        std::string label;
        const int n = static_cast<int>(rng() % 12);
        for (int j = 0; j < n; ++j) {
            if (i % 2) {
                label.push_back(static_cast<char>(rng() % 256));
            } else {
                label += pieces[rng() % pieces.size()];
            }
        }
        check(label);
    }
}

// GPT review 6 #13: a label that does not fit is measured, never built. A page stops before it
// (the answer's buffer stays the page's size: before, the label was appended, 6 bytes per control
// byte, and cut off again, its capacity kept), the unpaged list is refused before anything is
// written; a label that fits exactly is answered, one byte more is not; a first label longer
// than the cap is a page alone
TEST(ServerRequestColumnLabels, ALabelThatDoesNotFitIsNeverBuilt) {
    constexpr uint64_t kCap = 1 << 20;
    const std::vector<std::string> labels = { "a", std::string(kCap / 2, '\x01') };
    ASSERT_EQ(6 * (kCap / 2) + 2, json_string_length(labels[1]));
    const ColumnLabelList list({ &labels });
    ColumnLabelsPage p;
    p.paged = true;
    const std::string first = column_labels_page_text(list, p, kCap);
    EXPECT_EQ(R"x({"labels":["a"],"next_offset":1,"offset":0,"total":2})x", first);
    EXPECT_LT(first.capacity(), 4096u);
    EXPECT_FALSE(column_labels_text(list, kCap));
    // alone over the cap: the page of it alone, whole
    p.offset = 1;
    const std::string alone = column_labels_page_text(list, p, kCap);
    EXPECT_EQ(6 * (kCap / 2) + 2, alone.size() - std::string(
            R"x({"labels":[],"next_offset":null,"offset":1,"total":2})x").size());
    // an index label near the cap: kCap - 2 characters fit exactly (quotes counted), one more not
    for (uint64_t size : { kCap - 3, kCap - 2, kCap - 1 }) {
        const std::vector<std::string> near = { std::string(size, 'x') };
        const ColumnLabelList one({ &near });
        const std::optional<std::string> text = column_labels_text(one, kCap);
        EXPECT_EQ(size + 2 <= kCap, text.has_value()) << size;
        if (text) {
            EXPECT_EQ(json_text(array_of(near), false), *text);
            // built at its final size
            EXPECT_LE(text->capacity(), text->size() + 32);
        }
        const std::vector<std::string> two = { "b", std::string(size, 'x') };
        const ColumnLabelList pair({ &two });
        ColumnLabelsPage q;
        q.paged = true;
        // "b" takes 3 bytes: the long label fits beside it only when 3 + size + 2 <= kCap
        const std::string page = column_labels_page_text(pair, q, kCap);
        if (3 + size + 2 <= kCap) {
            EXPECT_EQ(R"x({"labels":["b",")x" + two[1] + R"x("],"next_offset":null,"offset":0,)x"
                      R"x("total":2})x", page);
        } else {
            EXPECT_EQ(R"x({"labels":["b"],"next_offset":1,"offset":0,"total":2})x", page) << size;
        }
    }
}

TEST(ServerRequestColumnLabels, Parameters) {
    using Query = std::vector<std::pair<std::string, std::string>>;
    EXPECT_FALSE(column_labels_parameters({}).paged);
    ColumnLabelsPage p = column_labels_parameters(Query { { "offset", "7" } });
    EXPECT_TRUE(p.paged);
    EXPECT_EQ(7u, p.offset);
    EXPECT_FALSE(p.limit);
    p = column_labels_parameters(Query { { "limit", "3" } });
    EXPECT_TRUE(p.paged);
    EXPECT_EQ(0u, p.offset);
    EXPECT_EQ(3u, p.limit.value());
    p = column_labels_parameters(Query { { "limit", "9007199254740991" },
                                         { "offset", "9007199254740991" } });
    EXPECT_EQ(9007199254740991u, p.offset);
    auto refused = [](const Query &q) {
        try {
            column_labels_parameters(q);
        } catch (const std::invalid_argument &e) {
            return std::string(e.what());
        }
        return std::string("accepted");
    };
    EXPECT_EQ("Bad request: unknown parameter 'graph' of GET /column_labels (offset, limit)",
              refused({ { "graph", "x" } }));
    EXPECT_EQ("Bad request: the parameter 'offset' is given twice",
              refused({ { "offset", "1" }, { "offset", "2" } }));
    EXPECT_EQ("Bad request: GET /column_labels: 'limit' must be an integer in [1, "
              "9007199254740991], got '0'", refused({ { "limit", "0" } }));
    // a client's key or value is quoted clipped past 1 KiB, as every error text quotes it (H1)
    const std::string long_key(5000, 'k');
    EXPECT_EQ("Bad request: unknown parameter '" + std::string(1024, 'k') + "... (5000 bytes)' of "
              "GET /column_labels (offset, limit)", refused({ { long_key, "1" } }));
    EXPECT_EQ("Bad request: GET /column_labels: 'limit' must be an integer in [1, "
              "9007199254740991], got '" + std::string(1024, '9') + "... (3000 bytes)'",
              refused({ { "limit", std::string(3000, '9') } }));
    for (const char *bad : { "", "-1", "+1", "1.0", "0x10", " 1", "9007199254740992",
                             "99999999999999999999" }) {
        EXPECT_EQ("Bad request: GET /column_labels: 'offset' must be an integer in [0, "
                  "9007199254740991], got '" + std::string(bad) + "'",
                  refused({ { "offset", bad } })) << bad;
    }
}

// request_intake states the connection cap (null without one)
TEST(ServerRequestIntake, CapabilitiesStateTheConnectionCap) {
    EXPECT_EQ("256", compact(request_intake_json(2, 262144, 256)["max_connections"]));
    EXPECT_TRUE(request_intake_json(2, 262144, 0)["max_connections"].isNull());
}
