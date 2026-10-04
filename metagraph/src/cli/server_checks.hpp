#ifndef __METAGRAPH_SERVER_CHECKS_HPP__
#define __METAGRAPH_SERVER_CHECKS_HPP__

#include <functional>
#include <string>

#include <json/json.h>


namespace mtg {
namespace cli {

// The parts of the server's request handling that need no HTTP library (server_utils.cpp),
// declared apart so that they are tested on their own

// Whether the peer of the connected TCP socket |fd| is gone: a non-blocking peek (consuming
// nothing) finds the connection closed — an orderly close or a half-close (the peer will send
// nothing more; nginx's default treats it the same) — or reset, or |fd| is no socket. Data
// waiting (a pipelined request, a trailing CRLF) hides a close behind it from the peek, so then
// the kernel's TCP state is read (Linux TCP_INFO, macOS TCP_CONNECTION_INFO): CLOSE_WAIT (the
// FIN arrived) or CLOSED (reset) is gone, otherwise connected. On another platform, or for a
// socket that is not TCP, data waiting reads "connected" (a close behind it is then not seen
// until the data is read). Nothing to read yet: connected. An error that says nothing about
// the peer reads "connected".
bool peer_closed(int fd);

// |value| as the server writes it: compact (no indentation) or the default writer's
// indentation, byte for byte what Json::writeString writes; with |check|, called every 64 KiB
// of text, whose exception stops the writing and reaches the caller
std::string json_text(const Json::Value &value, bool compact,
                      const std::function<void()> &check = nullptr);

} // namespace cli
} // namespace mtg

#endif // __METAGRAPH_SERVER_CHECKS_HPP__
