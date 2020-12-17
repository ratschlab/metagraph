#pragma once

#if _USE_FOLLY
#define FOLLY_HAVE_LIBGLOG 0
#define FOLLY_USE_JEMALLOC 1
#endif
