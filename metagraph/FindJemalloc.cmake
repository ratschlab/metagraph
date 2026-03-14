# from https://github.com/COMBINE-lab/salmon/blob/master/cmake/Modules/FindJemalloc.cmake
# From: https://raw.githubusercontent.com/STEllAR-GROUP/hpx/master/cmake/FindJemalloc.cmake
# Copyright (c)      2014 Thomas Heller
# Copyright (c) 2007-2012 Hartmut Kaiser
# Copyright (c) 2010-2011 Matt Anderson
# Copyright (c) 2011      Bryce Lelbach
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

find_package(PkgConfig)
pkg_check_modules(PC_JEMALLOC QUIET libjemalloc)

# Try to use jemalloc-config if available
find_program(JEMALLOC_CONFIG jemalloc-config)

if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "arm64")
  set(HOMEBREW_DIR /opt/homebrew)
elseif(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
  set(HOMEBREW_DIR /usr/local)
else()
  set(HOMEBREW_DIR ~/.linuxbrew)
endif()

unset(JEMALLOC_INCLUDE_DIR CACHE)
find_path(JEMALLOC_INCLUDE_DIR jemalloc/jemalloc.h
  HINTS
    ${JEMALLOC_ROOT} ENV JEMALLOC_ROOT
    ${PC_JEMALLOC_MINIMAL_INCLUDEDIR}
    ${PC_JEMALLOC_MINIMAL_INCLUDE_DIRS}
    ${PC_JEMALLOC_INCLUDEDIR}
    ${PC_JEMALLOC_INCLUDE_DIRS}
    ${HOMEBREW_DIR}/include
  PATH_SUFFIXES include)

unset(JEMALLOC_LIBRARY CACHE)
find_library(JEMALLOC_LIBRARY NAMES jemalloc libjemalloc
  HINTS
    ${JEMALLOC_ROOT} ENV JEMALLOC_ROOT
    ${PC_JEMALLOC_MINIMAL_LIBDIR}
    ${PC_JEMALLOC_MINIMAL_LIBRARY_DIRS}
    ${PC_JEMALLOC_LIBDIR}
    ${PC_JEMALLOC_LIBRARY_DIRS}
    ${HOMEBREW_DIR}/lib
  PATH_SUFFIXES lib lib64)

# If jemalloc-config is available, use it as additional hint
if(JEMALLOC_CONFIG)
  execute_process(COMMAND ${JEMALLOC_CONFIG} --prefix
    OUTPUT_VARIABLE JEMALLOC_PREFIX
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  # Add jemalloc-config prefix to search paths if not already found
  if(NOT JEMALLOC_INCLUDE_DIR)
    find_path(JEMALLOC_INCLUDE_DIR jemalloc/jemalloc.h
      HINTS ${JEMALLOC_PREFIX}
      PATH_SUFFIXES include)
  endif()

  if(NOT JEMALLOC_LIBRARY)
    find_library(JEMALLOC_LIBRARY NAMES jemalloc libjemalloc
      HINTS ${JEMALLOC_PREFIX}
      PATH_SUFFIXES lib lib64)
  endif()
endif()

if(JEMALLOC_INCLUDE_DIR)
  set(_version_regex "^#define[ \t]+JEMALLOC_VERSION[ \t]+\"([^\"]+)\".*")
  file(STRINGS "${JEMALLOC_INCLUDE_DIR}/jemalloc/jemalloc.h"
    JEMALLOC_VERSION REGEX "${_version_regex}")
  string(REGEX REPLACE "${_version_regex}" "\\1"
    JEMALLOC_VERSION "${JEMALLOC_VERSION}")
  unset(_version_regex)
endif()

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set JEMALLOC_FOUND to TRUE
# if all listed variables are TRUE and the requested version matches.
find_package_handle_standard_args(Jemalloc REQUIRED_VARS
                                  JEMALLOC_LIBRARY JEMALLOC_INCLUDE_DIR
                                  VERSION_VAR JEMALLOC_VERSION)


if(JEMALLOC_FOUND)
  set(JEMALLOC_LIBRARIES    ${JEMALLOC_LIBRARY})
  set(JEMALLOC_INCLUDE_DIRS ${JEMALLOC_INCLUDE_DIR})
endif()

mark_as_advanced(JEMALLOC_INCLUDE_DIR JEMALLOC_LIBRARY)
