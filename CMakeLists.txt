cmake_minimum_required(VERSION 2.8.12)
project(peng_motif CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(FindASan)

set (PENG_VERSION_MAJOR 1)
set (PENG_VERSION_MINOR 0)
set (PENG_VERSION_PATCH 0)
set (PENG_DATE "04-09-2017")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")

option(test "Build all tests." OFF)

if (test)
    add_subdirectory(lib/gtest-1.8.0)
    include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
    enable_testing()

    add_subdirectory(test)
endif()


INSTALL(FILES
  LICENSE
  README.md
  DESTINATION share/PEnG-motif
)

INSTALL(DIRECTORY
  scripts/
  DESTINATION bin
  PATTERN *
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE
)

include_directories(src)
add_subdirectory(src)
