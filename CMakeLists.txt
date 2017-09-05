cmake_minimum_required(VERSION 2.8.12)
project(peng_motif CXX)

set (PENG_VERSION_MAJOR 1)
set (PENG_VERSION_MINOR 0)
set (PENG_VERSION_PATCH 0)
set (PENG_DATE "04-09-2017")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall")

INSTALL(FILES
  LICENSE
  README.md
  DESTINATION share/PEnG-motif
)

INSTALL(DIRECTORY
  scripts/
  DESTINATION bin
)

include_directories(src)
add_subdirectory(src)
