cmake_minimum_required(VERSION 3.18)
project(dbscan)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(CMAKE_CXX_STANDARD 20)

find_package(Boost COMPONENTS program_options filesystem REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})
add_executable(dbscan main.cpp DataPoint.hpp ClusterPoint.hpp Cache.hpp)
target_link_libraries(dbscan ${Boost_LIBRARIES})
