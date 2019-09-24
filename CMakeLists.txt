cmake_minimum_required(VERSION 3.14)
project(parall)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)
include_directories(assignments)
include_directories(assignments/first)

add_executable(parall
        assignments/first/first.cpp
        assignments/first/random.h
        assignments/first/main.cpp)
