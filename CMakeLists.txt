cmake_minimum_required(VERSION 3.10)

project(LLL)
set(CMAKE_CXX_STANDARD 14)

#include("${CMAKE_BINARY_DIR}/conanbuildinfo.cmake")
#conan_basic_setup(NO_OUTPUT_DIRS)

include_directories(.)

add_executable(LLL
        lll.cpp
        Matrix.h
        Vector.h Manger.h)
