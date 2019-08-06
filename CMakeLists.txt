cmake_minimum_required(VERSION 3.10)

project(LLL)
set(CMAKE_CXX_STANDARD 14)

include(ExternalProject)
set(EXTERNAL_INSTALL_LOCATION ${CMAKE_BINARY_DIR}/external)

include_directories(${EXTERNAL_INSTALL_LOCATION}/include)
link_directories(${EXTERNAL_INSTALL_LOCATION}/lib)

ExternalProject_Add(googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG master
        SOURCE_DIR "${CMAKE_CURRENT_BINARY_DIR}/googletest-src"
        BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/googletest-build"
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION} -Dgtest_build_tests=OFF
        UPDATE_DISCONNECTED true
        )


include_directories(.)

add_library(LLLlib SHARED
        lll.cpp
        Matrix.h
        Vector.h
        Manger.h)

add_executable(LLL
        lll.cpp
        )

add_dependencies(LLL LLLlib)
target_link_libraries(LLL LLLlib)

add_executable(vector_tests vector_tests.cpp)
add_dependencies(vector_tests googletest LLLlib)
target_link_libraries(vector_tests gtest gtest_main LLLlib)