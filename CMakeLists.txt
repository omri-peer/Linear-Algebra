cmake_minimum_required(VERSION 3.10)

project(attacker)
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

add_library(LinearAlgebra SHARED
        lll.h
        Matrix.h
        Vector.h
        Manger.h
        Main.cpp)

add_executable(attacker
        Main.cpp
        )

add_dependencies(attacker LinearAlgebra)
target_link_libraries(attacker LinearAlgebra)

add_executable(tests vector_tests.cpp matrix_tests.cpp lll_tests.cpp Manger_tests.cpp Vector.cpp lll.cpp)
add_dependencies(tests googletest LinearAlgebra)
target_link_libraries(tests gtest gtest_main LinearAlgebra)