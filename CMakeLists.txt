cmake_minimum_required(VERSION 3.10)

project(ATtacker)
set(CMAKE_CXX_STANDARD 14)

include(ExternalProject)
set(EXTERNAL_INSTALL_LOCATION ${CMAKE_BINARY_DIR}/external)

include_directories(${EXTERNAL_INSTALL_LOCATION}/include)
link_directories(${EXTERNAL_INSTALL_LOCATION}/lib)

include_directories(.)

# Library module
add_library(LinearAlgebra SHARED
        lll.h
        Matrix.h
        Vector.h
        Manger.h
        Main.cpp)

target_link_libraries(LinearAlgebra gmp gmpxx)

# ATtacker module
add_executable(ATtacker
        Main.cpp)

add_dependencies(ATtacker LinearAlgebra)
# target_link_libraries(ATtacker LinearAlgebra)

# Tests config
include(cmake/googletest.cmake)
fetch_googletest(
        ${PROJECT_SOURCE_DIR}/cmake
        ${PROJECT_BINARY_DIR}/googletest
)
enable_testing()

# Test module
add_executable(unit_tests vector_tests.cpp matrix_tests.cpp lll_tests.cpp Manger_tests.cpp lll.cpp)
target_link_libraries(unit_tests LinearAlgebra gtest_main pthread)
# add_test(NAME unit COMMAND ${CMAKE_BINARY_DIR}/unit_tests)

include(GoogleTest)
gtest_add_tests(TARGET unit_tests)


