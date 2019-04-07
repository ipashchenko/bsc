cmake_minimum_required(VERSION 3.13)
project(bsc)

include_directories("src", "eigen", "DNest4/code")

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fext-numeric-literals")

set(SOURCE_FILES src/main.cpp src/Component.cpp src/Data.cpp src/Gain.cpp src/Gains.cpp src/SkyModel.cpp src/DNestModel.cpp)
add_executable(bsc ${SOURCE_FILES})

set(TEST_DIR tests)
set(TEST_FILES ${TEST_DIR}/simple_test.cpp)
add_executable(bsc_test ${TEST_FILES})