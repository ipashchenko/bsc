cmake_minimum_required(VERSION 3.9)
project(bsc)

include_directories("/usr/include/eigen3", "/home/ilya/github/DNest4/code", "src")


set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=c++11 -march=native -DNDEBUG -w -DEIGEN_MPL2_ONLY")

set(SOURCE_FILES src/main.cpp src/Component.cpp src/Data.cpp src/Gain.cpp src/Gains.cpp src/SkyModel.cpp src/DNestModel.cpp
        src/MyExceptions.h)
add_executable(bsc ${SOURCE_FILES})

TARGET_LINK_LIBRARIES(bsc /home/ilya/github/DNest4/code/libdnest4.a)

set(TEST_DIR tests)
set(TEST_FILES ${TEST_DIR}/simple_test.cpp)
add_executable(bsc_test ${TEST_FILES})