cmake_minimum_required(VERSION 3.7)
project(ImageCompression)

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES main.cpp image_util.h image_util.cpp)
add_executable(ImageCompression ${SOURCE_FILES})