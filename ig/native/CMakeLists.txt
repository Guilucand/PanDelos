cmake_minimum_required(VERSION 3.0)
project(native)

set(CMAKE_CXX_STANDARD 14)

include_directories(jni jni/linux)

add_library(native SHARED library.cpp pangene_native.h)
