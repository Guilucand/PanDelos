cmake_minimum_required(VERSION 3.2)
project(PanDelosNative)

set(CMAKE_CXX_STANDARD 14)

include_directories(include)

find_package(Threads)
add_executable(PanDelos src/main.cpp src/pangene_idata.cpp src/algorithm.cpp src/stats.cpp)
target_link_libraries(PanDelos Threads::Threads)