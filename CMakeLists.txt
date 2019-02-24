cmake_minimum_required(VERSION 3.13)
project(Coursework)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-O3")

find_package(BLAS REQUIRED)

add_executable(Coursework burgers.cpp Model.cpp Burgers.cpp)

target_link_libraries(Coursework ${BLAS_LIBRARIES})
