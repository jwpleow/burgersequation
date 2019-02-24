cmake_minimum_required(VERSION 3.13)
project(Coursework)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-O3")

find_package(BLAS REQUIRED)
find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})


add_executable(Coursework burgers.cpp Model.cpp Burgers.cpp)

target_link_libraries(Coursework ${BLAS_LIBRARIES} ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
    set_target_properties(Coursework PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(Coursework PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()