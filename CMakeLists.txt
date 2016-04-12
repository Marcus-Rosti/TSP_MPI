cmake_minimum_required(VERSION 2.6)
project(TSP_MPI)
# The version number.
set(TSP_MPI_VERSION_MAJOR 0)
set(TSP_MPI_VERSION_MINOR 1)

find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

if ("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
    # using Clang
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99 -Ofast -Weverything -pedantic")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    # using Gcc
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99 -O3 -Wall")
endif ()

# add the executable
add_executable(TSP_MPI main.c)

target_link_libraries(TSP_MPI ${MPI_LIBRARIES})

if (MPI_COMPILE_FLAGS)
    set_target_properties(TSP_MPI PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif ()

if (MPI_LINK_FLAGS)
    set_target_properties(TSP_MPI PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")
endif ()