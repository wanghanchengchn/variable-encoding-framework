cmake_minimum_required(VERSION 3.5)
project(variable-encoding-framework LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra") # Warning levels
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mbmi -mbmi2 -mavx2 -O3")  # Bit Manipulation Instructions for VECF, AVX2 for vqf

add_subdirectory(src)
add_subdirectory(test)
