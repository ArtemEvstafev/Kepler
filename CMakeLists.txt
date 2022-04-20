cmake_minimum_required(VERSION 3.17)
project(Kepler)


set(CMAKE_CXX_STANDARD 17)

enable_testing()
add_subdirectory(tests)
add_library(Kepler STATIC Kepler.cpp Kepler.hpp)