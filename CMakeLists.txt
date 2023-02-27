cmake_minimum_required(VERSION 3.0.0)
project(project VERSION 0.1.0)

include(CTest)
enable_testing()

add_library(project INTERFACE)
add_executable(SPE10 main.cpp vars.hpp f.hpp cell.hpp)
add_executable(PCG solver/pcg.cpp)

include_directories(${PROJECT_SOURCE_DIR})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
