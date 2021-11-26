cmake_minimum_required(VERSION 3.20)

project(digital_filters CXX)

set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wextra -pedantic -Wcast-align -Wcast-qual -Wconversion -Wdisabled-optimization -Wendif-labels -Wfloat-equal -Winit-self -Winline -Wmissing-include-dirs -Wnon-virtual-dtor -Wold-style-cast -Woverloaded-virtual -Wpacked -Wpointer-arith -Wredundant-decls -Wshadow -Wsign-promo -Wswitch-default -Wswitch-enum -Wvariadic-macros -Wwrite-strings ")
set(CMAKE_CXX_FLAGS_DEBUG "-g3 -O0 -pg")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -s -DNDEBUG -march=native")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g3 -Og -pg")
set(CMAKE_CXX_FLAGS_MINSIZEREL "-Os -s -DNDEBUG -march=native")

# include path
include_directories(include)

# sub directory
add_subdirectory(lib)

# libraries
add_library(digital_filters INTERFACE)
    target_link_libraries(digital_filters INTERFACE cascade_iir)

# executables
add_executable(develop_main develop_main.cpp)
target_link_libraries(develop_main digital_filters)

# tests
if(NOT without-test)
    enable_testing()
    add_subdirectory(test)
endif()
