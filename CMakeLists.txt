cmake_minimum_required(VERSION 3.16)

project(digital_filters CXX)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 11) # C++11
set(CMAKE_CXX_STANDARD_REQUIRED ON) #...is required...

set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic -Wcast-align -Wcast-qual -Wconversion -Wdisabled-optimization -Wendif-labels -Wfloat-equal -Winit-self -Winline -Wmissing-include-dirs -Wnon-virtual-dtor -Wold-style-cast -Woverloaded-virtual -Wpacked -Wpointer-arith -Wredundant-decls -Wshadow -Wsign-promo -Wswitch-default -Wswitch-enum -Wvariadic-macros -Wwrite-strings ")
set(CMAKE_CXX_FLAGS_DEBUG "-g3 -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g3 -Og")
set(CMAKE_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG -march=native")
if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pg")
        set(CMAKE_EXE_LINKER_FLAGS_DEBUG "-pg")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -pg")
        set(CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO "-pg")
endif(CMAKE_COMPILER_IS_GNUCXX)


# include path
include_directories(include)

# sub directory
add_subdirectory(lib)

# libraries
add_library(digital_filters INTERFACE)
    target_link_libraries(digital_filters INTERFACE cascade_iir)

# for development executable(optional activate)
option(DEVELOP_EXECUTABLE_DIGITAL_FILTERS "Build develop executable file on digital filters package" OFF)
if(DEVELOP_EXECUTABLE_DIGITAL_FILTERS)
# executables
    add_executable(develop_main develop_main.cpp)
    target_link_libraries(develop_main digital_filters)
endif()

# tests
if(NOT without-test)
    enable_testing()
    add_subdirectory(test)
endif()
