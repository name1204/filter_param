cmake_minimum_required(VERSION 3.16)

# project name
project(digital_filters CXX)

#--------------- CMake compile config ------------------------------
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

#----------------- clang-format config -------------------------
find_program(CLANG_FORMAT_EXE clang-format)

file(GLOB_RECURSE project_sources CONFIGURE_DEPENDS "[a-zA-Z]*.cpp" "[a-zA-Z]*.hpp")
list(FILTER project_sources EXCLUDE REGEX "build/*")

# add build target (clang-format lints)
# In use,
# `cmake --build . --target ${PROJECT_NAME}-format --verbose`
add_custom_target("${PROJECT_NAME}-format"
  COMMAND ${CLANG_FORMAT_EXE} --style=file --dry-run ${project_sources}
  COMMENT "Running clang-format(lints)"
  VERBATIM
  )

# add build target (clang-format warning overwrite)
# In use,
# `cmake --build . --target ${PROJECT_NAME}-format-fix`
add_custom_target("${PROJECT_NAME}-format-fix"
  COMMAND ${CLANG_FORMAT_EXE} --style=file -i ${project_sources}
  COMMENT "Running clang-format(overwrite)"
  VERBATIM
  )

#------------------ build source config ------------------

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

# for tests
if(NOT without-test)
    enable_testing()
    add_subdirectory(test)
endif()
