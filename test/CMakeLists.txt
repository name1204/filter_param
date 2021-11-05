cmake_minimum_required(VERSION 2.8)

add_executable(cscade-iir-test main.cpp)
target_link_libraries(cscade-iir-test
    digital_filters
    )

add_test(
    NAME cscade-iir
    COMMAND $<TARGET_FILE:cscade-iir-test>
    )

# run with: ctest -L lib
set_property(
    TEST cscade-iir
    PROPERTY LABELS lib cscade-iir
    )