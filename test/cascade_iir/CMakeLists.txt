cmake_minimum_required(VERSION 2.8)

add_executable(cascade-iir-test cascade_iir_test.cpp)
target_link_libraries(cascade-iir-test digital_filters)
add_custom_command(
    TARGET cascade-iir-test
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E
        copy ${CMAKE_CURRENT_SOURCE_DIR}/desire_filter.csv
            $<TARGET_FILE_DIR:cascade-iir-test>)

add_test(
    NAME cascade-iir-BandParam_new
    COMMAND $<TARGET_FILE:cascade-iir-test> BandParam_new
    )

add_test(
    NAME cascade-iir-Band_generator
    COMMAND $<TARGET_FILE:cascade-iir-test> Band_generator
    )

add_test(
    NAME cascade-iir-analyze_edges
    COMMAND $<TARGET_FILE:cascade-iir-test> analyze_edges
    )

add_test(
    NAME cascade-iir-FilterParam_read_csv
    COMMAND $<TARGET_FILE:cascade-iir-test> FilterParam_read_csv
    )

add_test(
    NAME cascade-iir-FilterParam_csw
    COMMAND $<TARGET_FILE:cascade-iir-test> FilterParam_csw
    )

add_test(
    NAME cascade-iir-FilterParam_desire_res
    COMMAND $<TARGET_FILE:cascade-iir-test> FilterParam_desire_res
    )

add_test(
    NAME cascade-iir-FilterParam_freq_res_se
    COMMAND $<TARGET_FILE:cascade-iir-test> FilterParam_freq_res_se
    )

add_test(
    NAME cascade-iir-FilterParam_freq_res_so
    COMMAND $<TARGET_FILE:cascade-iir-test> FilterParam_freq_res_so
    )

add_test(
    NAME cascade-iir-FilterParam_freq_res_no
    COMMAND $<TARGET_FILE:cascade-iir-test> FilterParam_freq_res_no
    )

add_test(
    NAME cascade-iir-FilterParam_freq_res_mo
    COMMAND $<TARGET_FILE:cascade-iir-test> FilterParam_freq_res_mo
    )

add_test(
    NAME cascade-iir-Filter_param_group_delay_se
    COMMAND $<TARGET_FILE:cascade-iir-test> Filter_param_group_delay_se
    )

add_test(
    NAME cascade-iir-Filter_param_group_delay_so
    COMMAND $<TARGET_FILE:cascade-iir-test> Filter_param_group_delay_so
    )

add_test(
    NAME cascade-iir-Filter_param_group_delay_no
    COMMAND $<TARGET_FILE:cascade-iir-test> Filter_param_group_delay_no
    )

add_test(
    NAME cascade-iir-Filter_param_group_delay_mo
    COMMAND $<TARGET_FILE:cascade-iir-test> Filter_param_group_delay_mo
    )

add_test(
    NAME cascade-iir-FilterParam_judge_stability_even
    COMMAND $<TARGET_FILE:cascade-iir-test> FilterParam_judge_stability_even
    )

add_test(
    NAME cascade-iir-FilterParam_judge_stability_odd
    COMMAND $<TARGET_FILE:cascade-iir-test> FilterParam_judge_stability_odd
    )

add_test(
    NAME cascade-iir-FilterParam_evaluate_objective_function
    COMMAND $<TARGET_FILE:cascade-iir-test> FilterParam_evaluate_objective_function
    )

add_test(
    NAME cascade-iir-FilterParam_init_coef
    COMMAND $<TARGET_FILE:cascade-iir-test> FilterParam_init_coef
    )

add_test(
    NAME cascade-iir-FilterParam_init_stable_coef
    COMMAND $<TARGET_FILE:cascade-iir-test> FilterParam_init_stable_coef
    )

add_test(
    NAME cascade-iir-FilterParam_gprint_amp
    COMMAND $<TARGET_FILE:cascade-iir-test> FilterParam_gprint_amp
    )

add_test(
    NAME cascade-iir-FilterParam_gprint_mag
    COMMAND $<TARGET_FILE:cascade-iir-test> FilterParam_gprint_mag
    )
