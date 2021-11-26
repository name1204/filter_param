cmake_minimum_required(VERSION 3.20)

add_executable(cascade-iir-test cascade_iir_test.cpp)
target_link_libraries(cascade-iir-test digital_filters)
# prepared necessary files for test
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
    set_property(
        TEST cascade-iir-BandParam_new
        PROPERTY LABELS lib cscade-iir cascade-iir-BandParam_new
        )

add_test(
    NAME cascade-iir-Band_generator
    COMMAND $<TARGET_FILE:cascade-iir-test> Band_generator
    )
    set_property(
        TEST cascade-iir-Band_generator
        PROPERTY LABELS lib cscade-iir cascade-iir-Band_generator
        )

add_test(
    NAME cascade-iir-analyze_edges
    COMMAND $<TARGET_FILE:cascade-iir-test> analyze_edges
    )
    set_property(
        TEST cascade-iir-analyze_edges
        PROPERTY LABELS lib cscade-iir cascade-iir-analyze_edges
        )

add_test(
    NAME cascade-iir-FilterParam_read_csv
    COMMAND $<TARGET_FILE:cascade-iir-test> FilterParam_read_csv
    )
    set_property(
        TEST cascade-iir-FilterParam_read_csv
        PROPERTY LABELS lib cscade-iir cascade-iir-FilterParam_read_csv
        )

add_test(
    NAME cascade-iir-FilterParam_csw
    COMMAND $<TARGET_FILE:cascade-iir-test> FilterParam_csw
    )
    set_property(
        TEST cascade-iir-FilterParam_csw
        PROPERTY LABELS lib cscade-iir cascade-iir-FilterParam_csw
        )

add_test(
    NAME cascade-iir-FilterParam_desire_res
    COMMAND $<TARGET_FILE:cascade-iir-test> FilterParam_desire_res
    )
    set_property(
        TEST cascade-iir-FilterParam_desire_res
        PROPERTY LABELS lib cscade-iir cascade-iir-FilterParam_desire_res
        )

add_test(
    NAME cascade-iir-FilterParam_freq_res_se
    COMMAND $<TARGET_FILE:cascade-iir-test> FilterParam_freq_res_se
    )
set_property(
        TEST cascade-iir-FilterParam_freq_res_se
        PROPERTY LABELS lib cscade-iir cascade-iir-FilterParam_freq_res_se
        )

add_test(
    NAME cascade-iir-FilterParam_freq_res_so
    COMMAND $<TARGET_FILE:cascade-iir-test> FilterParam_freq_res_so
    )
    set_property(
        TEST cascade-iir-FilterParam_freq_res_so
        PROPERTY LABELS lib cscade-iir cascade-iir-FilterParam_freq_res_so
        )

add_test(
    NAME cascade-iir-FilterParam_freq_res_no
    COMMAND $<TARGET_FILE:cascade-iir-test> FilterParam_freq_res_no
    )
    set_property(
        TEST cascade-iir-FilterParam_freq_res_no
        PROPERTY LABELS lib cscade-iir cascade-iir-FilterParam_freq_res_no
        )

add_test(
    NAME cascade-iir-FilterParam_freq_res_mo
    COMMAND $<TARGET_FILE:cascade-iir-test> FilterParam_freq_res_mo
    )
    set_property(
        TEST cascade-iir-FilterParam_freq_res_mo
        PROPERTY LABELS lib cscade-iir cascade-iir-FilterParam_freq_res_mo
        )

add_test(
    NAME cascade-iir-Filter_param_group_delay_se
    COMMAND $<TARGET_FILE:cascade-iir-test> Filter_param_group_delay_se
    )
    set_property(
        TEST cascade-iir-Filter_param_group_delay_se
        PROPERTY LABELS lib cscade-iir cascade-iir-Filter_param_group_delay_se
        )

add_test(
    NAME cascade-iir-Filter_param_group_delay_so
    COMMAND $<TARGET_FILE:cascade-iir-test> Filter_param_group_delay_so
    )
    set_property(
        TEST cascade-iir-Filter_param_group_delay_so
        PROPERTY LABELS lib cscade-iir cascade-iir-Filter_param_group_delay_so
        )

add_test(
    NAME cascade-iir-Filter_param_group_delay_no
    COMMAND $<TARGET_FILE:cascade-iir-test> Filter_param_group_delay_no
    )
    set_property(
        TEST cascade-iir-Filter_param_group_delay_no
        PROPERTY LABELS lib cscade-iir cascade-iir-Filter_param_group_delay_no
        )

add_test(
    NAME cascade-iir-Filter_param_group_delay_mo
    COMMAND $<TARGET_FILE:cascade-iir-test> Filter_param_group_delay_mo
    )
    set_property(
        TEST cascade-iir-Filter_param_group_delay_mo
        PROPERTY LABELS lib cscade-iir cascade-iir-Filter_param_group_delay_mo
        )

add_test(
    NAME cascade-iir-FilterParam_judge_stability_even
    COMMAND $<TARGET_FILE:cascade-iir-test> FilterParam_judge_stability_even
    )
    set_property(
        TEST cascade-iir-FilterParam_judge_stability_even
        PROPERTY LABELS lib cscade-iir cascade-iir-FilterParam_judge_stability_even
        )

add_test(
    NAME cascade-iir-FilterParam_judge_stability_odd
    COMMAND $<TARGET_FILE:cascade-iir-test> FilterParam_judge_stability_odd
    )
    set_property(
        TEST cascade-iir-FilterParam_judge_stability_odd
        PROPERTY LABELS lib cscade-iir cascade-iir-FilterParam_judge_stability_odd
        )

add_test(
    NAME cascade-iir-FilterParam_evaluate_objective_function
    COMMAND $<TARGET_FILE:cascade-iir-test> FilterParam_evaluate_objective_function
    )
    set_property(
        TEST cascade-iir-FilterParam_evaluate_objective_function
        PROPERTY LABELS lib cscade-iir cascade-iir-FilterParam_evaluate_objective_function
        )

add_test(
    NAME cascade-iir-FilterParam_init_coef
    COMMAND $<TARGET_FILE:cascade-iir-test> FilterParam_init_coef
    )
    set_property(
        TEST cascade-iir-FilterParam_init_coef
        PROPERTY LABELS lib cscade-iir cascade-iir-FilterParam_init_coef
        )

add_test(
    NAME cascade-iir-FilterParam_init_stable_coef
    COMMAND $<TARGET_FILE:cascade-iir-test> FilterParam_init_stable_coef
    )
    set_property(
        TEST cascade-iir-FilterParam_init_stable_coef
        PROPERTY LABELS lib cscade-iir cascade-iir-FilterParam_init_stable_coef
        )

add_test(
    NAME cascade-iir-FilterParam_gprint_amp
    COMMAND $<TARGET_FILE:cascade-iir-test> FilterParam_gprint_amp
    )
    set_property(
        TEST cascade-iir-FilterParam_gprint_amp
        PROPERTY LABELS lib cscade-iir cascade-iir-FilterParam_gprint_amp
        )

add_test(
    NAME cascade-iir-FilterParam_gprint_mag
    COMMAND $<TARGET_FILE:cascade-iir-test> FilterParam_gprint_mag
    )
    set_property(
        TEST cascade-iir-FilterParam_gprint_mag
        PROPERTY LABELS lib cscade-iir cascade-iir-FilterParam_gprint_mag
        )
