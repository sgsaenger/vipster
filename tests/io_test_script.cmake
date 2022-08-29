cmake_minimum_required(VERSION 3.14)

get_filename_component(FN ${FILE} NAME)
get_filename_component(EXT ${FILE} EXT)
string(SUBSTRING ${EXT} 1 -1 EXT)

# run first pass: convert from source to xyz
execute_process(COMMAND ${VIPSTER} convert ${EXT} ${TEST_DIR}/${FILE} xyz ${FN}.V1
                RESULT_VARIABLE ERROR)
if(ERROR)
    message(FATAL_ERROR "Converting file ${FN} failed: ${ERROR}")
endif()

# run second pass: convert new xyz once more
execute_process(COMMAND ${VIPSTER} convert xyz ${FN}.V1 xyz ${FN}.V2
                RESULT_VARIABLE ERROR)
if(ERROR)
    message(FATAL_ERROR "Reparsing converted file ${FN} failed: ${ERROR}")
endif()

# verify written file
execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files --ignore-eol ${FN}.V2 ${FN}.V1 #${TEST_DIR}/test.xyz
                RESULT_VARIABLE ERROR)
if(ERROR)
    message(FATAL_ERROR "File ${FN} differs after conversion: ${ERROR}")
endif()

# cleanup temp file
execute_process(COMMAND ${CMAKE_COMMAND} -E rm ${FN}.V1 ${FN}.V2
                RESULT_VARIABLE ERROR)
if(ERROR)
    message(FATAL_ERROR "Couldn't remove files: ${ERROR}")
endif()
