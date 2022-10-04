set(TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/files)
set(TEST_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/io_test_script.cmake)
set(TEST_FILES test.xyz
               test.pwi
               test.pwo
               test.lmp
               test.dmp
               test.cpi
               test.cube
               test.xsf
               test.orca
               test.pos
               test.json
               #test.mol moltemplate is not necessarily available, how can this be tested?
)

foreach(FILE ${TEST_FILES})
    add_test(NAME convert_${FILE}
             COMMAND ${CMAKE_COMMAND} -D FILE=${TEST_DIR}/${FILE} -D VIPSTER=$<TARGET_FILE:vipster> -P ${TEST_SCRIPT})
endforeach()
