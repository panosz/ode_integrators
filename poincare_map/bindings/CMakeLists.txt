add_subdirectory(core)


# configure_file(src/test_action_integration_ext.py ${CMAKE_CURRENT_BINARY_DIR}/test_action_integration_ext.py COPYONLY)
# configure_file(src/test_action_integration_numerics.py ${CMAKE_CURRENT_BINARY_DIR}/test_action_integration_numerics.py COPYONLY)

configure_file(setup.py ${CMAKE_CURRENT_BINARY_DIR}/setup.py COPYONLY)
