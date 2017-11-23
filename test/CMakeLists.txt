file(COPY test_data/ DESTINATION test_data/)
add_executable(runPengTestSuite main.h basic_tests.cpp test_base_pattern.cpp $<TARGET_OBJECTS:PENG_OBJECTS>)
target_link_libraries(runPengTestSuite gtest gtest_main shared)
add_test(NAME run_unit_tests COMMAND runPengTestSuite WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
