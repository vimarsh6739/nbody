# Create test executable
add_executable(cpu_test
    cpu-test.cpp
)

# Link against our core library and Google Test
target_link_libraries(cpu_test
    PRIVATE
    nbody_core
    gtest
    gtest_main
)

# Add the test
add_test(NAME cpu_test COMMAND cpu_test)
