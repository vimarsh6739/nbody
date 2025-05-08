# Create test executable
add_executable(cpu_test
    cpu-test.cpp
)

# Link against our core library and Google Test
target_link_libraries(cpu_test
    PRIVATE
    cpu_nbody_core
    gtest
    gtest_main
)

if(ENABLE_CILK)
  target_compile_definitions(cpu_test PRIVATE ENABLE_CILK)
  target_compile_options(cpu_test PRIVATE -fopencilk)
  target_link_options(cpu_test PRIVATE -fopencilk) 
endif()

# Add the test
add_test(NAME cpu_test COMMAND cpu_test)
