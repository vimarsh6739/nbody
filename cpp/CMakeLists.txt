# CPU implementation sources
set(CPU_SOURCES
    cpu-core.cpp
)

# Create CPU executable
add_executable(cpu_nbody ${CPU_SOURCES})
target_link_libraries(cpu_nbody PRIVATE 
    nbody_common
    Threads::Threads 
    m
)

if(ENABLE_CILK)
  target_compile_options(cpu_nbody PRIVATE -fopencilk -v)
endif()

# Create shared library for the core components (for tests to use)
add_library(cpu_nbody_core STATIC
    cpu-core.cpp
)
target_link_libraries(cpu_nbody_core PRIVATE nbody_common m)

