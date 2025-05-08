set(CPU_EXE_SRCS
    cpu-core.cpp
    cpu-main.cpp
)

# Create CPU executable
add_executable(cpu_nbody ${CPU_EXE_SRCS})

# Create shared library for the core components (for unit tests)
add_library(cpu_nbody_core STATIC cpu-core.cpp)

set(CPU_LIBS nbody_common Threads::Threads ${GPERFTOOLS_PROFILER})

target_link_libraries(cpu_nbody PRIVATE ${CPU_LIBS})
target_link_libraries(cpu_nbody_core PRIVATE ${CPU_LIBS})

if(ENABLE_CILK)
  target_compile_definitions(cpu_nbody PRIVATE ENABLE_CILK)
  target_compile_definitions(cpu_nbody_core PRIVATE ENABLE_CILK)
  
  target_compile_options(cpu_nbody PRIVATE -fopencilk)
  target_compile_options(cpu_nbody_core PRIVATE -fopencilk)
  
  target_link_options(cpu_nbody PRIVATE -fopencilk) 
  target_link_options(cpu_nbody_core PRIVATE -fopencilk) 
endif()
