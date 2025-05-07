# CPU implementation sources
set(CPU_SOURCES
    ../common/main.cpp
    ../common/octree.cpp
    ../common/philox_rng.cpp
    cpu-core.cpp
)


# Check for OpenCilk (clang with cilk support)
# Use the environment variable to find OpenCilk's clang++
if(DEFINED ENV{CILK_CLANG} AND ("$ENV{ENABLE_CILK}" STREQUAL "1"))
    set(CILK_COMPILER $ENV{CILK_CLANG})
    message(STATUS "OpenCilk compiler found: ${CILK_COMPILER}")
    
    # if(DEFINED ENV{CILK_INCLUDE_DIR})
    #     message(STATUS "CILK_INCLUDE_DIR: $ENV{CILK_INCLUDE_DIR}")
    #     include_directories($ENV{CILK_INCLUDE_DIR})
    # else()
    #     message(WARNING "CILK_INCLUDE_DIR environment variable not set. Please set it manually.")
    # endif()

    set(CMAKE_CXX_COMPILER ${CILK_COMPILER})
    add_compile_options(-fopencilk)  # Enable OpenCilk extensions
         # Link with OpenCilk runtime (cilkrts)
    
else()
    message(WARNING "OpenCilk compiler not found, building without Cilk")
endif()

# Create CPU executable
add_executable(cpu_nbody ${CPU_SOURCES})
target_link_libraries(cpu_nbody PRIVATE Threads::Threads)

if(DEFINED ENV{CILK_CLANG} AND ("$ENV{ENABLE_CILK}" STREQUAL "1"))
  target_link_options(cpu_nbody PRIVATE -fopencilk)
  target_compile_options(cpu_nbody PRIVATE -fopencilk)
  target_compile_options(cpu_nbody PRIVATE -DENABLE_CILK)
endif()


# Create shared library for the core components (for tests to use)
add_library(cpu_nbody_core STATIC
    ../common/octree.cpp
    ../common/philox_rng.cpp
    cpu-core.cpp
    ../common/main.cpp
)

