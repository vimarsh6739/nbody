# CPU implementation sources
set(CPU_SOURCES
    cpu-main.cpp
    octree.cpp
    philox_rng.cpp
    cpu-core.cpp
)


# Check for OpenCilk (clang with cilk support)
# Use the environment variable to find OpenCilk's clang++
if(DEFINED ENV{CILK_CLANG})
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
target_link_options(cpu_nbody PRIVATE -fopencilk)
target_compile_options(cpu_nbody PRIVATE -fopencilk)


# Create shared library for the core components (for tests to use)
add_library(nbody_core STATIC
    octree.cpp
    philox_rng.cpp
    cpu-core.cpp
)



# GPU implementation (if CUDA is available)
if(CUDA_AVAILABLE)
    
    # GPU sources (CUDA and C++)
    set(GPU_SOURCES
        gpu-main.cu
        octree.cpp
        philox_rng.cpp
    )
    
    # Create GPU executable
    add_executable(gpu_nbody ${GPU_SOURCES})
    set_target_properties(gpu_nbody PROPERTIES 
        CUDA_SEPARABLE_COMPILATION ON
    )
    target_link_libraries(gpu_nbody PRIVATE CUDA::cudart)
endif()
