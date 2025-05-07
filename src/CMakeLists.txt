# CPU implementation sources
set(CPU_SOURCES
    cpu-main.cpp
    octree.cpp
    philox_rng.cpp
    cpu-core.cpp
)

# Create CPU executable
add_executable(cpu_nbody ${CPU_SOURCES})
target_link_libraries(cpu_nbody PRIVATE Threads::Threads)

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
