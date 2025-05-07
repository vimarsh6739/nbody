# CPU implementation sources
set(GPU_SOURCES
    gpu-main.cpp
    gpu-core.cpp
    ../common/octree.cpp
    ../common/philox_rng.cpp
)


# Create shared library for the core components (for tests to use)
add_library(gpu_nbody_core STATIC
    ../common/octree.cpp
    ../common/philox_rng.cpp
    gpu-core.cpp
    gpu-main.cpp
)

if (NOT CUDA_AVAILABLE)
    message(FATAL_ERROR "CUDA is not available. GPU implementation cannot be built.")   
endif()


# Create GPU executable
add_executable(gpu_nbody ${GPU_SOURCES})
set_target_properties(gpu_nbody PROPERTIES 
    CUDA_SEPARABLE_COMPILATION ON
)
target_link_libraries(gpu_nbody PRIVATE CUDA::cudart)

