# CPU implementation sources
set(CPU_SOURCES
    cpu-main.cpp
    octree.cpp
    philox_rng.cpp
)

# Create CPU executable
add_executable(cpu_nbody ${CPU_SOURCES})
target_compile_options(cpu_nbody PRIVATE 
    $<$<CXX_COMPILER_ID:GNU>:-Wall -Wextra -O0 -g>
    $<$<CXX_COMPILER_ID:Clang>:-Wall -Wextra -O0 -g>
    $<$<CXX_COMPILER_ID:MSVC>:/W4 /Od /Zi>
)

# Create shared library for the core components (for tests to use)
add_library(nbody_core STATIC
    octree.cpp
    philox_rng.cpp
)
target_compile_options(nbody_core PRIVATE 
    $<$<CXX_COMPILER_ID:GNU>:-Wall -Wextra -O0 -g>
    $<$<CXX_COMPILER_ID:Clang>:-Wall -Wextra -O0 -g>
    $<$<CXX_COMPILER_ID:MSVC>:/W4 /Od /Zi>
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
    target_compile_options(gpu_nbody PRIVATE 
        $<$<COMPILE_LANGUAGE:CXX>:$<$<CXX_COMPILER_ID:GNU>:-Wall -Wextra -O0 -g>>
        $<$<COMPILE_LANGUAGE:CXX>:$<$<CXX_COMPILER_ID:Clang>:-Wall -Wextra -O0 -g>>
        $<$<COMPILE_LANGUAGE:CXX>:$<$<CXX_COMPILER_ID:MSVC>:/W4 /Od /Zi>>
    )
endif()
