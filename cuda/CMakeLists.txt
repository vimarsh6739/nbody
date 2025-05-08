set(GPU_SRCS
  gpu-core.cu
)

# Create CPU executable
add_executable(gpu_nbody ${GPU_SRCS} gpu-main.cu)

# Create shared library for the core components (for unit tests)
add_library(gpu_nbody_core STATIC ${GPU_SRCS})

set(GPU_LIBS 
  nbody_common 
  CUDA::cudart
)

set_target_properties(gpu_nbody PROPERTIES 
  CUDA_SEPARABLE_COMPILATION ON
)

# adds "#define ENABLE_CUDA" during compilation
target_compile_definitions(cpu_nbody PRIVATE ENABLE_CUDA)
target_compile_definitions(cpu_nbody_core PRIVATE ENABLE_CUDA)

# link against the libs we need here
target_link_libraries(gpu_nbody PRIVATE ${GPU_LIBS})
target_link_libraries(gpu_nbody_core PRIVATE ${GPU_LIBS})

