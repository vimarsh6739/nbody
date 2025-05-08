cmake_minimum_required(VERSION 3.18)


option(ENABLE_CILK "Enable Cilk based parallelization" ON)
# Check for OpenCilk (clang with cilk support)
if(ENABLE_CILK)
    if(DEFINED ENV{CILK_CLANG})
        set(CILK_COMPILER $ENV{CILK_CLANG})
        message(STATUS "OpenCilk compiler found: ${CILK_COMPILER}")
        set(CMAKE_CXX_COMPILER ${CILK_COMPILER})
    else()
        message(WARNING "OpenCilk compiler not found, building without Cilk")
    endif()
endif()

project(nbody_simulation LANGUAGES CXX C)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set default build type if not specified
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug")
endif()

# Set compiler flags for different build types
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -fno-omit-frame-pointer")

# system threads
find_package(Threads REQUIRED)

# Check for CUDA
find_package(CUDAToolkit)
if(CUDAToolkit_FOUND AND DEFINED ENV{ENABLE_CUDA})
  enable_language(CUDA)
    set(CUDA_AVAILABLE TRUE)
    
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    
    set(CMAKE_CUDA_ARCHITECTURES 80 86 87)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Wno-deprecated-gpu-targets")
    
    if(ENABLE_CUDA_DEBUG)
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -G") # enable cuda-gdb (may significantly affect performance on some targets)
    else()
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -lineinfo") # add line information to all builds for debug tools (exclusive to -G option)
    endif()
    
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --extended-lambda")
    
    message(STATUS "CUDA available, building GPU version")

else()
    set(CUDA_AVAILABLE FALSE)
    message(STATUS "CUDA not available (or not enabled. set ENABLE_CUDA), skipping GPU version")
endif()

# Options
option(BUILD_TESTING "Build tests" ON)

# Add Google Test
include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.17.0
)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/common)

# Add sources
add_subdirectory(common)
add_subdirectory(cpp)

if(CUDAToolkit_FOUND AND DEFINED ENV{ENABLE_CUDA})
  add_subdirectory(cuda)
endif()

# Testing
if(BUILD_TESTING)
    FetchContent_MakeAvailable(googletest)
    enable_testing()
    add_subdirectory(test)
endif()
