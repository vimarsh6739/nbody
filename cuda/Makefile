NVCC := nvcc
CXX := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -I. -O2
NVCCFLAGS := -std=c++17 -I.
LDFLAGS := -lcudart -ldl -lrt -lpthread

OBJ_DIR := obj
BENCH_BIN := nbody_gpu

# Source files
CPP_SRCS := $(wildcard *.cpp)
CUDA_SRCS := $(wildcard *.cu)

# Object files
CPP_OBJS := $(CPP_SRCS:%.cpp=$(OBJ_DIR)/%.o)
CUDA_OBJS := $(CUDA_SRCS:%.cu=$(OBJ_DIR)/%.o)
OBJS := $(CPP_OBJS) $(CUDA_OBJS)

# If no .cu files exist yet, create one for the nbody implementation
ifeq ($(strip $(CUDA_SRCS)),)
    CUDA_SRCS := nbody.cu
    CUDA_OBJS := $(OBJ_DIR)/nbody.o
    OBJS := $(CPP_OBJS) $(CUDA_OBJS)
endif

all: directories $(BENCH_BIN)

# Create necessary directories
directories:
	@mkdir -p $(OBJ_DIR)

# Build the final executable
$(BENCH_BIN): $(OBJS)
	@echo "Building CUDA benchmark program"
	$(NVCC) $(NVCCFLAGS) $^ -o $@ $(LDFLAGS)

# Compile C++ source files
$(OBJ_DIR)/%.o: %.cpp
	@echo "Compiling C++ file $<..."
	@$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile CUDA source files
$(OBJ_DIR)/%.o: %.cu
	@echo "Compiling CUDA file $<..."
	@$(NVCC) $(NVCCFLAGS) -c $< -o $@

# Create nbody.cu if it doesn't exist
$(CUDA_SRCS):
	@if [ ! -f $@ ]; then \
		echo "Creating $@..."; \
		cp nbody-gpu.cpp $@; \
	fi

# Clean build artifacts
clean:
	@echo "Cleaning build files..."
	@rm -rf $(OBJ_DIR) $(BENCH_BIN)

# Clean and rebuild
rebuild: clean all

.PHONY: all clean rebuild directories
