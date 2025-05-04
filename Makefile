CILK_C := $(realpath opencilk-build/bin/clang)
CILK_CXX := opencilk-build/bin/clang++
CILK_FLAGS := -fopencilk

CXX := g++
NVCC := nvcc
CXXFLAGS := -std=c++17 -Wall -Wextra -O2 -g -I.
LDFLAGS= 

# CPU implementation
CPU_SRCS := nbody.cpp \
	    octree.cpp \
	    philox_rng.cpp

CPU_OBJDIR := obj
CPU_OBJS := $(patsubst %.cpp,$(CPU_OBJDIR)/%.o,$(CPU_SRCS))

# GPU implementation
GPU_SRCS := nbody_gpu.cu \
	    octree.cpp \
	    philox_rng.cpp

GPU_OBJDIR := obj_gpu
GPU_OBJS := $(patsubst %.cpp,$(GPU_OBJDIR)/%.o,$(filter %.cpp,$(GPU_SRCS))) \
            $(patsubst %.cu,$(GPU_OBJDIR)/%.o,$(filter %.cu,$(GPU_SRCS)))

all: directories cpu_nbody

directories:
	@mkdir -p $(CPU_OBJDIR) $(GPU_OBJDIR)

cpu_nbody: $(CPU_OBJS) 
	@printf "Building %s > %s\n" "$^" $@
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(CPU_OBJDIR)/%.o: %.cpp
	@printf "Checking %s\n" $@
	$(if $(filter $(basename $<).cpp,$(CPU_SRCS)),\
		@printf "CPU build %-35s > %s\n" $< $@;\
		$(CXX) $(CXXFLAGS) -c $< -o $@ $(LDFLAGS),\
		@true)

$(GPU_OBJDIR)/%.o: %.cpp
	@printf "Checking %s\n" $@
	$(if $(filter $(basename $<).cpp,$(GPU_SRCS)),\
		@printf "GPU cpp-build %-35s > %s\n" $< $@";\
		$(CXX) $(CXXFLAGS) -c $< -o $@ $(LDFLAGS),\
		@true)

$(GPU_OBJDIR)/%.o: %.cu
	@echo "check:: $<..."
	$(if $(filter $(basename $<).cu,$(GPU_SRCS)),\
		@printf "GPU cuda-build %-35s > %s\n" $< $@";\
		$(NVCC) $(CXXFLAGS) -c $< -o $@ $(LDFLAGS),\
		@true)
clean:
	rm -f $(CPU_OBJS) $(GPU_OBJS)
	rm -rf $(CPU_OBJDIR) $(GPU_OBJDIR)
	rm -f cpu_nbody gpu_nbody

.PHONY: all clean rebuild directories
