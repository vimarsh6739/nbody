CXX := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -I. -O2
LDFLAGS= 

OBJ_DIR := obj
BENCH_BIN := nbody

SRCS := $(wildcard *.cpp)

OBJS := $(SRCS:%.cpp=$(OBJ_DIR)/%.o)

all: directories $(BENCH_BIN)

directories:
	@mkdir -p $(OBJ_DIR)

$(BENCH_BIN): $(OBJS)
	@echo "Building benchmark program"
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) 

$(OBJ_DIR)/%.o: %.cpp
	@echo "Compiling $<..."
	@$(CXX) $(CXXFLAGS) -c $< -o $@ 

clean:
	@echo "Cleaning build files..."
	@rm -rf $(OBJ_DIR) $(BENCH_BIN)

rebuild: clean all

.PHONY: all clean rebuild directories 
