# N body simulations

N body simulations implemented as a project for CS533 : Parallel Computer Architecture

## Build and run

### Basic cpp implementation
1. `mkdir build`
2. `cd build`
3. `cmake ..`
4. `make`

### Cilk implementation
1. Build and install OpenCILK by following the instructions at [OpenCILK infra](https://github.com/OpenCilk/infrastructure/blob/release/INSTALLING.md). Alternatively, just run `./install_cilk.sh`
2. `export CILK_CLANG=/path/to/opencilk/clang`
3. Build via cmake as above. Ensure that the cmake output has correctly identified and linked Cilk.
4. The nbody executable can be run via `CILK_NWORKERS=1 ./src/cpu_nbody`

### GPU implementation
TBD: this implementation is not done. However, can be built via:
1. `export ENABLE_CUDA=1`
2. Follow build instructions above. Ensure that the cmake output has correctly identified and linked CUDA.
