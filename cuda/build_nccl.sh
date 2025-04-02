nvcc -o nccl_hello nccl_hello.cu \
  -I$/home/vimarsh6739/.local/cuda/include -I/home/vimarsh6739/nccl/build/include \
  -L/home/vimarsh6739/.local/cuda/lib64 -L/home/vimarsh6739/nccl/build/lib \
  -lnccl -lcudart \
  -gencode=arch=compute_60,code=sm_60 \
  -gencode=arch=compute_61,code=sm_61 \
  -gencode=arch=compute_70,code=sm_70 \
  -std=c++11 -O3
