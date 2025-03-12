#include <stdio.h>

__global__ void helloCUDA()
{
    printf("Hello, CUDA from !(%d,%d)\n",blockIdx.x,threadIdx.x);
    __syncthreads();
    if(blockIdx.x == 0 && threadIdx.x==0)
    printf("Goodbye. from !(%d,%d)\n",blockIdx.x,threadIdx.x);
}

int main()
{
    helloCUDA<<<1, 32>>>();
    cudaDeviceSynchronize();
    return 0;
}
