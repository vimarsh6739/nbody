#ifndef CUDA_UTILS_H
#define CUDA_UTILS_H

#define CUDA_CALL(call)                                                        \
  do {                                                                         \
    cudaError_t err = call;                                                    \
    if (err != cudaSuccess) {                                                  \
      fprintf(stderr, "CUDA error %s@%d: %s (%d)\n", __FILE__, __LINE__,       \
              cudaGetErrorString(err), err);                                   \
      exit(EXIT_FAILURE);                                                      \
    }                                                                          \
  } while (0)


#endif // !CUDA_UTILS_H
