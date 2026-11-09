/* TEST ONLY minimal CUDA call declarations, not an NVIDIA runtime/header. */
#pragma once
#include <stddef.h>
typedef int cudaError_t;
enum {cudaSuccess=0};
typedef enum {cudaMemcpyHostToDevice=1,cudaMemcpyDeviceToHost=2} cudaMemcpyKind;
extern "C" {
cudaError_t cudaDeviceSynchronize(void);
cudaError_t cudaMalloc(void**,size_t);
cudaError_t cudaFree(void*);
cudaError_t cudaMemcpy(void*,const void*,size_t,cudaMemcpyKind);
cudaError_t cudaMemset(void*,int,size_t);
const char*cudaGetErrorString(cudaError_t);
}
