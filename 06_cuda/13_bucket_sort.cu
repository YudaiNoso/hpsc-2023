#include <cstdio>
#include <cstdlib>
#include <vector>

__global__ void initialize(int *bucket) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  bucket[i] = 0;
}

__global__ void reduction(int *bucket, int *key) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  atomicAdd(&bucket[key[i]], 1);
}

__global__ void makeoffset(int *bucket, int *tmp, int range) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  for(int j=1; j<range; j<<=1){
    tmp[i] = bucket[i];
    __syncthreads();
    bucket[i] += tmp[i-j];
    __syncthreads();
  }
}

__global__ void sort(int *bucket, int *key, int range) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  //int bucket[5] = {8, 18, 32, 41, 50};
  for(int j=range-1; j>=0; j--){
    if(i>=bucket[j]) return;
    key[i] = j;
  }

}

int main() {
  int n = 50;
  int range = 5;
  int *key;
  cudaMallocManaged(&key, n*sizeof(int));
  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]);
  }
  printf("\n");
  
  int *bucket;
  cudaMallocManaged(&bucket, range*sizeof(int));
  initialize<<<1, range>>>(bucket);
  cudaDeviceSynchronize();
  
  reduction<<<1, n>>>(bucket, key);
  cudaDeviceSynchronize();
  
  int *tmp;
  cudaMallocManaged(&tmp, range*sizeof(int));
  makeoffset<<<1, range>>>(bucket, tmp, range);
  cudaDeviceSynchronize();

  sort<<<1, n>>>(bucket, key, range);
  cudaDeviceSynchronize();

  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");
  cudaFree(key);
  cudaFree(bucket);
  cudaFree(tmp);
}
