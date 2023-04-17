#include <cstdio>
#include <cstdlib>
#include <vector>

int main() {
  int n = 50;
  int range = 5;
  std::vector<int> key(n);
  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]);
  }
  printf("\n");


  std::vector<int> bucket(range,0); 
  for (int i=0; i<n; i++)
    bucket[key[i]]++;
  std::vector<int> offset(bucket);
  std::vector<int> tmp_offset(range, 0);
  //for (int i=0; i<range; i++){
  //  printf("%d ",bucket[i]);
  //}
  //printf("\n");
#pragma omp prallel
  for (int k=1; k<range; k<<=1) {
#pragma omp for
    for (int i=0; i<range; i++)
      tmp_offset[i] = offset[i];
#pragma omp for 
    for (int i=k; i<range; i++)
      offset[i] += tmp_offset[i-k];
  }
  offset.insert(offset.begin(), 0);
  for (int i=0; i<range; i++) {
    int j = offset[i];
    for (; bucket[i]>0; bucket[i]--) {
      key[j++] = i;
    }
  }
  

  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");
}
