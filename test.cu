// File : test.cu
#include <stdio.h>
#include <cuda.h>

int main()
{
    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);
    if(deviceCount == 0){
        printf("no CUDA compatible GPU exitsts.\n");
    }
    else
    {
        cudaDeviceProp pr;
        for(int i = 0; i<deviceCount;i++){
	cudaGetDeviceProperties(&pr, i);
	printf("Dev #%lu is %lf \n", sizeof(char), (ceil((double)-1/2)));
	}
    }
    return 1;
}
// cache l1 size is 48kB, warp size 32
// used https://xmartlabs.github.io/cuda-calculator/ to check occupancy of SMs
// compute_sanitizer per trovare errori silenti di cuda
// max thread per multiprocessor 1024
// max registers per sm is 65536