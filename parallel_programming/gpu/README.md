# GPU

Manual build on lucia
```
module load Clang/16.0.6-GCCcore-11.3.0-CUDA-11.7.0
clang -O3 -fopenmp -fopenmp-targets=nvptx64-nvidia-cuda -o saxpy_gpu saxpy_gpu.c
```