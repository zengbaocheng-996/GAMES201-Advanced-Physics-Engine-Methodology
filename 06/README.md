**[Architecture-aware programming]** will become increasingly important in the future since processors and memory are stopping getting faster.

**[Parallelism: Computation is cheap]** Processors are more capable than you thought.

**[Locality: Communication is expensive]** Memory bandwidth is a more scarce resource nowadays

**Same rule apply to GPUs**

# Performance is key to high quality

Performance from algorithmic improvement (do less work)

- quick sort v.s. insertion sort; O(nlgn) v.s. O(n^2)
- Pre-conditioners for linear system solvers that make things converge faster

Performance from low-level programming (do work faster)

# Performance Engineering

### Hardware Architecture

###### Locality

Shrink the working set, so that data resides in lower-level (higher throughput, lower latency) memory

1. Spatial locality: try to access spatially neighboring data in main memory

   Higher cacheline utilization

   Fewer Cache/TLB misses

   Better harware prefetching on CPUs

2. Temporal locality: reuse the data as much as you can

   Higher cache-hit rates

   Lower main memory bandwidth pressure

###### Cachelines

memory.cpp

```cpp
constexpr int n = 256 * 1024 * 1024;
int a[n];
void benchmark()
{
    auto t = get_time();
    for(int i = 0; i < n; i += stride)
    {
        a[i] = i;
    }
    printf("%f\n", get_time() - t);
}
/* stride=1: 0.13s
   stride=2: 0.13s
   stride=4: 0.13s
   stride=8: 0.13s
   stride=16: 0.13s
   stride=32: 0.096s
   stride=64: 0.069s */
```

1.  different order

   0.244s

   ```cpp
   constexpr int n = 512;
   int A[n][n];
   int B[n][n];
   int C[n][n];
   for(int i = 0; i < n; i++)
       for(int j = 0; j < n; j++)
           for(int k = 0; k < n; k++)
               C[i][j] += A[i][k] * B[k][j]; // Spatial locality + Weak Spatial locality
   ```

   0.039s

   ```cpp
   constexpr int n = 512;
   int A[n][n];
   int B[n][n];
   int C[n][n];
   for(int i = 0; i < n; i++)
       for(int k = 0; k < n; k++)
           for(int j = 0; j < n; j++)
               C[i][j] += A[i][k] * B[k][j]; // Temporal locality + Spatial locality
   ```

2. different n (cache sets/cache tags)

   0.244s

   ```cpp
   constexpr int n = 512;
   int A[n][n];
   int B[n][n];
   int C[n][n];
   for(int i = 0; i < n; i++)
       for(int j = 0; j < n; j++)
           for(int k = 0; k < n; k++)
               C[i][j] += A[i][k] * B[k][j];
   ```

   0.135s

   ```cpp
   constexpr int n = 513;
   int A[n][n];
   int B[n][n];
   int C[n][n];
   for(int i = 0; i < n; i++)
       for(int j = 0; j < n; j++)
           for(int k = 0; k < n; k++)
               C[i][j] += A[i][k] * B[k][j];
   ```

###### CPU \muArch

1. Intel 64 and IA-32 Architectures Optimization Reference Manual
2. WikiChip

