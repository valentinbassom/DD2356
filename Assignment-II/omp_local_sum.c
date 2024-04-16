#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define ARRAY_SIZE	10000000
#define N_TIMES 15

double A[ARRAY_SIZE];

void generate_random(double *input, size_t size){
  for (size_t i = 0; i < size; i++) {
    input[i] = rand() / (double)(RAND_MAX);
  }
}

double omp_local_sum(double *x, size_t size) {
    int n_threads = omp_get_max_threads();
    double local_sums[n_threads];
    memset(local_sums, 0, n_threads * sizeof(double));

    #pragma omp parallel shared(local_sums)
    {
        int id = omp_get_thread_num();
        #pragma omp for
        for (size_t i = 0; i < size; i++) {
            local_sums[id] += x[i];
        }
    }

    double sum_val = 0.0;
    for (size_t i = 0; i < n_threads; i++) {
        sum_val += local_sums[i];
    }

    return sum_val;
}

double sample_stand_dev(double *x, size_t size){
  double sum_val = 0.0;
  for (size_t i = 0; i < size; i++) {
    sum_val += x[i];
  }
  double mean = sum_val/size;

  double sum_squared = 0.0;
  for (size_t i = 0; i < size; i++) {
    sum_squared += pow(x[i] - mean, 2);
  }
  double std = sqrt(sum_squared/(size-1));

  return std;
}

int main(int argc, char *argv[]){

    int coherent = 1;
    double start, end, first_sum, sum, avg_time, std;
    double times[N_TIMES], results[N_TIMES];
    generate_random(A, ARRAY_SIZE);

    first_sum = omp_local_sum(A,ARRAY_SIZE);

    for(size_t i = 0; i < N_TIMES; i++){
        start = omp_get_wtime();
        sum = omp_local_sum(A,ARRAY_SIZE);
        end = omp_get_wtime();
        times[i] = end - start;
        results[i] = sum;
    }

    avg_time = 0.0;
    for(size_t i = 0; i < N_TIMES; i++){
        if(fabs(results[i] - first_sum) > 1e-6){
          printf("Error = %.5e\n",fabs(results[i] - first_sum));
          coherent = 0;
        }
        avg_time += times[i];
    }
    avg_time = avg_time/N_TIMES;
    std = sample_stand_dev(times, N_TIMES);

    if(coherent) printf("Results are coherent across all tries : YES\n");
    else printf("Results are coherent across all tries : NO\n");

    printf("Average Time = %.6f        Standard deviation = %.6f\n",avg_time,std);
    return 0;
}