#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define SEED     921
#define NUM_ITER 1000000000

int main(int argc, char* argv[])
{
    int count = 0;
    double x, y, z, pi;
    
    clock_t t = clock();
    
    srand(SEED); // Important: Multiply SEED by "rank" when you introduce MPI!
    
    // Calculate PI following a Monte Carlo method
    for (int iter = 0; iter < NUM_ITER; iter++)
    {
        // Generate random (X,Y) points
        x = (double)random() / (double)RAND_MAX;
        y = (double)random() / (double)RAND_MAX;
        z = x*x + y*y;
        
        // Check if point is in unit circle
        if (z <= 1.0)
        {
            count++;
        }
    }
    
    // Estimate Pi and display the result
    pi = ((double)count / (double)NUM_ITER) * 4.0;
    
    t = clock() - t;
    printf("The result is %f\n", pi);
    printf("Execution time: %8.6f seconds\n", (double) t/CLOCKS_PER_SEC);
    
    return 0;
}

