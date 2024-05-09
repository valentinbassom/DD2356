#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define SEED     921
#define NUM_ITER 1000000000

int main(int argc, char* argv[])
{
    int count = 0;
    int sum = 0;
    int rank, size, provided;
    double x, y, z, pi;
    
    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double t = MPI_Wtime();
    
    srand(SEED*rank); // Important: Multiply SEED by "rank" when you introduce MPI!
    
    // Calculate PI following a Monte Carlo method
    for (int iter = rank; iter < NUM_ITER; iter+=size)
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
    MPI_Reduce(&count, &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    pi = ((double)sum / (double)NUM_ITER) * 4.0;
    
    t = MPI_Wtime() - t;
    if(rank==0){
        printf("The result is %f\n", pi);
        printf("Execution time: %8.6f seconds\n", t);
    }

    MPI_Finalize();
    
    return 0;
}

