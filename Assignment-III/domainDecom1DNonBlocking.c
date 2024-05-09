
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[]){

    int rank, size, i, provided;
    double commStart, commEnd, time;
    MPI_Status stats[2];
    MPI_Request reqs[2];
    
    // number of cells (global)
    int nxc = 128; // make sure nxc is divisible by size
    double L = 2*3.1415; // Length of the domain
    int tagRight = 1;
    int tagLeft = 2;
    

    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // number of nodes (local to the process): 0 and nxn_loc-1 are ghost cells 
    int nxn_loc = nxc/size + 3; // number of nodes is number cells + 1; we add also 2 ghost cells
    double L_loc = L/((double) size);
    double dx = L / ((double) nxc);
    
    // define out function
    double *f = calloc(nxn_loc, sizeof(double)); // allocate and fill with z
    double *dfdx = calloc(nxn_loc, sizeof(double)); // allocate and fill with z

    for (i=1; i<(nxn_loc-1); i++)
      f[i] = sin(L_loc*rank + (i-1) * dx);
    
    commStart = MPI_Wtime();
    // need to communicate and fill ghost cells f[0] and f[nxn_loc-1]
    // communicate ghost cells
    if (rank%2==0){
      // send right
      MPI_Isend(&f[nxn_loc-3], 1, MPI_DOUBLE, (rank+1)%size, tagRight, MPI_COMM_WORLD, reqs);
      MPI_Irecv(f, 1, MPI_DOUBLE, (rank > 0)?(rank-1):(size-1), tagRight, MPI_COMM_WORLD, reqs);
      // send left
      MPI_Isend(&f[2], 1, MPI_DOUBLE, (rank > 0)?(rank-1):(size-1), tagLeft, MPI_COMM_WORLD, reqs+1);
      MPI_Irecv(&f[nxn_loc-1], 1, MPI_DOUBLE, (rank+1)%size, tagLeft, MPI_COMM_WORLD, reqs+1);
    } else {
      // send right
      MPI_Irecv(f, 1, MPI_DOUBLE, (rank > 0)?(rank-1):(size-1), tagRight, MPI_COMM_WORLD, reqs);
      MPI_Isend(&f[nxn_loc-3], 1, MPI_DOUBLE, (rank+1)%size, tagRight, MPI_COMM_WORLD, reqs);
      // send left
      MPI_Irecv(&f[nxn_loc-1], 1, MPI_DOUBLE, (rank+1)%size, tagLeft, MPI_COMM_WORLD, reqs+1);
      MPI_Isend(&f[2], 1, MPI_DOUBLE, (rank > 0)?(rank-1):(size-1), tagLeft, MPI_COMM_WORLD, reqs+1);
    }
    MPI_Waitall(2, reqs, stats);
    commEnd = MPI_Wtime();
    time = commEnd - commStart;


    // here we finish the calculations

    // calculate first order derivative using central difference
    // here we need to correct value of the ghost cells!
    for (i=1; i<(nxn_loc-1); i++)
      dfdx[i] = (f[i+1] - f[i-1])/(2*dx);

    
    // Print f values
    if (rank==0){ // print only rank 0 for convenience
        printf("My rank %d of %d\n", rank, size );
        printf("Here are my values for f including ghost cells:\n");
        printf("[");
        for (i=0; i<nxn_loc-1; i++)
	        printf("%f, ", f[i]);
        printf("%f", f[nxn_loc-1]);
        printf("]\n");
        printf("Here are my values for dfdx:\n");
        printf("[");
        for (i=1; i<nxn_loc-2; i++)
	        printf("%f, ", dfdx[i]);
        printf("%f", dfdx[nxn_loc-2]);
        printf("]\n");
        printf("Communication time : %f\n",time);
    }

    free(f);
    free(dfdx);

    MPI_Finalize();
}






