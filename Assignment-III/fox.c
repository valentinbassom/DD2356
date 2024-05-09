#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <mpi.h>

#define M 120

double **Allocate_Matrix(int m){
    double **matrix = (double **)malloc(m * sizeof(double *));
    double *data = (double *)malloc(m * m * sizeof(double));

    for (int i = 0; i < m; i++) {
        matrix[i] = data + i*m;
    }
    return matrix;
}

void Free_Matrix(double **matrix){
    free(matrix[0]);
    free(matrix);
}

void Initialize_Matrix(double **matrix, int m, int coords[2]){
    for(int i=0; i < m; i++){
        for(int j=0; j < m; j++){
            matrix[i][j] = (coords[1]*M*m + coords[0]*m + i*M + j);
        }
    }
}

void Initialize_To_0(double **matrix, int m){
    for(int i=0; i < m; i++){
        for(int j=0; j < m; j++){
            matrix[i][j] = 0.0;
        }
    }
}

void Multiply_Matrices(double **A, double **B, double **C, int m){
    for(int i = 0 ; i < m ; i++){
        for(int k = 0 ; k < m ; k++){
            for(int j = 0 ; j < m ; j++){
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void print_matrix(double **matrix, int m){
    for(int i = 0; i < m; i++){
        printf("[ ");
        for(int j = 0; j < m; j++){
            if(matrix[i][j] < 0) printf(" %8.3f ", matrix[i][j]);
            else printf("  %8.4f ", matrix[i][j]);
        } printf(" ]\n");
    } printf("\n");
}

int main(int argc, char* argv[])
{
    int rank, sub_rank, size, provided, grid_size, m, rank_up, rank_down;
    int coords[2];
    double time;
    MPI_Comm CartComm, RowComm;
    MPI_Status stat;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double sqrtN = sqrt(size);
    if(sqrtN != floor(sqrtN)){
        if(rank==0){
            printf("The requested amount of processes (%d) is not a square! Decomposition in subblocks not feasible. Exiting...\n", size);
        }
        MPI_Finalize();
        exit(0);
    }
    grid_size = (int) sqrtN;
    double M_sqrtN = (double) M / sqrtN;
    if(M_sqrtN != floor(M_sqrtN)){
        if(rank==0){
            printf("The requested amount of processes (%d) and associated grid size (%d) do not match the matrix size (%d)! Exiting...\n", size, grid_size, M);
        }
        MPI_Finalize();
        exit(0);
    }
    m = (int) M_sqrtN;

    double **A_loc, **A_temp, **B_even, **B_odd, **C;
    A_loc  = Allocate_Matrix(m);
    A_temp = Allocate_Matrix(m);
    B_even = Allocate_Matrix(m);
    B_odd  = Allocate_Matrix(m);
    C      = Allocate_Matrix(m);

    Initialize_To_0(C, m);

    //Create Cartesian communicator
    int dims[2] = {grid_size,grid_size};
    int periods[2] = {1,1};
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &CartComm);
    MPI_Comm_rank(CartComm, &rank);
    MPI_Cart_coords(CartComm, rank, 2, coords);
    MPI_Cart_shift(CartComm, 1, 1, &rank_up, &rank_down);
    if(rank==0) printf("Matrix size: %d*%d, Number of processes: %d, Grid size: %d*%d, Block size: %d*%d\n",M,M,size,grid_size,grid_size,m,m);

    int kept_dims[2] = {1,0};
    MPI_Cart_sub(CartComm, kept_dims, &RowComm);
    MPI_Comm_rank(RowComm, &sub_rank);
    if(sub_rank != coords[0]) printf("Error : sub_rank (%d) not equal to column index (%d)\n",sub_rank,coords[0]);

    Initialize_Matrix(A_loc, m, coords);
    Initialize_Matrix(B_even, m, coords);


    time = MPI_Wtime();
    for(int ite=0; ite < grid_size; ite++){
        //Broadcast diagonal+ite
        if(coords[0] == (coords[1]+ite)%grid_size){
            memcpy(A_temp[0], A_loc[0], m*m*sizeof(double));
        }
        MPI_Bcast(A_temp[0], m*m, MPI_DOUBLE, (coords[1]+ite)%grid_size, RowComm);
        
        //Matrix multiply
        if(ite%2==0){
            Multiply_Matrices(A_temp,B_even,C,m);
        } else {
            Multiply_Matrices(A_temp,B_odd,C,m);
        }

        //Roll B
        // /!\ For practical applications it is probably better to remove the first condition (ite != grid_size-1):
        // That way, at the end of execution B is stored as it was previously (at the beginning of execution), which is not the case now.
        // If you do so, pay attention also to odd and even iterations, the original B will be in B_odd instead of B_even if grid_size is odd.
        if(ite != grid_size-1){
            if(ite%2==0){
                MPI_Sendrecv(B_even[0], m*m, MPI_DOUBLE, rank_up, 1, B_odd[0], m*m, MPI_DOUBLE, rank_down, 1, CartComm, &stat);
            } else {
                MPI_Sendrecv(B_odd[0], m*m, MPI_DOUBLE, rank_up, 1, B_even[0], m*m, MPI_DOUBLE, rank_down, 1, CartComm, &stat);
            }
        }
    }

    time = MPI_Wtime() - time;
    //Remove if statement and uncomment the sleep operation if you want to print the full results.
    if(coords[0]==0 && coords[1] == 0){
        //sleep(coords[1] + coords[0]*grid_size);
        printf("The results for the block C_{%d,%d} are:\n",coords[1],coords[0]);
        print_matrix(C,m);
        printf("Execution time: %f seconds\n",time);
    }

    Free_Matrix(A_loc);
    Free_Matrix(A_temp);
    Free_Matrix(B_even);
    Free_Matrix(B_odd);
    Free_Matrix(C);
    
    MPI_Finalize();
    return 0;
}