#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int main(int argc, char *argv[]){

    if(argc > 1){
        omp_set_num_threads(atoi(argv[1]));
    }
    #pragma omp parallel
    {
        int id = omp_get_thread_num();
        printf("Hello World from Thread %d!\n",id);
    }

    return 0;
}