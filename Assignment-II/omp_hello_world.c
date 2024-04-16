#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

int main(int argc, char *argv[]){

    if(argc > 1 && strcmp(argv[1],"-n_threads") == 0){
        omp_set_num_threads(atoi(argv[2]));
    }
    #pragma omp parallel
    {
        int id = omp_get_thread_num();
        printf("Hello World from Thread %d!\n",id);
    }

    return 0;
}