/*
     PARALLEL AND DISTRIBUTED SYSTEMS

  Main Code for V3 openmp of the exercise
  Aikaterini Prokou/Eleftherios Mourelatos


  Calculates the number of triangles in a
  given undirected graph.The graphs were
  taken from Matrix market:
  (See : https://math.nist.gov/MatrixMarket/mmio-c.html)
  and have been formatted according to the
  CSC Format



*/



#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <omp.h>
#include "mmio.h"




//  function to revert a matric from a COO format
//  to a CSC format
//  taken from: https://github.com/AUTh-csal/pds-codebase/blob/main/coo2csc/coo2csc.c
void coo2csc(
  u_int32_t       * const row,       /*!< CSC row start indices */
  u_int32_t       * const col,       /*!< CSC column indices */
  u_int32_t const * const row_coo,   /*!< COO row indices */
  u_int32_t const * const col_coo,   /*!< COO column indices */
  u_int32_t const         nnz,       /*!< Number of nonzero elements */
  u_int32_t const         n,         /*!< Number of rows/columns */
  u_int32_t const         isOneBased /*!< Whether COO is 0- or 1-based */
) {
  // ----- cannot assume that input is already 0!
  for (u_int32_t l = 0; l < n+1; l++) col[l] = 0;


  // ----- find the correct column sizes
  for (u_int32_t l = 0; l < nnz; l++)
    col[col_coo[l] - isOneBased]++;

  // ----- cumulative sum
  for (u_int32_t i = 0, cumsum = 0; i < n; i++) {
    u_int32_t temp = col[i];
    col[i] = cumsum;
    cumsum += temp;
  }
  col[n] = nnz;
  // ----- copy the row indices to the correct place
  for (u_int32_t l = 0; l < nnz; l++) {
    u_int32_t col_l;
    col_l = col_coo[l] - isOneBased;

    u_int32_t dst = col[col_l];
    row[dst] = row_coo[l] - isOneBased;

    col[col_l]++;
  }
  // ----- revert the column pointers
  for (u_int32_t i = 0, last = 0; i < n; i++) {
    u_int32_t temp = col[i];
    col[i] = last;
    last = temp;
  }
}

void triangles_parallel_openMP_count(u_int32_t* indices,u_int32_t* indptr,int M){
  u_int32_t triangles=0;
  u_int32_t Omp_array[M];
  #pragma omp parallel for
  for(int j=0;j<M; j++)
    Omp_array[j]=0;
  #pragma omp parallel for
  for(int j=0; j<M; j++){
    int temp=indptr[j+1]-indptr[j];
    //i is the index of the row
    //for the a[i][j] element
    for(int o=0; o<temp; o++){
      int i=indices[indptr[j]+o];
      if (i<=j)
        continue;
      if ((temp>=2)){
      for(int k=0; k<temp; k++){
        int z=indices[indptr[j]+k];
        if((z!=i)&&(z!=j)){
          int temp1=indptr[z+1]-indptr[z];
          for(int w=0; w<temp1; w++)
          if (indices[indptr[z]+w]==i)
            Omp_array[j]++;
     }
    }
   }
  }
 }
 for(int j=0; j<M; j++)
    triangles+=Omp_array[j];
  printf("The number of triangles is %d \n",triangles);
}

int main(int argc, char *argv[]){
  // code to get the .mtx file
  // and take the matrix in COO format
  // the code is the same as the
  // example_read.c
  // VARIABLES NEEDED FOR EXAMPLE_READ
  int ret_code;
  MM_typecode matcode;
  FILE *f;
  int M, N, nz;
  int *I, *J;
  double *val;

  // VARIABLES NEEDED FOR THE EXERCISE
  u_int32_t *indices,*indptr;
  struct timespec ts_start,ts_end;
  if (argc < 2){
  	fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
  	exit(1);
  }
  else
    if ((f = fopen(argv[1], "r")) == NULL)
              exit(1);

  if (mm_read_banner(f, &matcode) != 0){
          printf("Could not process Matrix Market banner.\n");
          exit(1);
  }
  if(argc<3){
    printf("No thread number was declared...exiting");
    exit(1);
  }
  int NTHREADS=(int)(atoi(argv[2]));
  printf("NTHREADS IS %d \n",NTHREADS);
  omp_set_num_threads(NTHREADS);

  /*  This is how one can screen matrix types if their application */
  /*  only supports a subset of the Matrix Market data types.      */

  if (mm_is_complex(matcode) && mm_is_matrix(matcode) && mm_is_sparse(matcode) ){
      printf("Sorry, this application does not support ");
      printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
      exit(1);
  }

  /* find out size of sparse matrix .... */

  if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
      exit(1);


  /* reseve memory for matrices */

  I = (int *) malloc(nz * sizeof(int));
  J = (int *) malloc(nz * sizeof(int));
  val = (double *) malloc(nz * sizeof(double));


  /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
  /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
  /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

  for (int i=0; i<nz; i++){
    fscanf(f, "%d %d\n", &I[i], &J[i]);
    I[i]--;  /* adjust from 1-based to 0-based */
    J[i]--;
  }
  // val has all the nz elements
  // I has all the rows
  // J has all the columns
  if (f !=stdin) fclose(f);

  /************************/
  /* now write out matrix */
  /************************/
  // Matrix must be squared
  if (M!=N){
    printf("The given matrix is not squared");
    exit(1);
  }
  //u_int32_t array_of_cols[M];
  indices = (u_int32_t *) malloc(nz * sizeof(u_int32_t));
  // indptr has M+1 size due to its definition
  indptr  = (u_int32_t *) malloc((M+1) * sizeof(u_int32_t));

  // code above adjusts the matrix given
  // from 1 base to 0 based
  // function to covert the matrix from coo format
  // to csc format
  coo2csc(indices,indptr,I,J,nz,N,0);
  mm_write_banner(stdout, matcode);
  mm_write_mtx_crd_size(stdout, M, N, nz);

  // ***********************************************
  //
  // TIME FOR CALCULATING TRIANGLES IN PARALLEL
  // with OpenMP
  clock_gettime(CLOCK_MONOTONIC,&ts_start);
  // algorithm to count the triangles
  // in a csc formatted Matrix
  triangles_parallel_openMP_count(indices,indptr,M);
  clock_gettime(CLOCK_MONOTONIC,&ts_end);
  printf("Time for parallel(openMP) imlpementation in secs : %lf \n",( (double)ts_end.tv_sec+(double)ts_end.tv_nsec*(0.000000001) - (double)ts_start.tv_sec-(double)ts_start.tv_nsec*(0.000000001)));

  // ***********************************************
  free(indices);
  free(indptr);
  return 0;
}
