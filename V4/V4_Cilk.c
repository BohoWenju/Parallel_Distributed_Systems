/*
     PARALLEL AND DISTRIBUTED SYSTEMS

  Main Code for V4 cilk of the exercise
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
#include <cilk/cilk.h>
#include <cilk/reducer_opadd.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include "mmio.h"



//binary search
u_int32_t binary_match(u_int32_t* a,u_int32_t* b,u_int32_t l1,u_int32_t l2){
  int matches=0;
  for(int i=0; i<l1; i++)
    for(int j=0; j<l2; j+=2){
      if ( b[j]==a[i])
        matches++;
      else if (b[j]>a[i]){
        if ((b[j-1]==a[i]) && ((j-1)>=0))
          matches++;
        else
          break;
      }
      else if ((b[j]<a[i]) &&(j+2==l2))
        if(b[j+1]==a[i])
          matches++;
      }
    return matches;
}
//  funtion to calulate the A.*(A^2) of the given array
//  in CSC format
u_int32_t parallel_cilk_Hadamard_squared(u_int32_t* indices, u_int32_t* indptr,int M,
u_int32_t* p_indices,u_int32_t* p_indptr,u_int32_t* p_data){
  CILK_C_REDUCER_OPADD(sum,int, 0);    // Declare the reducer
  CILK_C_REGISTER_REDUCER(sum);
  p_indptr[0]=0;
  int temp2=0;
  int previous=0;
  cilk_for (int i = 0; i<M; ++i){
    u_int32_t temp=indptr[i+1]-indptr[i];
    for (u_int32_t j=0; j<temp; j++){
      u_int32_t row=indices[indptr[i]+j];
      u_int32_t temp1=indptr[row+1]-indptr[row];
      REDUCER_VIEW(sum)=binary_match(indices+indptr[i],indices+indptr[row],temp,temp1);
      if(REDUCER_VIEW(sum)>0){
        temp2++;
        p_data[previous]=REDUCER_VIEW(sum);
        p_indices[previous]=row;
        previous++;
        REDUCER_VIEW(sum)=0;
      }

  }
  p_indptr[i+1]=temp2;
}
 CILK_C_UNREGISTER_REDUCER(sum);


  p_indices=(u_int32_t*)realloc(p_indices,sizeof(u_int32_t)*(previous+1));
  p_data=(u_int32_t*)realloc(p_data,sizeof(u_int32_t)*(previous+1));
  return (previous);
}

void sort(int* a,int l){
  for(int i=0;i<l; i++)
    for(int j=0; j<l; j++)
      if ( (a[i]<a[j]) && (i>j) ){
        int temp=a[i];
        a[i]=a[j];
        a[j]=temp;
      }

}

void c_cilk_arr(u_int32_t * c3,int M,u_int32_t* indptr_H,u_int32_t* indices_H,u_int32_t* dataH){
CILK_C_REDUCER_OPADD(c3,uint, 0);
  CILK_C_REGISTER_REDUCER(c3);
  cilk_for(int j = 0; j < M; j++){
    c[j]=0;
    int temp = indptr_H[j+1] - indptr_H[j];
    for(int o = 0; o < temp; o++){
      int i = indices_H[indptr_H[j]+o];
      REDUCER_VIEW(c3)+= dataH[indptr_H[j]+o];
    }
    c[j]=REDUCER_VIEW(c3)/2;
    REDUCER_VIEW(c3)=0;
  }
  CILK_C_UNREGISTER_REDUCER(c3);
}
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
  //arrays for this exercise are symmetrical
  int* Isym = (int *)malloc((2*nz)*sizeof(int));
  int* Jsym =(int *)malloc((2*nz)*sizeof(int));


  for(int i=0; i<nz; i++){
    Isym[i]=I[i];
    Isym[nz+i]=J[i];
    Jsym[i]=J[i];
    Jsym[nz+i]=I[i];
  }


  nz= nz*2;
  //u_int32_t array_of_cols[M];
  indices = (u_int32_t *) malloc(nz * sizeof(u_int32_t));
  // indptr has M+1 size due to its definition
  indptr  = (u_int32_t *) malloc((M+1) * sizeof(u_int32_t));

  // code above adjusts the matrix given
  // from 1 base to 0 based
  // function to covert the matrix from coo format
  // to csc format
  coo2csc(indices,indptr,Isym,Jsym,nz,N,0);

  mm_write_banner(stdout, matcode);
  mm_write_mtx_crd_size(stdout, M, N, nz);

  for(int i=0; i<M; i++){
    u_int32_t l=indptr[i+1]-indptr[i];
    sort(indices+indptr[i],l);
  }
/*
  printf("indices are :");
  for (int i=0; i<nz; i++)
    printf("%d ",indices[i]);
  printf("\n");
  printf("indptr are :");
  for (int i=0; i<(M+1); i++)
    printf("%d ",indptr[i]);
  printf("\n");
*/
  /*printing the coo format of the matrix

  for (int i=0; i<nz; i++)
    fprintf(stdout, "%d %d %20.19g\n", I[i]+1, J[i]+1);

  */
  // checking if the arary of indices is the same as I
  //u_int32_t pass=1;
  //for (u_int32_t i = 0; i < nz; i++) {
      //pass &= (indices[i] == Isym[i]);
   //}
  //if (pass)
    //printf(" all good \n");

  //*********************************************
  int numWorkers = __cilkrts_get_nworkers();
  printf("The number of workers to %d.\n",numWorkers);

   // ***********************************************
  //
  // TIME FOR CALCULATING TRIANGLES IN SERIAL
  clock_gettime(CLOCK_MONOTONIC,&ts_start);

  // counting triangles with the (A!*(A*A)*e/2) formula

  // size of array C,where C is the product
  // of the hadamard operation,
  // must be the same as the size
  // of array A,again that's just in our case
  u_int32_t* indices_H=(u_int32_t *) calloc(nz , sizeof(u_int32_t));
  u_int32_t* dataH=(u_int32_t *) calloc(nz , sizeof(u_int32_t));
  u_int32_t* indptr_H=(u_int32_t*)malloc((M+1)*sizeof(u_int32_t));
  u_int32_t new_size=parallel_cilk_Hadamard_squared(indices,indptr,M,indices_H,indptr_H,dataH);

  u_int32_t* c=(u_int32_t*)calloc(M,sizeof(u_int32_t));
  c_cilk_arr(c,M,indptr_H,indices_H,dataH);

  clock_gettime(CLOCK_MONOTONIC,&ts_end);
  printf("Time for cilk imlpementation in secs : %lf secs\n", (double)ts_end.tv_sec +(double)ts_end.tv_nsec*(0.000000001) - (double)ts_start.tv_sec-(double)ts_start.tv_nsec*(0.000000001) );

/*
for(int i=0; i<M; i++)
    printf("%d %d \n", i, c[i]);
*/
  // ***********************************************
  free(indices_H);
  free(dataH);
  free(indptr_H);
  free(indices);
  free(indptr);
  return 0;
}
