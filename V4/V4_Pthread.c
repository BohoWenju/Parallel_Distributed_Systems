/*
     PARALLEL AND DISTRIBUTED SYSTEMS

  Main Code for V4 pthread of the exercise
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
#include <pthread.h>
#include "mmio.h"

pthread_mutex_t mutex=PTHREAD_MUTEX_INITIALIZER;
u_int32_t previous=0;
int thread_to_run=0;
int NTHREADS=1;
//  structures needed for the pthreading
typedef struct{
  u_int32_t* indices;
  u_int32_t* indptr;
  u_int32_t* p_indices;
  u_int32_t* p_indptr;
  u_int32_t* p_data;
  int nz;
  int id;
  int M;
}param_H;

typedef struct{
  u_int32_t* indptr;
  u_int32_t* indices;
  u_int32_t* data;
  u_int32_t* c3;
  int id;
  int M;
}param_c3;


//binary search
u_int32_t binary_match(u_int32_t* a,u_int32_t* b,u_int32_t l1,u_int32_t l2){
  u_int32_t matches=0;
  int flag=0;
  for(int i=0; i<l1; i++)
    if(flag==0)
      for(int j=0; j<l2; j+=2){
        if ( b[j]==a[i]){
          matches++;
          break;
        }
        else if (b[j]>a[i]){
          if ((b[j-1]==a[i]) && ((j-1)>=0)){
            matches++;
            break;
          }
          else
            break;
        }
        else
          if ((b[j+1]==a[i])&&(j+2==l2)){
            matches++;
            break;
          }
          else if ((b[j]<a[i])&&(j==(l2-1))){
            flag=1;
            break;
        }
      }
      else
        break;
    return matches;
}

//  funtion to calulate the A.*(A^2) of the given array
//  in CSC format with pthreads
void* Hadamard_squared_pthread(void* arg){
  param_H* p=(param_H*)arg;
  int nz= p->nz;
  int loopiter=0; //index to know in which repetition i am
  u_int32_t previous_t=0;
  u_int32_t* t_indices=(u_int32_t*)malloc(nz*sizeof(u_int32_t));
  u_int32_t* t_data=(u_int32_t*)malloc(nz*sizeof(u_int32_t));
  u_int32_t* t_indptr=(u_int32_t*)malloc((p->M/NTHREADS + 1)*sizeof(u_int32_t));
  t_indptr[loopiter]=0;
  //flag to show the thread when to stop calculating elements due to
  //the fact that no more elements are there for that thread to calculate
  int flag_H=1;
  while (thread_to_run < p->M){
    int i= p->id + loopiter*NTHREADS; //which column i am calculating
    if ( ((loopiter==0)||((p->id)!=(thread_to_run%NTHREADS))) && (flag_H==1)){
      if (i >= p->M)
        flag_H=0;
      else{
        int temp2=0;
        int temp=p->indptr[i+1]-p->indptr[i];
        for(int j=0; j<temp; j++){
          int row=p->indices[p->indptr[i]+j];
          int temp1=p->indptr[row+1]-p->indptr[row];
          u_int32_t sum=binary_match(p->indices+p->indptr[i],p->indices+p->indptr[row],temp,temp1);
          if (sum!=0){
            temp2++;
            t_data[previous_t]=sum;
            t_indices[previous_t]=row;
            previous_t++;
            sum=0;
          }
        }
        t_indptr[loopiter+1]=t_indptr[loopiter]+temp2;
        loopiter++;
      }
    }
    else if(((p->id)==(thread_to_run%NTHREADS)) && (thread_to_run < p->M)){
      int index=thread_to_run/NTHREADS ;
      int temp=t_indptr[index+1]-t_indptr[index];
      for(int j=0; j<temp; j++){
        p->p_indices[previous]=t_indices[t_indptr[index]+j];
        p->p_data[previous]=t_data[t_indptr[index]+j];
        previous++;
      }
      p->p_indptr[thread_to_run+1]=p->p_indptr[thread_to_run]+temp;
      thread_to_run++;
    }
  }
  free(t_indptr);
  free(t_indices);
  free(t_data);
  pthread_exit(NULL);
}

void* C_x_E(void* arg){
  param_c3* p=(param_c3*)arg;
  for(int j = (p->id); j < (p->M); j=j+NTHREADS){
    int temp = p->indptr[j+1] - p->indptr[j];
    for(int o = 0; o < temp; o++){
      int i=p->indices[p->indptr[j]+o];
      p->c3[i] += ((p->data[p->indptr[j]+o]));
    }
  }
  pthread_exit(NULL);
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

void* C_E_2(void* arg){
  param_c3* p=(param_c3*)arg;
  for(int j = (p->id); j < (p->M); j=j+NTHREADS){
    p->c3[j] = p->c3[j]/2;
    }
    pthread_exit(NULL);
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
  // Input number of threads
  if (argc<3){
    printf("No thread number was declared...exiting \n");
  	exit(1);
  }
  NTHREADS=(int)(atoi(argv[2]));
  printf("NTHREADS IS %d \n",NTHREADS);
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

   // ***********************************************
  //
  // TIME FOR CALCULATING TRIANGLES IN SERIAL
  clock_gettime(CLOCK_MONOTONIC,&ts_start);
  param_H param_H[NTHREADS];
  param_c3 param_c3[NTHREADS];
  pthread_t thread_arr[NTHREADS];

  u_int32_t* c3=(u_int32_t*)calloc(M,sizeof(u_int32_t));
  u_int32_t* indices_H=(u_int32_t *) malloc(nz * sizeof(u_int32_t));
  u_int32_t* dataH=(u_int32_t *) malloc(nz * sizeof(u_int32_t));
  u_int32_t* indptr_H=(u_int32_t*)malloc((M+1)*sizeof(u_int32_t));
  indptr_H[0]=0;
  for(int i=0; i<NTHREADS; i++){
    param_H[i].indices=indices;
    param_H[i].indptr=indptr;
    param_H[i].p_indices=indices_H;
    param_H[i].p_indptr=indptr_H;
    param_H[i].p_data=dataH;
    param_H[i].id=i;
    param_H[i].nz=nz;
    param_H[i].M=M;
    pthread_create(&thread_arr[i],NULL,Hadamard_squared_pthread,(void*)(param_H+i));
  }

  for(int i=0; i<NTHREADS; i++)
    pthread_join(thread_arr[i],NULL);
  indices_H=(u_int32_t*)realloc(indices_H,sizeof(u_int32_t)*(previous));
  dataH=(u_int32_t*)realloc(dataH,sizeof(u_int32_t*)*(previous));
  // counting triangles with the (A.*(A*A)*e/2) formula

  // size of array C,where C is the product
  // of the hadamard operation,
  // must be the same as the size
  // of array A,again that's just in our case
  for(int i=0; i<NTHREADS; i++){
    param_c3[i].indptr=indptr_H;
    param_c3[i].data=dataH;
    param_c3[i].indices=indices_H;
    param_c3[i].c3=c3;
    param_c3[i].id=i;
    param_c3[i].M=M;
    pthread_create(&thread_arr[i],NULL,C_x_E,(void*)(param_c3+i));
  }
  for(int i=0; i<NTHREADS; i++)
   pthread_join(thread_arr[i],NULL);

  for(int i=0; i<NTHREADS; i++){
    pthread_create(&thread_arr[i],NULL,C_E_2,(void*)(param_c3+i));
  }
  for(int i=0; i<NTHREADS; i++)
    pthread_join(thread_arr[i],NULL);

  clock_gettime(CLOCK_MONOTONIC,&ts_end);
  printf(" \n Time for parallel(pthreads) imlpementation  : %lf seconds \n",(double)ts_end.tv_sec +(double)ts_end.tv_nsec*(0.000000001) - (double)ts_start.tv_sec-(double)ts_start.tv_nsec*(0.000000001) );
  // ***********************************************
  free(indices_H);
  free(dataH);
  free(indptr_H);
  free(indices);
  free(indptr);
  return 0;
}
