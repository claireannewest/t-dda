//Determine whether 'target' is located in 'arr'
int find_neighbor(int *arr, int N, int target);

//BLAS matrix-vector multiplication
extern void sgemv_(char *TRANS, int *M, int *N, float *ALPHA, float *A, int *LDA, float *X, int *INCX, float *BETA, float *Y, int *INCY);

//BLAS matrix-matrix multiplication
extern void sgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K, float *ALPHA, float *A, int *LDA, float *B, int *LDB, float *BETA, float *C, int *LDC);

//LAPACK linear system solver
extern void dgesv_(int *n, int *nrhs, double *a, int *lda,
                   int *ipiv, double *b, int *ldb, int *info);
