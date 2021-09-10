const char* dgemm_desc = "Basic implementation, three-loop dgemm.";

/*
 * This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are n-by-n matrices stored in column-major format.
 * On exit, A and B maintain their input values.
 */
void square_dgemm(int n, double* A, double* B, double* C) 
{
   // This implementation is the triply nested loop as discussed in class
   for (int i=0; i<n; i++){
      for (int j=0; j<n; j++){
         for(int k=0; k<n; k++){
            // C[i,j] += A[i,k] * B[k,j]
            C[i*n+j] += A[i*n+k] * B[k*n+j];
         }
      }
   }
}
