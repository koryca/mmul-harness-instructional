#include <math>
#include <vector>
#include <cstring>
const char* dgemm_desc = "Blocked dgemm.";

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are n-by-n matrices stored in column-major format.
 * On exit, A and B maintain their input values. */
void square_dgemm_blocked(int n, int block_size, double* A, double* B, double* C) 
{
   std::vector<double> buf(3 * n * n);
   double * Clocal = buf.data() + 0;
   double * Alocal = Clocal + n * n;
   double * Blocal = Alocal + n * n;

   for (int p=0; p<n; p+=block_size){
      for (int q=0; q<n; q+=block_size){
         //copy to local
         memcpy((void *)Clocal, (const void *)C, sizeof(double)*block_size*block_size);
         // for(int k=0; k<n; k+=block_size){
            memcpy((void *)Alocal, (const void *)A, sizeof(double)*block_size*block_size);
            memcpy((void *)Blocal, (const void *)B, sizeof(double)*block_size*block_size);
            for (int i=0; i<n; i++){
               for (int j=p; j<p+block_size; j++){
                  for(int k=0; k<q+block_size; k++){
                     // C[i,j] += A[i,k] * B[k,j]
                     Clocal[i+j*block_size] += Alocal[i+k*block_size] * Blocal[k+j*block_size];
                  }
               }
            }
            memcpy((void *)C, (const void *)Clocal, sizeof(double)*block_size*block_size);
         // }
      }
   }
}
