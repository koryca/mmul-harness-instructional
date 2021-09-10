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

   for (int i=0; i<n; i+=block_size){
      for (int j=0; j<n; j+=block_size){
         //copy to local
         memcpy((void *)Clocal, (const void *)C, sizeof(double)*block_size*block_size);
         for(int k=0; k<n; k+=block_size){
            memcpy((void *)Alocal, (const void *)A, sizeof(double)*block_size*block_size);
            memcpy((void *)Blocal, (const void *)B, sizeof(double)*block_size*block_size);
            for (int p=i; p<i+block_size; p++){
               for (int q=j; q<j+block_size; q++){
                  double * temp = new double[block_size * block_size];
                  for(int m=k; m<k+block_size; m++){
                     // C[i,j] += A[i,k] * B[k,j]
                     temp[p*block_size+q*block_size]  += Alocal[p*block_size+m] * Blocal[m+q*block_size];
                  }
                  Clocal[p+q*block_size] = temp[p*block_size+q*block_size];
               }
            }
            memcpy((void *)C, (const void *)Clocal, sizeof(double)*block_size*block_size);
         }
      }
   }
}
