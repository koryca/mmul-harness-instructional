#include <algorithm>
#include <vector>
#include <cstring>
const char* dgemm_desc = "Blocked dgemm.";

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are n-by-n matrices stored in column-major format.
 * On exit, A and B maintain their input values. */
void square_dgemm_blocked(int n, int block_size, double* A, double* B, double* C) 
{

   for (int i=0; i<n; i+=block_size){
      for (int j=0; j<n; j+=block_size){
         //copy to local
         std::vector<double> buf(12 * block_size * block_size);
         double * Clocal = buf.data() + block_size * block_size;
         // double Clocal[1024];
         // std::copy(C, C + block_size, Clocal);
         memcpy((void *)Clocal, (const void *)C, sizeof(double)*block_size*block_size);
         for(int k=0; k<n; k+=block_size){
            double * Alocal = Clocal + block_size * block_size;
            double * Blocal = Alocal + block_size * block_size;
            // long double * Alocal = new long double[block_size];
            // long double * Blocal = new long double[block_size];
            // std::copy(A, A + block_size, Alocal);
            // std::copy(B, B + block_size, Blocal);
            memcpy((void *)Alocal, (const void *)A, sizeof(double)*block_size*block_size);
            memcpy((void *)Blocal, (const void *)B, sizeof(double)*block_size*block_size);
            for (int p=i; p<i+block_size; p++){
               for (int q=j; q<j+block_size; q++){
                  for(int m=k; m<k+block_size; m++){
                     // C[i,j] += A[i,k] * B[k,j]
                     Clocal[p+q*block_size] += Alocal[p+m*block_size] * Blocal[m+q*block_size];
                  }
               }
            }
            // std::copy(Clocal, Clocal + block_size, C);
            memcpy((void *)C, (const void *)Clocal, sizeof(double)*block_size*block_size);
         }
      }
   }
}
