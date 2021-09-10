#include <vector>
#include <iostream>
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
   // double * temp = Blocal + n * n;

   for (int i=0; i<n; i+=block_size){
      for (int j=0; j<n; j+=block_size){
         //copy to local
         memcpy((void *)Clocal, (const void *)C, sizeof(double)*block_size*block_size);
         for(int k=0; k<n; k+=block_size){
            memcpy((void *)Alocal, (const void *)A, sizeof(double)*block_size*block_size);
            memcpy((void *)Blocal, (const void *)B, sizeof(double)*block_size*block_size);
            for (int ii=i; ii<i+block_size; ii++){
               for (int jj=j; jj<j+block_size; jj++){
                  for(int kk=k; kk<k+block_size; kk++){
                     // C[i,j] += A[i,k] * B[k,j]
                     Clocal[ii + jj * block_size] += Alocal[ii + kk * block_size] * Blocal[kk + jj * block_size];
                  }
               }
            }
            memcpy((void *)C, (const void *)Clocal, sizeof(double)*block_size*block_size);
         }
      }
      std::cout << *A << " " << *B << " " << *C << std::endl;
   }
}
