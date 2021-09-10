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
   std::vector<double> buf(4 * n * n);
   double * Clocal = buf.data() + 0;
   double * Alocal = Clocal + n * n;
   double * Blocal = Alocal + n * n;
   double * temp = Blocal + n * n;

   int nb = n/block_size;
   for (int i=0; i<n; i+=nb){
      for (int j=0; j<n; j+=nb){
         //copy to local
         // memcpy((void *)Clocal, (const void *)C, sizeof(double)*block_size*block_size);
         memcpy((void *)temp, (const void *)C, sizeof(double)*block_size*block_size);
         for(int k=0; k<n; k+=nb){
            memcpy((void *)Alocal, (const void *)A, sizeof(double)*block_size*block_size);
            memcpy((void *)Blocal, (const void *)B, sizeof(double)*block_size*block_size);
            // temp[i + j * n] += A[i + k * n] * B[k + j *n];
            for (int jj=j; jj<j+block_size; jj++){
               for (int ii=i; ii<i+block_size; ii++){
                  for(int kk=k; kk<k+block_size; kk++){
                     // C[i,j] += A[i,k] * B[k,j]
                     C[ii + jj * block_size] += A[ii + kk * block_size] * B[kk + jj * block_size];
                     temp[ii + jj * block_size] += Alocal[ii + kk * block_size] * Blocal[kk + jj * block_size];
                     // std::cout << "After: " << Clocal[ii + jj] << " " << Alocal[ii + kk] << " " << Blocal[kk + jj] << std::endl;
                  }
               }
            }
            // memcpy((void *)C, (const void *)Clocal, sizeof(double)*block_size*block_size);
         }
      }
   }
   // std::cout << *A << " " << *B << " " << *C << std::endl;
   std::cout << *A << " " << *B << " " << *C << " " << *temp << std::endl;
}
