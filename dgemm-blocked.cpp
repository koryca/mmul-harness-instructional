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
   std::vector<double> buffer(4 * n * n);
   double * Alocal = buffer.data() + 0;
   double * Blocal = Alocal + n * n;
   double * Clocal = Blocal + n * n;
        
   for (int i=0; i<n; i+=block_size){
      for (int j=0; j<n; j+=block_size){ 
         //copy C
         for(int ic = i; ic < i + block_size; ic++){
            for(int jc = j; jc < j + block_size; jc++){
               memcpy(&Clocal[ic + jc * n], &C[ic + jc * n], sizeof(double)*block_size*block_size);
            }
         }
         for(int k=0; k<n; k+=block_size){ 
            //copy A
            for(int ia = i; ia < i + block_size; ia++){
               for(int ka = k; ka < k + block_size; ka++){
                  memcpy(&Alocal[ia + ka * n], &A[ia + ka * n], sizeof(double)*block_size*block_size);
               }
            }
            //copy B
            for(int kb = k; kb < k + block_size; kb++){
               for(int jb = j; jb < j + block_size; jb++){
                  memcpy(&Blocal[kb + jb * n], &B[kb + jb* n], sizeof(double)*block_size*block_size);
               }
            }
            for (int ii=i; ii<i+block_size; ii++){
               for (int jj=j; jj<j+block_size; jj++){ 
                  for(int kk=k; kk<k+block_size; kk++){ 
                     Clocal[ii + jj * n] += Alocal[ii + kk * n] * Blocal[kk + jj * n];
                  }
               }
            }
         }
         //write back to C
         for(int iic = i; iic < i + block_size; iic++){
            for(int jjc = j; jjc < j + block_size; jjc++){
               memcpy(&C[iic + jjc * n], &Clocal[iic + jjc * n], sizeof(double)*block_size*block_size);
               
            }
         }
      }
   }
}
