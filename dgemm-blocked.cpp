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
   // std::vector<double> buffer(3 * n * n);
   // double * Alocal = buffer.data() + 0;
   // double * Blocal = Alocal + n * n;
   // double * Clocal = Blocal + n * n;
        
   double * Alocal = (double*) malloc(block_size * sizeof(double));
   double * Blocal = (double*) malloc(block_size * sizeof(double));
   double * Clocal = (double*) malloc(block_size * sizeof(double));

   for (int i=0; i<n; i+=block_size){
      for (int j=0; j<n; j+=block_size){ 
         //copy C
         for(int ic = i; ic < i + block_size; ic++){
         // for(int jc = j; jc < j + block_size; jc++){
               memcpy(&Clocal[ic + j * n], &C[ic + j * n], sizeof(double)*block_size);
               std::cout << "Clocal at copy: " << Clocal[ic + j * n] 
                        << " C at copy: " << C[ic + j * n] 
                        << " C[" << ic << "][" << j << "]"<< std::endl;
            // }
         }
         for(int k=0; k<n; k+=block_size){ 
            //copy A
            for(int ia = i; ia < i + block_size; ia++){
            // for(int ka = k; ka < k + block_size; ka++){
               memcpy(&Alocal[ia + k * block_size], &A[ia * n + k], sizeof(double)*block_size);
               std::cout << "Alocal at copy: " << Alocal[ia + k * block_size] << " " << Alocal[ia + k * block_size + 1]
                      << " A at copy: " << A[ia * n + k] << A[ia * n + k + 1]
                      << " A[" << ia << "][" << k << "]" << std::endl;
               // }
            }
            //copy B
            for(int kb = k; kb < k + block_size; kb++){
            // for(int jb = j; jb < j + block_size; jb++){
               memcpy(&Blocal[kb + j * block_size], &B[kb * n + j], sizeof(double)*block_size);
               std::cout << "Blocal at copy: " << Blocal[kb + j * block_size] << " " << Blocal[kb + j * block_size + 1]
                      << " B at copy: " << B[kb * n + j] << " " << B[kb * n + j + 1]
                      << " B[" << kb << "][" << j << "]"<< std::endl;
               // }
            }  
            for (int ii=i; ii<i+block_size; ii++){
               for (int jj=j; jj<j+block_size; jj++){ 
                  double temp = 0.0;
                  for(int kk=k; kk<k+block_size; kk++){ 
                     temp += Alocal[ii + kk * block_size] * Blocal[kk + jj * block_size];
                     std::cout << "A[" << ii << "][" << kk << "] " << Alocal[ii + kk * block_size] 
                              << " B[" << kk << "][" << jj << "] " << Blocal[kk + jj * block_size]
                              << " ->C[" << ii << "][" << jj << "] += " << Alocal[ii + kk * block_size] 
                                               << " * " << Blocal[kk + jj * block_size] << std::endl;
                  }
                  Clocal[ii + jj * block_size] += temp;
               }
            }
         }
         //write back to C
         for(int iic = i; iic < i + block_size; iic++){
         // for(int jjc = j; jjc < j + block_size; jjc++){
               memcpy(&C[iic + j * n], &Clocal[iic + j * n], sizeof(double)*block_size);
            // std::cout << "Clocal at write: " << Clocal[iic + j * block_size] 
            //           << "C at write: " << C[iic + j * block_size] << std::endl;      
            // }
         }
      }
   }
   // free(Alocal);
   // free(Blocal);
   // free(Clocal);

   // delete[] Alocal;
   // delete[] Blocal;
   // delete[] Clocal;
}
