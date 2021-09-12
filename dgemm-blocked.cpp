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
         // for(int ic = i; ic < i + block_size; ic++){
         for(int jc = j; jc < j + block_size; jc++){
               memcpy(&Clocal[i + jc * n], &C[i + jc * n], sizeof(double)*block_size);
               std::cout << "Clocal at copy: " << Clocal[i + jc * n] 
                        << " C at copy: " << C[i + jc * n] 
                        << " C[" << i << "][" << jc << "]"<< std::endl;
            // }
         }
         for(int k=0; k<n; k+=block_size){ 
            //copy A
            // for(int ia = i; ia < i + block_size; ia++){
            for(int ka = k; ka < k + block_size; ka++){
               memcpy(&Alocal[i + ka * block_size], &A[i * n + ka], sizeof(double)*block_size);
               std::cout << "Alocal at copy: " << Alocal[i * block_size + ka] 
                      << " A at copy: " << A[i * n + ka] 
                      << " A[" << i << "][" << ka << "]" << std::endl;
               // }
            }
            //copy B
            // for(int kb = k; kb < k + block_size; kb++){
            for(int jb = j; jb < j + block_size; jb++){
               memcpy(&Blocal[k + jb * block_size], &B[k + j * n], sizeof(double)*block_size);
               std::cout << "Blocal at copy: " << Blocal[k + j * block_size] 
                      << " B at copy: " << B[k + j * n] 
                      << " B[" << k << "][" << j << "]"<< std::endl;
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
         // for(int iic = i; iic < i + block_size; iic++){
         for(int jjc = j; jjc < j + block_size; jjc++){
               memcpy(&C[i + jjc * n], &Clocal[i + jjc * n], sizeof(double)*block_size);
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
