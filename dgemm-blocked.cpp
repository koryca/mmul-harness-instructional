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
        
   double * Alocal = (double*) malloc(n * n * sizeof(double));
   double * Blocal = (double*) malloc(n * n * sizeof(double));
   double * Clocal = (double*) malloc(n * n * sizeof(double));

   for (int i=0; i<n; i+=block_size){
      for (int j=0; j<n; j+=block_size){ 
         //copy C
         // for(int ic = i; ic < i + block_size; ic++){
         for(int jc = j; jc < j + block_size; jc++){
               memcpy(&Clocal[i + jc * block_size], &C[i + jc * n], sizeof(double)*block_size);
               // std::cout << "Clocal at copy: " << Clocal[ic + j * block_size] << " " << Clocal[ic + j * block_size + 1]
               //          << " C at copy: " << C[ic * n + j] << " " << C[ic * n + j + 1]
               //          << " C[" << ic << "][" << j << "]"<< std::endl;
               std::cout << "C: " << Clocal[i + jc * block_size] << " " << i + jc * block_size << " "
                        << Clocal[i + jc * block_size + 1] << " " << i + jc * block_size + 1 << std::endl;
            // }
         }
         for(int k=0; k<n; k+=block_size){ 
            //copy A
            // for(int ia = i; ia < i + block_size; ia++){
            for(int ka = k; ka < k + block_size; ka++){
               memcpy(&Alocal[i + ka * block_size], &A[i + ka * n], sizeof(double)*block_size);
               // std::cout << "Alocal at copy: " << Alocal[ia + k * block_size] << " " << Alocal[ia + k * block_size + 1]
               //        << " A at copy: " << A[ia * n + k] << " " << A[ia * n + k + 1]
               //        << " A[" << ia << "][" << k << "]" << std::endl;
               std::cout << "A: " << Alocal[i + ka * block_size] << " " << i + ka * block_size << " "
                        << Alocal[i + ka * block_size + 1] << " " << i + ka * block_size + 1 << std::endl;
               // }
            }
            //copy B
            // for(int kb = k; kb < k + block_size; kb++){
            for(int jb = j; jb < j + block_size; jb++){
               memcpy(&Blocal[k * block_size + jb * block_size], &B[k + jb * n], sizeof(double)*block_size);
               // std::cout << "Blocal at copy: " << Blocal[kb + j * block_size] << " " << Blocal[kb + j * block_size + 1]
               //        << " B at copy: " << B[kb * n + j] << " " << B[kb * n + j + 1]
               //        << " B[" << kb << "][" << j << "]"<< std::endl;
               std::cout <<"B: " << Blocal[k * block_size + jb * block_size] << " " << k * block_size + jb * block_size<< " "
                        << Blocal[k * block_size + jb * block_size + 1] << " " << k * block_size + jb * block_size + 1 << std::endl;
               // }
            }  
            for (int ii=i; ii<i+block_size; ii++){
               for (int jj=j; jj<j+block_size; jj++){ 
                  double temp = 0.0;
                  for(int kk=k; kk<k+block_size; kk++){ 
                     temp += Alocal[ii * block_size + kk] * Blocal[kk + jj * block_size];
                     std::cout << "A[" << ii << "][" << kk << "] " << Alocal[ii + kk * block_size] << " " << ii + kk * block_size
                              << " B[" << kk << "][" << jj << "] " << Blocal[kk + jj * block_size]<< " " << kk + jj * block_size
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
               memcpy(&C[i + jjc * n], &Clocal[i + jjc * block_size], sizeof(double)*block_size);
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
