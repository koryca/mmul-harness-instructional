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
   std::vector<double> buf(3 * block_size * block_size);
   double * Clocal = buf.data() + 0;
   double * Alocal = Clocal + block_size * block_size;
   double * Blocal = Alocal + block_size * block_size;
   
   // test with block_size 2 BxB=2x2 MM size 64x64
   // i from 0 to n, increment by block_size, so n/block_size = 64/2 = 32 times of iteration
   for (int i=0; i<n; i+=block_size){
      for (int j=0; j<n; j+=block_size){ //same as i
         //original C, column
         for(int ci = 0; ci<n; ci++){
            for(int cj = 0; cj<n; cj++){
               std::cout << "C[" << cj <<"][" << ci << "] = " << C[ci+cj*n] 
               << " col index: " << ci+cj*n << std::endl;
            }
         }
         //copy C
         for(int ic = i; ic < i + block_size; ic++){
            for(int jc = j; jc < j + block_size; jc++){
               memcpy(&Clocal[ic + jc * n], &C[ic + jc * n], sizeof(double)*block_size*block_size);
               std::cout << "local: " << Clocal[ic + jc * n] << " " 
                        << "C: " << C[ic + jc * n] 
                        << " i[" << jc << "], j[" << ic << "]" << std::endl;
            }
         }
         for(int k=0; k<n; k+=block_size){ // same as i
            //copy A
            for(int ia = i; ia < i + block_size; ia++){
               for(int ka = k; ka < k + block_size; ka++){
                  memcpy(&Alocal[ia + ka * n], &A[ia + ka * n], sizeof(double)*block_size*block_size);
                  std::cout << "local: " << Alocal[ia + ka * n] << " " 
                        << "A: " << A[ia + ka* n] 
                        << " i[" << ka << "], j[" << ia << "]" << std::endl;
               }
            }
            //copy B
            for(int kb = k; kb < k + block_size; kb++){
               for(int jb = j; jb < j + block_size; jb++){
                  memcpy(&Blocal[kb + jb * n], &B[kb + jb* n], sizeof(double)*block_size*block_size);
                  std::cout << "local: " << Blocal[kb+jb* n] << " " 
                        << "B: " << B[kb + jb* n] 
                        << " i[" << jb << "], j[" << kb << "]"<< std::endl;
               }
            }
            
            // init ii=i, ii from i to i+block_size ->block_size times of iteration
            for (int ii=i; ii<i+block_size; ii++){
               for (int jj=j; jj<j+block_size; jj++){ //same as ii
                  for(int kk=k; kk<k+block_size; kk++){ //same as kk
                      // reference from basic : C[i + j * n] += A[i + k * n] * B[k + j * n];
                     Clocal[ii + jj * n] += Blocal[ii + kk * n] * Alocal[kk + jj * n];
                     // std::cout << "C[" << ii << "][" << jj <<"]: " << Clocal[ii + jj * block_size] 
                     //    << " A[" << ii << "][" << kk <<"]: " << Alocal[ii + kk * block_size] 
                     //    << " B[" << kk << "][" << jj <<"]: " << Clocal[ii + jj * block_size] << std::endl;
                  }
               }
            }
         }
         //write back to C
         for(int iic = i; iic < i + block_size; iic++){
            for(int jjc = j; jjc < j + block_size; jjc++){
               memcpy(&C[iic + jjc * n], &Clocal[iic + jjc * n], sizeof(double)*block_size*block_size);
               std::cout << "C: " << C[iic + jjc * n] << " " 
                        << "local: " << Clocal[iic + jjc * n] << std::endl;
            }
         }
      }
   }
   // std::cout << *A << " " << *B << " " << *C << std::endl;
}
