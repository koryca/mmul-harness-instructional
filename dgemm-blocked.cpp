const char* dgemm_desc = "Blocked dgemm.";

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are n-by-n matrices stored in column-major format.
 * On exit, A and B maintain their input values. */
void square_dgemm_blocked(int n, int block_size, double* A, double* B, double* C) 
{
   // the block size is b=n/Nb where Nb is the # of blocks in each row/col
   
   int nb = n/block_size;
   for (int j=0; j<nb; j++){
      for (int i=0; i<nb; i++){
         //copy to local
         std::vector<double> buf(6 * n * n);
         double * Clocal = buf.data() + nb * nb;
         memcpy((void *)Clocal, (const void *)C, sizeof(double)*nb*nb);
         for(int k=0; k<nb; k++){
            // C[i,j] += A[i,k] * B[k,j]
            // C[i*nb+j] += A[i*nb+k] * B[k*nb+j];
            double * Alocal = Clocal + nb * nb;
            double * Blocal = Alocal + nb * nb;
            memcpy((void *)Alocal, (const void *)A, sizeof(double)*nb*nb);
            memcpy((void *)Blocal, (const void *)B, sizeof(double)*nb*nb);
            for (int m=0; m<block_size; m++){
               for (int p=0; p<block_size; p++){
                  for(int l=0; l<block_size; l++){
                     // C[i,j] += A[i,k] * B[k,j]
                     Clocal[p*block_size+m] += Alocal[p*block_size+l] * B[l*block_size+m];
                  }
               }
            }
            memcpy((void *)C, (const void *)Clocal, sizeof(double)*nb*nb);
         }
      }
   }
}
