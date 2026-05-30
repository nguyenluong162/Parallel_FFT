#ifndef MPI_FFT_H
#define MPI_FFT_H

#include <mpi.h>
#include <stddef.h>
#include "complex_number.h"
#include "fft.h"

void mpi_fft(ComplexNumber *values,
             size_t n,
             FftDirection direction,
             MPI_Comm comm);
void mpi_fft_dnc(ComplexNumber *values,
                 size_t n,
                 FftDirection direction,
                 MPI_Comm comm);

void mpi_complex_bcast(ComplexNumber *values,
                       size_t count,
                       int root,
                       MPI_Comm comm);

#endif
