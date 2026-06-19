#include "mpi_fft.h"

#include <math.h>
#include <stdlib.h>

#ifndef FFT_PI
#define FFT_PI 3.141592653589793238462643383279502884
#endif

static double profile_communication_seconds = 0.0;

void mpi_fft_profile_reset(void) {
    profile_communication_seconds = 0.0;
}

double mpi_fft_profile_communication_seconds(void) {
    return profile_communication_seconds;
}

void mpi_fft_profile_barrier(MPI_Comm comm) {
    double start = MPI_Wtime();
    MPI_Barrier(comm);
    profile_communication_seconds += MPI_Wtime() - start;
}

void mpi_complex_bcast(ComplexNumber *values,
                       size_t count,
                       int root,
                       MPI_Comm comm) {
    if (values == NULL || count == 0) {
        return;
    }
    double start = MPI_Wtime();
    MPI_Bcast((double *)values, (int)(2 * count), MPI_DOUBLE, root, comm);
    profile_communication_seconds += MPI_Wtime() - start;
}

static void fft_serial_unscaled(ComplexNumber *values,
                                size_t n,
                                FftDirection direction) {
    if (values == NULL || !is_power_of_two_size(n)) {
        return;
    }

    fft_bit_reverse_permute(values, n);
    for (size_t len = 2; len <= n; len <<= 1) {
        size_t half = len >> 1;
        double angle = (direction == FFT_FORWARD ? -2.0 : 2.0) * FFT_PI / (double)len;
        ComplexNumber w_len = complex_from_polar(1.0, angle);
        for (size_t block = 0; block < n; block += len) {
            ComplexNumber w = complex_one();
            for (size_t j = 0; j < half; ++j) {
                ComplexNumber u = values[block + j];
                ComplexNumber v = complex_mul(values[block + j + half], w);
                values[block + j] = complex_add(u, v);
                values[block + j + half] = complex_sub(u, v);
                w = complex_mul(w, w_len);
            }
        }
    }
}

static void split_even_odd_for_branch(const ComplexNumber *values,
                                      ComplexNumber *branch_values,
                                      size_t n,
                                      int use_odd_branch) {
    size_t half = n >> 1;
    size_t start = use_odd_branch ? 1U : 0U;
    for (size_t i = 0; i < half; ++i) {
        branch_values[i] = values[start + 2U * i];
    }
}

static void combine_dnc_results(ComplexNumber *values,
                                const ComplexNumber *even_values,
                                const ComplexNumber *odd_values,
                                size_t n,
                                FftDirection direction) {
    size_t half = n >> 1;
    double angle_unit = (direction == FFT_FORWARD ? -2.0 : 2.0) * FFT_PI / (double)n;
    for (size_t k = 0; k < half; ++k) {
        ComplexNumber w = complex_from_polar(1.0, angle_unit * (double)k);
        ComplexNumber odd_rotated = complex_mul(w, odd_values[k]);
        values[k] = complex_add(even_values[k], odd_rotated);
        values[k + half] = complex_sub(even_values[k], odd_rotated);
    }
}

static void mpi_fft_dnc_kernel(ComplexNumber *values,
                               size_t n,
                               FftDirection direction,
                               MPI_Comm comm) {
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if (values == NULL || n <= 1) {
        return;
    }

    if (size <= 1) {
        fft_serial_unscaled(values, n, direction);
        return;
    }

    int left_size = size / 2;
    int color = rank < left_size ? 0 : 1;
    int parent_root_even = 0;
    int parent_root_odd = left_size;
    MPI_Comm branch_comm = MPI_COMM_NULL;
    double communication_start = MPI_Wtime();
    MPI_Comm_split(comm, color, rank, &branch_comm);
    profile_communication_seconds += MPI_Wtime() - communication_start;

    size_t half = n >> 1;
    ComplexNumber *branch_values = complex_array_alloc(half);
    ComplexNumber *even_values = complex_array_alloc(half);
    ComplexNumber *odd_values = complex_array_alloc(half);
    int local_ok = branch_values != NULL && even_values != NULL && odd_values != NULL;
    int all_ok = 0;
    communication_start = MPI_Wtime();
    MPI_Allreduce(&local_ok, &all_ok, 1, MPI_INT, MPI_MIN, comm);
    profile_communication_seconds += MPI_Wtime() - communication_start;
    if (!all_ok) {
        complex_array_free(branch_values);
        complex_array_free(even_values);
        complex_array_free(odd_values);
        communication_start = MPI_Wtime();
        MPI_Comm_free(&branch_comm);
        profile_communication_seconds += MPI_Wtime() - communication_start;
        return;
    }

    split_even_odd_for_branch(values, branch_values, n, color == 1);
    mpi_fft_dnc_kernel(branch_values, half, direction, branch_comm);

    if (rank == parent_root_even) {
        complex_array_copy(even_values, branch_values, half);
    }
    if (rank == parent_root_odd) {
        complex_array_copy(odd_values, branch_values, half);
    }

    mpi_complex_bcast(even_values, half, parent_root_even, comm);
    mpi_complex_bcast(odd_values, half, parent_root_odd, comm);
    combine_dnc_results(values, even_values, odd_values, n, direction);

    complex_array_free(branch_values);
    complex_array_free(even_values);
    complex_array_free(odd_values);
    communication_start = MPI_Wtime();
    MPI_Comm_free(&branch_comm);
    profile_communication_seconds += MPI_Wtime() - communication_start;
}

void mpi_fft_dnc(ComplexNumber *values,
                 size_t n,
                 FftDirection direction,
                 MPI_Comm comm) {
    if (values == NULL || !is_power_of_two_size(n)) {
        return;
    }

    mpi_fft_dnc_kernel(values, n, direction, comm);
    if (direction == FFT_INVERSE) {
        fft_normalize_inverse(values, n);
    }
}

void mpi_fft(ComplexNumber *values,
             size_t n,
             FftDirection direction,
             MPI_Comm comm) {
    mpi_fft_dnc(values, n, direction, comm);
}
