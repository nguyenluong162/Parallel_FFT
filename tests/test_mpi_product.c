#include <mpi.h>
#include <math.h>
#include <stdio.h>

#include "polynomial.h"

static int build_case(Polynomial *a, Polynomial *b) {
    if (!polynomial_init(a, 15)) {
        return 0;
    }
    if (!polynomial_init(b, 15)) {
        polynomial_destroy(a);
        return 0;
    }
    polynomial_fill_demo(a, 7.0);
    polynomial_fill_demo(b, 11.0);
    return 1;
}

static int check_mpi_dnc(const Polynomial *a,
                         const Polynomial *b,
                         const Polynomial *naive,
                         int rank) {
    PolynomialProduct product;
    product.result = polynomial_empty();
    product.fft_size = 0;
    product.elapsed_seconds = 0.0;
    product.compute_seconds = 0.0;
    product.communication_seconds = 0.0;
    product.max_rounding_error = 0.0;

    int ok = polynomial_multiply_mpi(a, b, &product, MPI_COMM_WORLD);
    if (ok) {
        double timing_error = fabs(product.elapsed_seconds -
                                   product.compute_seconds -
                                   product.communication_seconds);
        ok = product.elapsed_seconds >= 0.0 &&
             product.compute_seconds >= 0.0 &&
             product.communication_seconds >= 0.0 &&
             timing_error <= 1.0e-9;
        if (!ok) {
            fprintf(stderr, "Rank %d timing accounting error %.6e\n",
                    rank, timing_error);
        }
    }
    if (rank == 0 && ok) {
        double error = polynomial_max_abs_error(naive, &product.result);
        ok = error <= 1.0e-6;
        if (!ok) {
            fprintf(stderr, "MPI divide-and-conquer product error %.6e\n", error);
        }
    }

    int all_ok = 0;
    MPI_Allreduce(&ok, &all_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    polynomial_product_destroy(&product);
    return all_ok;
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Polynomial a = polynomial_empty();
    Polynomial b = polynomial_empty();
    Polynomial naive = polynomial_empty();
    int ok = build_case(&a, &b);
    if (rank == 0) {
        ok = ok && polynomial_multiply_naive(&a, &b, &naive);
    }

    int all_ok = 0;
    MPI_Allreduce(&ok, &all_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    ok = all_ok;
    ok = ok && check_mpi_dnc(&a, &b, &naive, rank);

    all_ok = 0;
    MPI_Allreduce(&ok, &all_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    polynomial_destroy(&a);
    polynomial_destroy(&b);
    polynomial_destroy(&naive);
    MPI_Finalize();

    if (rank == 0) {
        if (all_ok) {
            printf("test_mpi_product: PASS\n");
        } else {
            fprintf(stderr, "test_mpi_product: FAIL\n");
        }
    }
    return all_ok ? 0 : 1;
}
