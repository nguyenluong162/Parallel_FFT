#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <mpi.h>
#include <stddef.h>
#include "complex_number.h"
#include "fft.h"

typedef struct Polynomial {
    size_t degree;
    size_t length;
    double *coefficients;
} Polynomial;

typedef struct PolynomialProduct {
    Polynomial result;
    size_t fft_size;
    double elapsed_seconds;
    double compute_seconds;
    double communication_seconds;
    double max_rounding_error;
} PolynomialProduct;

int polynomial_init(Polynomial *poly, size_t degree);
void polynomial_destroy(Polynomial *poly);
int polynomial_resize(Polynomial *poly, size_t degree);
Polynomial polynomial_empty(void);

int polynomial_copy(Polynomial *dst, const Polynomial *src);
int polynomial_from_array(Polynomial *dst, const double *coeffs, size_t length);
void polynomial_fill_demo(Polynomial *poly, double seed);
void polynomial_print(const Polynomial *poly, const char *name);
void polynomial_print_expression(const Polynomial *poly, const char *name);
void polynomial_trim_small(Polynomial *poly, double eps);

int polynomial_multiply_serial(const Polynomial *a,
                               const Polynomial *b,
                               PolynomialProduct *product);
int polynomial_multiply_naive(const Polynomial *a,
                              const Polynomial *b,
                              Polynomial *result);
int polynomial_multiply_mpi(const Polynomial *a,
                            const Polynomial *b,
                            PolynomialProduct *product,
                            MPI_Comm comm);

double polynomial_max_abs_error(const Polynomial *a, const Polynomial *b);
int polynomial_almost_equal(const Polynomial *a, const Polynomial *b, double eps);
void polynomial_product_destroy(PolynomialProduct *product);
size_t polynomial_required_fft_size(const Polynomial *a, const Polynomial *b);

#endif
