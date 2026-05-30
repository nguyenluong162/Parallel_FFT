#include "polynomial.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi_fft.h"

static double now_mpi_seconds(void) {
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (initialized) {
        return MPI_Wtime();
    }
    return 0.0;
}

Polynomial polynomial_empty(void) {
    Polynomial poly;
    poly.degree = 0;
    poly.length = 0;
    poly.coefficients = NULL;
    return poly;
}

int polynomial_init(Polynomial *poly, size_t degree) {
    if (poly == NULL) {
        return 0;
    }
    poly->degree = degree;
    poly->length = degree + 1;
    poly->coefficients = (double *)calloc(poly->length, sizeof(double));
    if (poly->coefficients == NULL) {
        poly->degree = 0;
        poly->length = 0;
        return 0;
    }
    return 1;
}

void polynomial_destroy(Polynomial *poly) {
    if (poly == NULL) {
        return;
    }
    free(poly->coefficients);
    poly->coefficients = NULL;
    poly->degree = 0;
    poly->length = 0;
}

int polynomial_resize(Polynomial *poly, size_t degree) {
    if (poly == NULL) {
        return 0;
    }
    size_t new_length = degree + 1;
    double *new_coefficients = (double *)calloc(new_length, sizeof(double));
    if (new_coefficients == NULL) {
        return 0;
    }
    size_t copy_count = poly->length < new_length ? poly->length : new_length;
    for (size_t i = 0; i < copy_count; ++i) {
        new_coefficients[i] = poly->coefficients[i];
    }
    free(poly->coefficients);
    poly->coefficients = new_coefficients;
    poly->degree = degree;
    poly->length = new_length;
    return 1;
}

int polynomial_copy(Polynomial *dst, const Polynomial *src) {
    if (dst == NULL || src == NULL || src->coefficients == NULL) {
        return 0;
    }
    if (!polynomial_init(dst, src->degree)) {
        return 0;
    }
    memcpy(dst->coefficients, src->coefficients, src->length * sizeof(double));
    return 1;
}

int polynomial_from_array(Polynomial *dst, const double *coeffs, size_t length) {
    if (dst == NULL || coeffs == NULL || length == 0) {
        return 0;
    }
    if (!polynomial_init(dst, length - 1)) {
        return 0;
    }
    memcpy(dst->coefficients, coeffs, length * sizeof(double));
    return 1;
}

void polynomial_fill_demo(Polynomial *poly, double seed) {
    if (poly == NULL || poly->coefficients == NULL) {
        return;
    }
    for (size_t i = 0; i < poly->length; ++i) {
        double x = (double)(i + 1);
        double sign = (i % 2 == 0) ? 1.0 : -1.0;
        poly->coefficients[i] = sign * (seed + fmod(x * x + 3.0 * x, 17.0));
    }
}

void polynomial_print(const Polynomial *poly, const char *name) {
    if (poly == NULL || poly->coefficients == NULL) {
        printf("%s = <null>\n", name != NULL ? name : "poly");
        return;
    }
    printf("%s degree=%zu length=%zu\n", name != NULL ? name : "poly", poly->degree, poly->length);
    for (size_t i = 0; i < poly->length; ++i) {
        printf("  [%zu] %.1f\n", i, poly->coefficients[i]);
    }
}

void polynomial_print_expression(const Polynomial *poly, const char *name) {
    const char *label = name != NULL ? name : "P";
    if (poly == NULL || poly->coefficients == NULL) {
        printf("%s(x) = <null>\n", label);
        return;
    }

    printf("%s(x) = ", label);
    int printed = 0;
    for (size_t offset = 0; offset < poly->length; ++offset) {
        size_t power = poly->degree - offset;
        double coeff = poly->coefficients[power];
        if (fabs(coeff) < 1.0e-9) {
            continue;
        }

        long long abs_coeff = llround(fabs(coeff));
        if (!printed) {
            if (coeff < 0.0) {
                printf("-");
            }
        } else {
            printf(" %c ", coeff < 0.0 ? '-' : '+');
        }

        if (power == 0) {
            printf("%lld", abs_coeff);
        } else if (power == 1) {
            printf("%lldx", abs_coeff);
        } else {
            printf("%lldx^%zu", abs_coeff, power);
        }
        printed = 1;
    }

    if (!printed) {
        printf("0");
    }
    printf("\n");
}

void polynomial_trim_small(Polynomial *poly, double eps) {
    if (poly == NULL || poly->coefficients == NULL) {
        return;
    }
    for (size_t i = 0; i < poly->length; ++i) {
        if (fabs(poly->coefficients[i]) < eps) {
            poly->coefficients[i] = 0.0;
        }
    }
    while (poly->length > 1 && fabs(poly->coefficients[poly->length - 1]) < eps) {
        --poly->length;
    }
    poly->degree = poly->length - 1;
}

size_t polynomial_required_fft_size(const Polynomial *a, const Polynomial *b) {
    if (a == NULL || b == NULL) {
        return 0;
    }
    size_t result_length = a->length + b->length - 1;
    return next_power_of_two_size(result_length);
}

static void load_coefficients(ComplexNumber *dst,
                              size_t n,
                              const Polynomial *poly) {
    complex_array_zero(dst, n);
    for (size_t i = 0; i < poly->length; ++i) {
        dst[i] = complex_make(poly->coefficients[i], 0.0);
    }
}

static void extract_coefficients(Polynomial *poly,
                                 const ComplexNumber *values,
                                 size_t result_length,
                                 double *max_rounding_error) {
    double max_error = 0.0;
    for (size_t i = 0; i < result_length; ++i) {
        double rounded = nearbyint(values[i].real);
        double error = fabs(values[i].real - rounded);
        if (error > max_error) {
            max_error = error;
        }
        poly->coefficients[i] = fabs(error) < 1.0e-6 ? rounded : values[i].real;
    }
    if (max_rounding_error != NULL) {
        *max_rounding_error = max_error;
    }
}

int polynomial_multiply_naive(const Polynomial *a,
                              const Polynomial *b,
                              Polynomial *result) {
    if (a == NULL || b == NULL || result == NULL ||
        a->coefficients == NULL || b->coefficients == NULL) {
        return 0;
    }
    size_t result_degree = a->degree + b->degree;
    if (!polynomial_init(result, result_degree)) {
        return 0;
    }
    for (size_t i = 0; i < a->length; ++i) {
        for (size_t j = 0; j < b->length; ++j) {
            result->coefficients[i + j] += a->coefficients[i] * b->coefficients[j];
        }
    }
    return 1;
}

int polynomial_multiply_serial(const Polynomial *a,
                               const Polynomial *b,
                               PolynomialProduct *product) {
    if (a == NULL || b == NULL || product == NULL ||
        a->coefficients == NULL || b->coefficients == NULL) {
        return 0;
    }
    product->result = polynomial_empty();
    product->fft_size = polynomial_required_fft_size(a, b);
    product->elapsed_seconds = 0.0;
    product->max_rounding_error = 0.0;

    size_t n = product->fft_size;
    size_t result_length = a->length + b->length - 1;
    if (!polynomial_init(&product->result, result_length - 1)) {
        return 0;
    }

    ComplexNumber *fa = complex_array_alloc(n);
    ComplexNumber *fb = complex_array_alloc(n);
    if (fa == NULL || fb == NULL) {
        complex_array_free(fa);
        complex_array_free(fb);
        polynomial_product_destroy(product);
        return 0;
    }

    double start = now_mpi_seconds();
    load_coefficients(fa, n, a);
    load_coefficients(fb, n, b);
    fft_serial(fa, n, FFT_FORWARD);
    fft_serial(fb, n, FFT_FORWARD);
    for (size_t i = 0; i < n; ++i) {
        fa[i] = complex_mul(fa[i], fb[i]);
    }
    fft_serial(fa, n, FFT_INVERSE);
    product->elapsed_seconds = now_mpi_seconds() - start;
    extract_coefficients(&product->result,
                         fa,
                         result_length,
                         &product->max_rounding_error);

    complex_array_free(fa);
    complex_array_free(fb);
    return 1;
}

static void multiply_pointwise_replicated(ComplexNumber *fa,
                                          const ComplexNumber *fb,
                                          size_t n) {
    for (size_t i = 0; i < n; ++i) {
        fa[i] = complex_mul(fa[i], fb[i]);
    }
}

int polynomial_multiply_mpi(const Polynomial *a,
                            const Polynomial *b,
                            PolynomialProduct *product,
                            MPI_Comm comm) {
    if (a == NULL || b == NULL || product == NULL ||
        a->coefficients == NULL || b->coefficients == NULL) {
        return 0;
    }

    product->result = polynomial_empty();
    product->fft_size = polynomial_required_fft_size(a, b);
    product->elapsed_seconds = 0.0;
    product->max_rounding_error = 0.0;
    size_t n = product->fft_size;
    size_t result_length = a->length + b->length - 1;

    if (!polynomial_init(&product->result, result_length - 1)) {
        return 0;
    }

    ComplexNumber *fa = complex_array_alloc(n);
    ComplexNumber *fb = complex_array_alloc(n);
    if (fa == NULL || fb == NULL) {
        complex_array_free(fa);
        complex_array_free(fb);
        polynomial_product_destroy(product);
        return 0;
    }

    double start = MPI_Wtime();
    load_coefficients(fa, n, a);
    load_coefficients(fb, n, b);
    mpi_fft(fa, n, FFT_FORWARD, comm);
    mpi_fft(fb, n, FFT_FORWARD, comm);
    multiply_pointwise_replicated(fa, fb, n);
    mpi_fft(fa, n, FFT_INVERSE, comm);
    product->elapsed_seconds = MPI_Wtime() - start;

    extract_coefficients(&product->result,
                         fa,
                         result_length,
                         &product->max_rounding_error);

    complex_array_free(fa);
    complex_array_free(fb);
    return 1;
}

double polynomial_max_abs_error(const Polynomial *a, const Polynomial *b) {
    if (a == NULL || b == NULL || a->coefficients == NULL || b->coefficients == NULL) {
        return INFINITY;
    }
    size_t max_length = a->length > b->length ? a->length : b->length;
    double max_error = 0.0;
    for (size_t i = 0; i < max_length; ++i) {
        double av = i < a->length ? a->coefficients[i] : 0.0;
        double bv = i < b->length ? b->coefficients[i] : 0.0;
        double error = fabs(av - bv);
        if (error > max_error) {
            max_error = error;
        }
    }
    return max_error;
}

int polynomial_almost_equal(const Polynomial *a, const Polynomial *b, double eps) {
    return polynomial_max_abs_error(a, b) <= eps;
}

void polynomial_product_destroy(PolynomialProduct *product) {
    if (product == NULL) {
        return;
    }
    polynomial_destroy(&product->result);
    product->fft_size = 0;
    product->elapsed_seconds = 0.0;
    product->max_rounding_error = 0.0;
}
