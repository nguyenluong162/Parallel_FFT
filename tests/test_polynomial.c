#include <stdio.h>

#include "polynomial.h"

static int check_polynomial_case(const double *a_values,
                                 size_t a_length,
                                 const double *b_values,
                                 size_t b_length) {
    Polynomial a = polynomial_empty();
    Polynomial b = polynomial_empty();
    Polynomial naive = polynomial_empty();
    PolynomialProduct product;
    product.result = polynomial_empty();
    product.fft_size = 0;
    product.elapsed_seconds = 0.0;
    product.compute_seconds = 0.0;
    product.communication_seconds = 0.0;
    product.max_rounding_error = 0.0;

    int ok = polynomial_from_array(&a, a_values, a_length) &&
             polynomial_from_array(&b, b_values, b_length) &&
             polynomial_multiply_naive(&a, &b, &naive) &&
             polynomial_multiply_serial(&a, &b, &product);

    if (ok) {
        double error = polynomial_max_abs_error(&naive, &product.result);
        ok = error <= 1.0e-6;
        if (!ok) {
            fprintf(stderr, "polynomial error %.6e\n", error);
        }
    }

    polynomial_destroy(&a);
    polynomial_destroy(&b);
    polynomial_destroy(&naive);
    polynomial_product_destroy(&product);
    return ok;
}

static int test_small_cases(void) {
    double a1[] = {1.0};
    double b1[] = {5.0};
    double a2[] = {1.0, 2.0, 3.0};
    double b2[] = {4.0, 5.0, 6.0};
    double a3[] = {-2.0, 0.0, 7.0, 1.0};
    double b3[] = {3.0, -1.0};
    int ok = 1;
    ok = ok && check_polynomial_case(a1, 1, b1, 1);
    ok = ok && check_polynomial_case(a2, 3, b2, 3);
    ok = ok && check_polynomial_case(a3, 4, b3, 2);
    return ok;
}

static int test_demo_case(void) {
    Polynomial a = polynomial_empty();
    Polynomial b = polynomial_empty();
    Polynomial naive = polynomial_empty();
    PolynomialProduct product;
    product.result = polynomial_empty();
    product.fft_size = 0;
    product.elapsed_seconds = 0.0;
    product.compute_seconds = 0.0;
    product.communication_seconds = 0.0;
    product.max_rounding_error = 0.0;

    int ok = polynomial_init(&a, 31) && polynomial_init(&b, 31);
    if (ok) {
        polynomial_fill_demo(&a, 3.0);
        polynomial_fill_demo(&b, 5.0);
        ok = polynomial_multiply_naive(&a, &b, &naive) &&
             polynomial_multiply_serial(&a, &b, &product);
    }
    if (ok) {
        double error = polynomial_max_abs_error(&naive, &product.result);
        ok = error <= 1.0e-6;
        if (!ok) {
            fprintf(stderr, "demo polynomial error %.6e\n", error);
        }
    }

    polynomial_destroy(&a);
    polynomial_destroy(&b);
    polynomial_destroy(&naive);
    polynomial_product_destroy(&product);
    return ok;
}

int main(void) {
    int ok = test_small_cases() && test_demo_case();
    if (ok) {
        printf("test_polynomial: PASS\n");
        return 0;
    }
    fprintf(stderr, "test_polynomial: FAIL\n");
    return 1;
}
