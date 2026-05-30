#ifndef COMPLEX_NUMBER_H
#define COMPLEX_NUMBER_H

#include <stddef.h>

typedef struct ComplexNumber {
    double real;
    double imag;
} ComplexNumber;

ComplexNumber complex_make(double real, double imag);
ComplexNumber complex_zero(void);
ComplexNumber complex_one(void);
ComplexNumber complex_i(void);
ComplexNumber complex_from_polar(double radius, double angle);

ComplexNumber complex_add(ComplexNumber a, ComplexNumber b);
ComplexNumber complex_sub(ComplexNumber a, ComplexNumber b);
ComplexNumber complex_mul(ComplexNumber a, ComplexNumber b);
ComplexNumber complex_div(ComplexNumber a, ComplexNumber b);
ComplexNumber complex_scale(ComplexNumber a, double scale);
ComplexNumber complex_neg(ComplexNumber a);
ComplexNumber complex_conj(ComplexNumber a);

double complex_abs(ComplexNumber a);
double complex_abs2(ComplexNumber a);
double complex_distance(ComplexNumber a, ComplexNumber b);
int complex_almost_equal(ComplexNumber a, ComplexNumber b, double eps);

void complex_array_zero(ComplexNumber *items, size_t count);
void complex_array_copy(ComplexNumber *dst, const ComplexNumber *src, size_t count);
void complex_array_scale(ComplexNumber *items, size_t count, double scale);
void complex_array_from_double(ComplexNumber *dst, const double *src, size_t count);
void complex_array_to_double(double *dst, const ComplexNumber *src, size_t count);

ComplexNumber *complex_array_alloc(size_t count);
ComplexNumber *complex_array_clone(const ComplexNumber *src, size_t count);
void complex_array_free(ComplexNumber *items);

double complex_array_max_imag_abs(const ComplexNumber *items, size_t count);
double complex_array_max_error(const ComplexNumber *a,
                               const ComplexNumber *b,
                               size_t count);

#endif
