#include "complex_number.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

ComplexNumber complex_make(double real, double imag) {
    ComplexNumber z;
    z.real = real;
    z.imag = imag;
    return z;
}

ComplexNumber complex_zero(void) {
    return complex_make(0.0, 0.0);
}

ComplexNumber complex_one(void) {
    return complex_make(1.0, 0.0);
}

ComplexNumber complex_i(void) {
    return complex_make(0.0, 1.0);
}

ComplexNumber complex_from_polar(double radius, double angle) {
    return complex_make(radius * cos(angle), radius * sin(angle));
}

ComplexNumber complex_add(ComplexNumber a, ComplexNumber b) {
    return complex_make(a.real + b.real, a.imag + b.imag);
}

ComplexNumber complex_sub(ComplexNumber a, ComplexNumber b) {
    return complex_make(a.real - b.real, a.imag - b.imag);
}

ComplexNumber complex_mul(ComplexNumber a, ComplexNumber b) {
    return complex_make(a.real * b.real - a.imag * b.imag,
                        a.real * b.imag + a.imag * b.real);
}

ComplexNumber complex_div(ComplexNumber a, ComplexNumber b) {
    double denom = b.real * b.real + b.imag * b.imag;
    return complex_make((a.real * b.real + a.imag * b.imag) / denom,
                        (a.imag * b.real - a.real * b.imag) / denom);
}

ComplexNumber complex_scale(ComplexNumber a, double scale) {
    return complex_make(a.real * scale, a.imag * scale);
}

ComplexNumber complex_neg(ComplexNumber a) {
    return complex_make(-a.real, -a.imag);
}

ComplexNumber complex_conj(ComplexNumber a) {
    return complex_make(a.real, -a.imag);
}

double complex_abs(ComplexNumber a) {
    return hypot(a.real, a.imag);
}

double complex_abs2(ComplexNumber a) {
    return a.real * a.real + a.imag * a.imag;
}

double complex_distance(ComplexNumber a, ComplexNumber b) {
    return complex_abs(complex_sub(a, b));
}

int complex_almost_equal(ComplexNumber a, ComplexNumber b, double eps) {
    return complex_distance(a, b) <= eps;
}

void complex_array_zero(ComplexNumber *items, size_t count) {
    if (items == NULL || count == 0) {
        return;
    }
    memset(items, 0, count * sizeof(ComplexNumber));
}

void complex_array_copy(ComplexNumber *dst, const ComplexNumber *src, size_t count) {
    if (dst == NULL || src == NULL || count == 0) {
        return;
    }
    memcpy(dst, src, count * sizeof(ComplexNumber));
}

void complex_array_scale(ComplexNumber *items, size_t count, double scale) {
    if (items == NULL) {
        return;
    }
    for (size_t i = 0; i < count; ++i) {
        items[i] = complex_scale(items[i], scale);
    }
}

void complex_array_from_double(ComplexNumber *dst, const double *src, size_t count) {
    if (dst == NULL || src == NULL) {
        return;
    }
    for (size_t i = 0; i < count; ++i) {
        dst[i] = complex_make(src[i], 0.0);
    }
}

void complex_array_to_double(double *dst, const ComplexNumber *src, size_t count) {
    if (dst == NULL || src == NULL) {
        return;
    }
    for (size_t i = 0; i < count; ++i) {
        dst[i] = src[i].real;
    }
}

ComplexNumber *complex_array_alloc(size_t count) {
    if (count == 0) {
        return NULL;
    }
    return (ComplexNumber *)calloc(count, sizeof(ComplexNumber));
}

ComplexNumber *complex_array_clone(const ComplexNumber *src, size_t count) {
    if (src == NULL || count == 0) {
        return NULL;
    }
    ComplexNumber *dst = complex_array_alloc(count);
    if (dst == NULL) {
        return NULL;
    }
    complex_array_copy(dst, src, count);
    return dst;
}

void complex_array_free(ComplexNumber *items) {
    free(items);
}

double complex_array_max_imag_abs(const ComplexNumber *items, size_t count) {
    double max_value = 0.0;
    if (items == NULL) {
        return max_value;
    }
    for (size_t i = 0; i < count; ++i) {
        double value = fabs(items[i].imag);
        if (value > max_value) {
            max_value = value;
        }
    }
    return max_value;
}

double complex_array_max_error(const ComplexNumber *a,
                               const ComplexNumber *b,
                               size_t count) {
    double max_error = 0.0;
    if (a == NULL || b == NULL) {
        return INFINITY;
    }
    for (size_t i = 0; i < count; ++i) {
        double error = complex_distance(a[i], b[i]);
        if (error > max_error) {
            max_error = error;
        }
    }
    return max_error;
}
