#include "fft.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef FFT_PI
#define FFT_PI 3.141592653589793238462643383279502884
#endif

static double monotonic_seconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1.0e-9;
}

int is_power_of_two_size(size_t n) {
    return n != 0 && (n & (n - 1)) == 0;
}

size_t next_power_of_two_size(size_t n) {
    if (n <= 1) {
        return 1;
    }
    size_t value = 1;
    while (value < n) {
        value <<= 1;
    }
    return value;
}

unsigned int log2_size_exact(size_t n) {
    unsigned int bits = 0;
    while (n > 1) {
        n >>= 1;
        ++bits;
    }
    return bits;
}

size_t bit_reverse_size(size_t value, unsigned int bits) {
    size_t reversed = 0;
    for (unsigned int i = 0; i < bits; ++i) {
        reversed <<= 1;
        reversed |= (value & 1U);
        value >>= 1;
    }
    return reversed;
}

void fft_bit_reverse_permute(ComplexNumber *values, size_t n) {
    if (values == NULL || n <= 1) {
        return;
    }
    unsigned int bits = log2_size_exact(n);
    for (size_t i = 0; i < n; ++i) {
        size_t j = bit_reverse_size(i, bits);
        if (j > i) {
            ComplexNumber tmp = values[i];
            values[i] = values[j];
            values[j] = tmp;
        }
    }
}

void fft_normalize_inverse(ComplexNumber *values, size_t n) {
    if (values == NULL || n == 0) {
        return;
    }
    double scale = 1.0 / (double)n;
    for (size_t i = 0; i < n; ++i) {
        values[i] = complex_scale(values[i], scale);
    }
}

void fft_serial(ComplexNumber *values, size_t n, FftDirection direction) {
    fft_serial_with_stats(values, n, direction, NULL);
}

void fft_serial_with_stats(ComplexNumber *values,
                           size_t n,
                           FftDirection direction,
                           FftStats *stats) {
    double start = monotonic_seconds();
    if (stats != NULL) {
        stats->n = n;
        stats->stages = 0;
        stats->butterflies = 0;
        stats->elapsed_seconds = 0.0;
    }
    if (values == NULL || !is_power_of_two_size(n)) {
        return;
    }

    fft_bit_reverse_permute(values, n);

    for (size_t len = 2; len <= n; len <<= 1) {
        double angle = (direction == FFT_FORWARD ? -2.0 : 2.0) * FFT_PI / (double)len;
        ComplexNumber w_len = complex_from_polar(1.0, angle);
        for (size_t block = 0; block < n; block += len) {
            ComplexNumber w = complex_one();
            size_t half = len >> 1;
            for (size_t j = 0; j < half; ++j) {
                ComplexNumber u = values[block + j];
                ComplexNumber v = complex_mul(values[block + j + half], w);
                values[block + j] = complex_add(u, v);
                values[block + j + half] = complex_sub(u, v);
                w = complex_mul(w, w_len);
                if (stats != NULL) {
                    ++stats->butterflies;
                }
            }
        }
        if (stats != NULL) {
            ++stats->stages;
        }
    }

    if (direction == FFT_INVERSE) {
        fft_normalize_inverse(values, n);
    }
    if (stats != NULL) {
        stats->elapsed_seconds = monotonic_seconds() - start;
    }
}

void dft_slow(const ComplexNumber *input,
              ComplexNumber *output,
              size_t n,
              FftDirection direction) {
    if (input == NULL || output == NULL || n == 0) {
        return;
    }
    double sign = direction == FFT_FORWARD ? -1.0 : 1.0;
    for (size_t k = 0; k < n; ++k) {
        ComplexNumber sum = complex_zero();
        for (size_t t = 0; t < n; ++t) {
            double angle = sign * 2.0 * FFT_PI * (double)(k * t) / (double)n;
            ComplexNumber w = complex_from_polar(1.0, angle);
            sum = complex_add(sum, complex_mul(input[t], w));
        }
        if (direction == FFT_INVERSE) {
            sum = complex_scale(sum, 1.0 / (double)n);
        }
        output[k] = sum;
    }
}

double fft_compare_with_dft(const ComplexNumber *input, size_t n) {
    if (input == NULL || n == 0 || !is_power_of_two_size(n)) {
        return INFINITY;
    }
    ComplexNumber *fft_values = complex_array_clone(input, n);
    ComplexNumber *dft_values = complex_array_alloc(n);
    if (fft_values == NULL || dft_values == NULL) {
        complex_array_free(fft_values);
        complex_array_free(dft_values);
        return INFINITY;
    }
    fft_serial(fft_values, n, FFT_FORWARD);
    dft_slow(input, dft_values, n, FFT_FORWARD);
    double error = complex_array_max_error(fft_values, dft_values, n);
    complex_array_free(fft_values);
    complex_array_free(dft_values);
    return error;
}
