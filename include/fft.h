#ifndef FFT_H
#define FFT_H

#include <stddef.h>
#include "complex_number.h"

typedef enum FftDirection {
    FFT_FORWARD = 1,
    FFT_INVERSE = -1
} FftDirection;

typedef struct FftStats {
    size_t n;
    size_t stages;
    size_t butterflies;
    double elapsed_seconds;
} FftStats;

int is_power_of_two_size(size_t n);
size_t next_power_of_two_size(size_t n);
size_t bit_reverse_size(size_t value, unsigned int bits);
unsigned int log2_size_exact(size_t n);

void fft_bit_reverse_permute(ComplexNumber *values, size_t n);
void fft_serial(ComplexNumber *values, size_t n, FftDirection direction);
void fft_serial_with_stats(ComplexNumber *values,
                           size_t n,
                           FftDirection direction,
                           FftStats *stats);
void dft_slow(const ComplexNumber *input,
              ComplexNumber *output,
              size_t n,
              FftDirection direction);

double fft_compare_with_dft(const ComplexNumber *input, size_t n);
void fft_normalize_inverse(ComplexNumber *values, size_t n);

#endif
