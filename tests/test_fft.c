#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "complex_number.h"
#include "fft.h"

static int check_close(const char *label, double actual, double expected, double eps) {
    double error = fabs(actual - expected);
    if (error > eps) {
        fprintf(stderr,
                "%s failed: actual=%.12f expected=%.12f error=%.6e\n",
                label,
                actual,
                expected,
                error);
        return 0;
    }
    return 1;
}

static int test_power_of_two_helpers(void) {
    int ok = 1;
    ok = ok && is_power_of_two_size(1);
    ok = ok && is_power_of_two_size(2);
    ok = ok && is_power_of_two_size(1024);
    ok = ok && !is_power_of_two_size(0);
    ok = ok && !is_power_of_two_size(3);
    ok = ok && next_power_of_two_size(0) == 1;
    ok = ok && next_power_of_two_size(1) == 1;
    ok = ok && next_power_of_two_size(9) == 16;
    ok = ok && log2_size_exact(8) == 3;
    ok = ok && bit_reverse_size(3, 3) == 6;
    return ok;
}

static int test_fft_matches_dft(void) {
    const size_t n = 8;
    ComplexNumber input[8];
    for (size_t i = 0; i < n; ++i) {
        input[i] = complex_make((double)(i + 1), (double)(i % 3) - 1.0);
    }
    double error = fft_compare_with_dft(input, n);
    if (error > 1.0e-9) {
        fprintf(stderr, "FFT/DFT comparison failed: %.6e\n", error);
        return 0;
    }
    return 1;
}

static int test_inverse_round_trip(void) {
    const size_t n = 16;
    ComplexNumber values[16];
    ComplexNumber original[16];
    for (size_t i = 0; i < n; ++i) {
        values[i] = complex_make(sin((double)i), cos((double)(2 * i)));
        original[i] = values[i];
    }
    fft_serial(values, n, FFT_FORWARD);
    fft_serial(values, n, FFT_INVERSE);
    double error = complex_array_max_error(values, original, n);
    if (error > 1.0e-9) {
        fprintf(stderr, "FFT round trip failed: %.6e\n", error);
        return 0;
    }
    return 1;
}

static int test_known_transform(void) {
    ComplexNumber values[4];
    values[0] = complex_make(1.0, 0.0);
    values[1] = complex_make(2.0, 0.0);
    values[2] = complex_make(3.0, 0.0);
    values[3] = complex_make(4.0, 0.0);
    fft_serial(values, 4, FFT_FORWARD);
    int ok = 1;
    ok = ok && check_close("X0.real", values[0].real, 10.0, 1.0e-9);
    ok = ok && check_close("X0.imag", values[0].imag, 0.0, 1.0e-9);
    ok = ok && check_close("X1.real", values[1].real, -2.0, 1.0e-9);
    ok = ok && check_close("X1.imag", values[1].imag, 2.0, 1.0e-9);
    ok = ok && check_close("X2.real", values[2].real, -2.0, 1.0e-9);
    ok = ok && check_close("X2.imag", values[2].imag, 0.0, 1.0e-9);
    ok = ok && check_close("X3.real", values[3].real, -2.0, 1.0e-9);
    ok = ok && check_close("X3.imag", values[3].imag, -2.0, 1.0e-9);
    return ok;
}

int main(void) {
    int ok = 1;
    ok = ok && test_power_of_two_helpers();
    ok = ok && test_fft_matches_dft();
    ok = ok && test_inverse_round_trip();
    ok = ok && test_known_transform();
    if (ok) {
        printf("test_fft: PASS\n");
        return 0;
    }
    fprintf(stderr, "test_fft: FAIL\n");
    return 1;
}
