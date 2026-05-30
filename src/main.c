#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "io_utils.h"
#include "polynomial.h"

static int root_rank(void) {
    return 0;
}

static int should_show_help(int argc, char **argv) {
    for (int i = 1; i < argc; ++i) {
        if (argv[i] != NULL &&
            (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0)) {
            return 1;
        }
    }
    return 0;
}

static InputStatus load_or_generate_input(const ProgramOptions *options,
                                          Polynomial *a,
                                          Polynomial *b) {
    if (options == NULL || a == NULL || b == NULL) {
        return INPUT_ERR_ARGS;
    }
    if (options->input_path != NULL) {
        return read_polynomials_from_file(options->input_path, a, b);
    }
    if (!polynomial_init(a, (size_t)options->demo_degree)) {
        return INPUT_ERR_MEMORY;
    }
    if (!polynomial_init(b, (size_t)options->demo_degree)) {
        polynomial_destroy(a);
        return INPUT_ERR_MEMORY;
    }
    polynomial_fill_demo(a, 1.0);
    polynomial_fill_demo(b, 2.0);
    return INPUT_OK;
}

static void print_summary_header(const ProgramOptions *options,
                                 const Polynomial *a,
                                 const Polynomial *b,
                                 int comm_size) {
    if (options == NULL || options->quiet) {
        return;
    }
    printf("Cooley-Tukey FFT polynomial multiplication\n");
    printf("  processors: %d\n", comm_size);
    if (options->use_serial) {
        printf("  mode: serial FFT\n");
    } else {
        printf("  mode: MPI FFT (divide-and-conquer)\n");
    }
    printf("  degree(A): %zu\n", a->degree);
    printf("  degree(B): %zu\n", b->degree);
}

static int verify_result(const Polynomial *a,
                         const Polynomial *b,
                         const Polynomial *candidate,
                         int quiet) {
    Polynomial expected = polynomial_empty();
    if (!polynomial_multiply_naive(a, b, &expected)) {
        fprintf(stderr, "Verification failed: could not compute naive product.\n");
        return 0;
    }
    double error = polynomial_max_abs_error(&expected, candidate);
    int ok = error <= 1.0e-6;
    if (!quiet) {
        printf("  verification: %s (max abs error %.6e)\n", ok ? "PASS" : "FAIL", error);
    } else {
        printf("verification=%s error=%.6e\n", ok ? "PASS" : "FAIL", error);
    }
    polynomial_destroy(&expected);
    return ok;
}

static int run_root_serial(const Polynomial *a,
                           const Polynomial *b,
                           PolynomialProduct *product) {
    return polynomial_multiply_serial(a, b, product);
}

static int run_mpi_product(const Polynomial *a,
                           const Polynomial *b,
                           PolynomialProduct *product,
                           MPI_Comm comm) {
    return polynomial_multiply_mpi(a, b, product, comm);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank = 0;
    int comm_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    if (should_show_help(argc, argv)) {
        if (rank == root_rank()) {
            print_usage(argv[0]);
        }
        MPI_Finalize();
        return 0;
    }

    ProgramOptions options;
    InputStatus status = parse_program_options(argc, argv, &options, rank);
    if (status != INPUT_OK) {
        if (rank == root_rank()) {
            fprintf(stderr, "Error: %s\n", input_status_message(status));
            print_usage(argv[0]);
        }
        MPI_Finalize();
        return (int)status;
    }

    Polynomial a = polynomial_empty();
    Polynomial b = polynomial_empty();
    InputStatus load_status = INPUT_OK;

    if (rank == root_rank()) {
        load_status = load_or_generate_input(&options, &a, &b);
        if (load_status != INPUT_OK) {
            fprintf(stderr, "Error loading input: %s\n", input_status_message(load_status));
        }
    }

    int load_code = (int)load_status;
    MPI_Bcast(&load_code, 1, MPI_INT, root_rank(), MPI_COMM_WORLD);
    if (load_code != (int)INPUT_OK) {
        polynomial_destroy(&a);
        polynomial_destroy(&b);
        MPI_Finalize();
        return load_code;
    }

    status = broadcast_polynomial(&a, root_rank(), MPI_COMM_WORLD);
    if (status == INPUT_OK) {
        status = broadcast_polynomial(&b, root_rank(), MPI_COMM_WORLD);
    }
    int broadcast_code = (int)status;
    MPI_Allreduce(MPI_IN_PLACE, &broadcast_code, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (broadcast_code != (int)INPUT_OK) {
        if (rank == root_rank()) {
            fprintf(stderr, "Error broadcasting input: %s\n",
                    input_status_message((InputStatus)broadcast_code));
        }
        polynomial_destroy(&a);
        polynomial_destroy(&b);
        MPI_Finalize();
        return broadcast_code;
    }

    if (rank == root_rank()) {
        print_summary_header(&options, &a, &b, comm_size);
        if (options.print_vectors) {
            polynomial_print_expression(&a, "A");
            polynomial_print_expression(&b, "B");
        }
    }

    PolynomialProduct product;
    product.result = polynomial_empty();
    product.fft_size = 0;
    product.elapsed_seconds = 0.0;
    product.max_rounding_error = 0.0;

    int ok = 1;
    if (options.use_serial) {
        if (rank == root_rank()) {
            ok = run_root_serial(&a, &b, &product);
        }
        MPI_Bcast(&ok, 1, MPI_INT, root_rank(), MPI_COMM_WORLD);
    } else {
        ok = run_mpi_product(&a, &b, &product, MPI_COMM_WORLD);
        int all_ok = 0;
        MPI_Allreduce(&ok, &all_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        ok = all_ok;
    }

    if (!ok) {
        if (rank == root_rank()) {
            fprintf(stderr, "Error: multiplication failed.\n");
        }
        polynomial_product_destroy(&product);
        polynomial_destroy(&a);
        polynomial_destroy(&b);
        MPI_Finalize();
        return 10;
    }

    int verified = 1;
    if (rank == root_rank()) {
        polynomial_trim_small(&product.result, 1.0e-7);
        if (!options.quiet) {
            printf("  fft size: %zu\n", product.fft_size);
            printf("  elapsed: %.6f seconds\n", product.elapsed_seconds);
            printf("  max rounding error: %.6e\n", product.max_rounding_error);
        } else {
            printf("fft_size=%zu elapsed=%.6f rounding=%.6e\n",
                   product.fft_size,
                   product.elapsed_seconds,
                   product.max_rounding_error);
        }

        if (options.verify) {
            verified = verify_result(&a, &b, &product.result, options.quiet);
        }

        if (options.print_vectors) {
            polynomial_print_expression(&product.result, "C");
        }

        if (options.output_path != NULL) {
            InputStatus write_status = write_polynomial_to_file(options.output_path,
                                                                &product.result);
            if (write_status != INPUT_OK) {
                fprintf(stderr, "Error writing output: %s\n",
                        input_status_message(write_status));
                verified = 0;
            } else if (!options.quiet) {
                printf("  output: %s\n", options.output_path);
            }
        }
    }

    polynomial_product_destroy(&product);
    polynomial_destroy(&a);
    polynomial_destroy(&b);
    MPI_Finalize();
    return verified ? 0 : 20;
}
