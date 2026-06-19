#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

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

static int append_timing_csv(const ProgramOptions *options,
                             const Polynomial *a,
                             const Polynomial *b,
                             const PolynomialProduct *product,
                             double program_start,
                             MPI_Comm comm) {
    if (options == NULL || options->timing_csv_path == NULL) {
        return 1;
    }

    int rank = 0;
    int comm_size = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);

    MPI_Barrier(comm);
    double local_values[4] = {
        product->elapsed_seconds,
        product->compute_seconds,
        product->communication_seconds,
        MPI_Wtime() - program_start
    };

    char processor_name[MPI_MAX_PROCESSOR_NAME] = {0};
    int processor_name_length = 0;
    MPI_Get_processor_name(processor_name, &processor_name_length);
    (void)processor_name_length;

    double *all_values = NULL;
    char *all_names = NULL;
    int export_ok = 1;
    if (rank == root_rank()) {
        all_values = (double *)calloc((size_t)comm_size * 4U, sizeof(double));
        all_names = (char *)calloc((size_t)comm_size,
                                   (size_t)MPI_MAX_PROCESSOR_NAME);
        export_ok = all_values != NULL && all_names != NULL;
    }
    MPI_Bcast(&export_ok, 1, MPI_INT, root_rank(), comm);
    if (!export_ok) {
        free(all_values);
        free(all_names);
        return 0;
    }

    MPI_Gather(local_values, 4, MPI_DOUBLE,
               all_values, 4, MPI_DOUBLE, root_rank(), comm);
    MPI_Gather(processor_name, MPI_MAX_PROCESSOR_NAME, MPI_CHAR,
               all_names, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, root_rank(), comm);

    if (rank == root_rank()) {
        FILE *file = fopen(options->timing_csv_path, "a+");
        if (file == NULL) {
            export_ok = 0;
        } else {
            if (fseek(file, 0L, SEEK_END) != 0) {
                export_ok = 0;
            } else if (ftell(file) == 0L) {
                fprintf(file,
                        "run_id,mode,input_size,degree_a,degree_b,fft_size,"
                        "processes,rank,host,kernel_total_s,compute_s,"
                        "communication_s,program_s\n");
            }

            char generated_run_id[64];
            const char *run_id = options->run_id;
            if (run_id == NULL) {
                snprintf(generated_run_id, sizeof(generated_run_id),
                         "run-%lld", (long long)time(NULL));
                run_id = generated_run_id;
            }
            size_t input_size = a->length > b->length ? a->length : b->length;
            const char *mode = options->use_serial ? "serial" : "mpi";
            for (int i = 0; export_ok && i < comm_size; ++i) {
                const double *values = &all_values[4 * i];
                const char *name = &all_names[(size_t)i * MPI_MAX_PROCESSOR_NAME];
                if (fprintf(file,
                            "%s,%s,%zu,%zu,%zu,%zu,%d,%d,%s,"
                            "%.9f,%.9f,%.9f,%.9f\n",
                            run_id, mode, input_size, a->degree, b->degree,
                            product->fft_size, comm_size, i, name,
                            values[0], values[1], values[2], values[3]) < 0) {
                    export_ok = 0;
                }
            }
            if (fclose(file) != 0) {
                export_ok = 0;
            }
        }
    }

    MPI_Bcast(&export_ok, 1, MPI_INT, root_rank(), comm);
    free(all_values);
    free(all_names);
    return export_ok;
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    double program_start = MPI_Wtime();

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
    product.compute_seconds = 0.0;
    product.communication_seconds = 0.0;
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
            printf("  compute: %.6f seconds\n", product.compute_seconds);
            printf("  communication: %.6f seconds\n", product.communication_seconds);
            printf("  max rounding error: %.6e\n", product.max_rounding_error);
        } else {
            printf("fft_size=%zu elapsed=%.6f compute=%.6f communication=%.6f "
                   "rounding=%.6e\n",
                   product.fft_size,
                   product.elapsed_seconds,
                   product.compute_seconds,
                   product.communication_seconds,
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

    if (!append_timing_csv(&options, &a, &b, &product,
                           program_start, MPI_COMM_WORLD)) {
        if (rank == root_rank()) {
            fprintf(stderr, "Error writing timing CSV: %s\n",
                    options.timing_csv_path);
        }
        verified = 0;
    }

    polynomial_product_destroy(&product);
    polynomial_destroy(&a);
    polynomial_destroy(&b);
    MPI_Finalize();
    return verified ? 0 : 20;
}
