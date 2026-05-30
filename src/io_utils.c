#include "io_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void program_options_default(ProgramOptions *options) {
    if (options == NULL) {
        return;
    }
    options->input_path = NULL;
    options->output_path = NULL;
    options->demo_degree = 8;
    options->verify = 1;
    options->quiet = 0;
    options->print_vectors = 0;
    options->use_serial = 0;
}

static int option_is(const char *arg, const char *long_name, const char *short_name) {
    if (arg == NULL) {
        return 0;
    }
    return strcmp(arg, long_name) == 0 || (short_name != NULL && strcmp(arg, short_name) == 0);
}

InputStatus parse_program_options(int argc,
                                  char **argv,
                                  ProgramOptions *options,
                                  int rank) {
    (void)rank;
    if (options == NULL) {
        return INPUT_ERR_ARGS;
    }
    program_options_default(options);

    for (int i = 1; i < argc; ++i) {
        if (option_is(argv[i], "--input", "-i")) {
            if (i + 1 >= argc) {
                return INPUT_ERR_ARGS;
            }
            options->input_path = argv[++i];
        } else if (option_is(argv[i], "--output", "-o")) {
            if (i + 1 >= argc) {
                return INPUT_ERR_ARGS;
            }
            options->output_path = argv[++i];
        } else if (option_is(argv[i], "--demo", "-d")) {
            if (i + 1 >= argc) {
                return INPUT_ERR_ARGS;
            }
            options->demo_degree = atoi(argv[++i]);
            if (options->demo_degree < 0) {
                return INPUT_ERR_ARGS;
            }
        } else if (strcmp(argv[i], "--no-verify") == 0) {
            options->verify = 0;
        } else if (option_is(argv[i], "--quiet", "-q")) {
            options->quiet = 1;
        } else if (strcmp(argv[i], "--print") == 0) {
            options->print_vectors = 1;
        } else if (strcmp(argv[i], "--serial") == 0) {
            options->use_serial = 1;
        } else if (option_is(argv[i], "--help", "-h")) {
            return INPUT_ERR_ARGS;
        } else {
            return INPUT_ERR_ARGS;
        }
    }
    return INPUT_OK;
}

void print_usage(const char *program_name) {
    const char *name = program_name != NULL ? program_name : "parallel_fft_poly";
    printf("Usage:\n");
    printf("  mpirun -np <p> %s [options]\n", name);
    printf("\n");
    printf("Options:\n");
    printf("  -i, --input <file>     Read two polynomials from a text file.\n");
    printf("  -o, --output <file>    Write result polynomial to a text file.\n");
    printf("  -d, --demo <degree>    Generate deterministic demo polynomials.\n");
    printf("      --serial           Run serial FFT instead of MPI FFT.\n");
    printf("      --no-verify        Skip naive O(n^2) verification.\n");
    printf("      --print            Print all coefficient vectors.\n");
    printf("  -q, --quiet            Print only essential lines.\n");
    printf("  -h, --help             Show this help message.\n");
    printf("\n");
    printf("Input format:\n");
    printf("  degree_of_A\n");
    printf("  a1 a2 ... a(n+1)    where a1 is coefficient of x^n\n");
    printf("  degree_of_B\n");
    printf("  b1 b2 ... b(n+1)    where b1 is coefficient of x^n\n");
}

static InputStatus read_one_polynomial(FILE *file, Polynomial *poly) {
    size_t degree = 0;
    if (fscanf(file, "%zu", &degree) != 1) {
        return INPUT_ERR_FORMAT;
    }
    if (!polynomial_init(poly, degree)) {
        return INPUT_ERR_MEMORY;
    }
    for (size_t input_index = 0; input_index < poly->length; ++input_index) {
        size_t power = degree - input_index;
        if (fscanf(file, "%lf", &poly->coefficients[power]) != 1) {
            polynomial_destroy(poly);
            return INPUT_ERR_FORMAT;
        }
    }
    return INPUT_OK;
}

InputStatus read_polynomials_from_file(const char *path,
                                       Polynomial *a,
                                       Polynomial *b) {
    if (path == NULL || a == NULL || b == NULL) {
        return INPUT_ERR_ARGS;
    }
    FILE *file = fopen(path, "r");
    if (file == NULL) {
        return INPUT_ERR_IO;
    }
    InputStatus status = read_one_polynomial(file, a);
    if (status == INPUT_OK) {
        status = read_one_polynomial(file, b);
    }
    fclose(file);
    return status;
}

InputStatus write_polynomial_to_file(const char *path,
                                     const Polynomial *poly) {
    if (path == NULL || poly == NULL || poly->coefficients == NULL) {
        return INPUT_ERR_ARGS;
    }
    FILE *file = fopen(path, "w");
    if (file == NULL) {
        return INPUT_ERR_IO;
    }
    fprintf(file, "%zu\n", poly->degree);
    for (size_t output_index = 0; output_index < poly->length; ++output_index) {
        size_t power = poly->degree - output_index;
        fprintf(file, "%.12g", poly->coefficients[power]);
        if (output_index + 1 < poly->length) {
            fputc(' ', file);
        }
    }
    fputc('\n', file);
    fclose(file);
    return INPUT_OK;
}

InputStatus broadcast_polynomial(Polynomial *poly, int root, MPI_Comm comm) {
    if (poly == NULL) {
        return INPUT_ERR_ARGS;
    }
    int rank = 0;
    MPI_Comm_rank(comm, &rank);

    size_t degree = poly->degree;
    MPI_Bcast(&degree, 1, MPI_UNSIGNED_LONG, root, comm);
    if (rank != root) {
        if (!polynomial_init(poly, degree)) {
            return INPUT_ERR_MEMORY;
        }
    }
    MPI_Bcast(poly->coefficients,
              (int)poly->length,
              MPI_DOUBLE,
              root,
              comm);
    return INPUT_OK;
}

InputStatus broadcast_options(ProgramOptions *options, int root, MPI_Comm comm) {
    if (options == NULL) {
        return INPUT_ERR_ARGS;
    }
    int payload[5];
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    if (rank == root) {
        payload[0] = options->demo_degree;
        payload[1] = options->verify;
        payload[2] = options->quiet;
        payload[3] = options->print_vectors;
        payload[4] = options->use_serial;
    }
    MPI_Bcast(payload, 5, MPI_INT, root, comm);
    if (rank != root) {
        options->demo_degree = payload[0];
        options->verify = payload[1];
        options->quiet = payload[2];
        options->print_vectors = payload[3];
        options->use_serial = payload[4];
    }
    return INPUT_OK;
}

const char *input_status_message(InputStatus status) {
    switch (status) {
        case INPUT_OK:
            return "ok";
        case INPUT_ERR_ARGS:
            return "invalid arguments";
        case INPUT_ERR_IO:
            return "input/output error";
        case INPUT_ERR_FORMAT:
            return "invalid input format";
        case INPUT_ERR_MEMORY:
            return "memory allocation failure";
        default:
            return "unknown error";
    }
}
