#ifndef IO_UTILS_H
#define IO_UTILS_H

#include <mpi.h>
#include "polynomial.h"

typedef struct ProgramOptions {
    const char *input_path;
    const char *output_path;
    int demo_degree;
    int verify;
    int quiet;
    int print_vectors;
    int use_serial;
} ProgramOptions;

typedef enum InputStatus {
    INPUT_OK = 0,
    INPUT_ERR_ARGS = 1,
    INPUT_ERR_IO = 2,
    INPUT_ERR_FORMAT = 3,
    INPUT_ERR_MEMORY = 4
} InputStatus;

void program_options_default(ProgramOptions *options);
InputStatus parse_program_options(int argc,
                                  char **argv,
                                  ProgramOptions *options,
                                  int rank);
void print_usage(const char *program_name);

InputStatus read_polynomials_from_file(const char *path,
                                       Polynomial *a,
                                       Polynomial *b);
InputStatus write_polynomial_to_file(const char *path,
                                     const Polynomial *poly);

InputStatus broadcast_polynomial(Polynomial *poly, int root, MPI_Comm comm);
InputStatus broadcast_options(ProgramOptions *options, int root, MPI_Comm comm);

const char *input_status_message(InputStatus status);

#endif
