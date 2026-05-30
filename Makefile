HOMEBREW_MPICC := /opt/homebrew/bin/mpicc
HOMEBREW_MPIRUN := /opt/homebrew/bin/mpirun
MPICC ?= $(if $(wildcard $(HOMEBREW_MPICC)),$(HOMEBREW_MPICC),mpicc)
MPIRUN ?= $(if $(wildcard $(HOMEBREW_MPIRUN)),$(HOMEBREW_MPIRUN),mpirun)
CC ?= cc
RM ?= rm -f
MKDIR_P ?= mkdir -p

BUILD_DIR := build
BIN_DIR := bin
SRC_DIR := src
INC_DIR := include
TEST_DIR := tests

TARGET := $(BIN_DIR)/parallel_fft_poly
TEST_FFT := $(BIN_DIR)/test_fft
TEST_POLY := $(BIN_DIR)/test_polynomial
TEST_MPI := $(BIN_DIR)/test_mpi_product

CFLAGS ?= -std=c11 -O2 -Wall -Wextra -pedantic
CPPFLAGS := -I$(INC_DIR)
LDFLAGS ?=
LDLIBS := -lm

COMMON_SRCS := \
	$(SRC_DIR)/complex_number.c \
	$(SRC_DIR)/fft.c \
	$(SRC_DIR)/mpi_fft.c \
	$(SRC_DIR)/polynomial.c \
	$(SRC_DIR)/io_utils.c

APP_SRCS := $(COMMON_SRCS) $(SRC_DIR)/main.c
APP_OBJS := $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(APP_SRCS))

COMMON_OBJS := $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(COMMON_SRCS))
TEST_FFT_OBJ := $(BUILD_DIR)/test_fft.o
TEST_POLY_OBJ := $(BUILD_DIR)/test_polynomial.o
TEST_MPI_OBJ := $(BUILD_DIR)/test_mpi_product.o

.PHONY: all clean run demo test test-serial test-mpi line-count format-directories

all: $(TARGET)

format-directories:
	$(MKDIR_P) $(BUILD_DIR) $(BIN_DIR)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | format-directories
	$(MPICC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/test_%.o: $(TEST_DIR)/test_%.c | format-directories
	$(MPICC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(TARGET): $(APP_OBJS) | format-directories
	$(MPICC) $(CFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(TEST_FFT): $(COMMON_OBJS) $(TEST_FFT_OBJ) | format-directories
	$(MPICC) $(CFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(TEST_POLY): $(COMMON_OBJS) $(TEST_POLY_OBJ) | format-directories
	$(MPICC) $(CFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(TEST_MPI): $(COMMON_OBJS) $(TEST_MPI_OBJ) | format-directories
	$(MPICC) $(CFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

run: $(TARGET)
	$(MPIRUN) -np 4 $(TARGET) --demo 8 --print

demo: $(TARGET)
	$(MPIRUN) -np 4 $(TARGET) --input examples/input_small.txt --output build/output_small.txt --print

test: test-serial test-mpi

test-serial: $(TEST_FFT) $(TEST_POLY)
	$(TEST_FFT)
	$(TEST_POLY)

test-mpi: $(TEST_MPI) $(TARGET)
	$(MPIRUN) -np 4 $(TEST_MPI)
	$(MPIRUN) -np 4 $(TARGET) --input examples/input_small.txt --output build/output_small.txt --quiet

line-count:
	find include src tests scripts -type f \( -name '*.h' -o -name '*.c' -o -name '*.sh' \) -print0 | xargs -0 wc -l

clean:
	$(RM) $(BUILD_DIR)/*.o
	$(RM) $(TARGET) $(TEST_FFT) $(TEST_POLY) $(TEST_MPI)
	$(RM) $(BUILD_DIR)/output_small.txt
