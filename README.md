# Parallel FFT Polynomial Multiplication

Project này cài đặt thuật toán nhân hai đa thức bằng Cooley-Tukey FFT.
Phiên bản MPI hiện tại sử dụng hướng **parallel divide-and-conquer**:
mỗi tầng chia bài toán FFT thành hai nửa `even` và `odd`, chia communicator
MPI bằng `MPI_Comm_split`, sau đó ghép kết quả bằng butterfly.

## Cấu Trúc Folder

```text
.
├── Makefile
├── README.md
├── examples/
│   ├── input_small.txt
│   └── input_degree8.txt
├── include/
│   ├── complex_number.h
│   ├── fft.h
│   ├── io_utils.h
│   ├── mpi_fft.h
│   └── polynomial.h
├── scripts/
│   ├── count_code_lines.sh
│   ├── run_demo.sh
│   └── run_file.sh
├── src/
│   ├── complex_number.c
│   ├── fft.c
│   ├── io_utils.c
│   ├── main.c
│   ├── mpi_fft.c
│   └── polynomial.c
└── tests/
    ├── test_fft.c
    ├── test_mpi_product.c
    └── test_polynomial.c
```

## Yêu Cầu

Cần có:

```text
make
OpenMPI
mpicc
mpirun
```

Trên máy hiện tại, nên dùng OpenMPI của Homebrew:

```sh
/opt/homebrew/bin/mpicc
/opt/homebrew/bin/mpirun
```

Nếu `which mpirun` đang trỏ tới Anaconda, ví dụ:

```text
/opt/anaconda3/bin/mpirun
```

thì nên chạy bằng đường dẫn đầy đủ:

```sh
/opt/homebrew/bin/mpirun ...
```

hoặc đưa Homebrew lên trước trong `PATH`:

```sh
export PATH="/opt/homebrew/bin:$PATH"
```

## Build

Biên dịch chương trình chính:

```sh
make all
```

File executable được tạo tại:

```text
bin/parallel_fft_poly
```

Xoá file build:

```sh
make clean
```

## Chạy Với File Input

Ví dụ:

```sh
/opt/homebrew/bin/mpirun -np 4 ./bin/parallel_fft_poly \
  --input examples/input_small.txt \
  --output build/output_small.txt \
  --print
```

Hoặc dùng script:

```sh
scripts/run_file.sh 4 examples/input_small.txt build/output_small.txt
```

Trong đó:

```text
-np 4                  chạy với 4 process MPI
--input <file>          file chứa hai đa thức A và B
--output <file>         file ghi đa thức kết quả C
--print                 in A(x), B(x), C(x) ra terminal
```

## Chạy Với Dữ Liệu Demo

Tự sinh hai đa thức bậc `16`:

```sh
/opt/homebrew/bin/mpirun -np 4 ./bin/parallel_fft_poly --demo 16 --print
```

Hoặc dùng script:

```sh
scripts/run_demo.sh 4 16
```

## Chạy Bản Serial

Để chạy FFT tuần tự, thêm option `--serial`:

```sh
/opt/homebrew/bin/mpirun -np 1 ./bin/parallel_fft_poly \
  --demo 100000 \
  --serial \
  --quiet \
  --no-verify
```

So sánh với bản MPI divide-and-conquer:

```sh
/opt/homebrew/bin/mpirun -np 4 ./bin/parallel_fft_poly \
  --demo 100000 \
  --quiet \
  --no-verify
```

Lưu ý: dùng `--no-verify` khi degree lớn vì verify naive có độ phức tạp
`O(n^2)`.

## Format Input

Input gồm hai đa thức. Với đa thức bậc `n`, hệ số được viết theo thứ tự
từ bậc cao xuống bậc thấp:

```text
degree_of_A
a1 a2 ... a(n+1)
degree_of_B
b1 b2 ... b(n+1)
```

Ý nghĩa:

```text
A(x) = a1*x^n + a2*x^(n-1) + ... + a(n+1)
B(x) = b1*x^n + b2*x^(n-1) + ... + b(n+1)
```

Ví dụ:

```text
3
1 2 3 4
3
5 6 7 8
```

tương ứng:

```text
A(x) = 1x^3 + 2x^2 + 3x + 4
B(x) = 5x^3 + 6x^2 + 7x + 8
```

## Format Output

Output cũng ghi hệ số theo thứ tự từ bậc cao xuống bậc thấp:

```text
degree_of_C
c1 c2 ... c(m+1)
```

Ví dụ output:

```text
6
5 16 34 60 61 52 32
```

tương ứng:

```text
C(x) = 5x^6 + 16x^5 + 34x^4 + 60x^3 + 61x^2 + 52x + 32
```

## Test

Chạy toàn bộ test:

```sh
make test
```

Chỉ chạy test tuần tự:

```sh
make test-serial
```

Chỉ chạy test MPI:

```sh
make test-mpi
```

## Đếm Dòng Code

```sh
make line-count
```

hoặc:

```sh
scripts/count_code_lines.sh
```

## Các Option Chính

```text
--input <file>      đọc input từ file
--output <file>     ghi output ra file
--demo <degree>     sinh input demo với đa thức bậc degree
--serial            chạy FFT tuần tự
--no-verify         bỏ kiểm tra bằng nhân naive O(n^2)
--quiet             chỉ in thông tin ngắn
--print             in A(x), B(x), C(x)
--help              in hướng dẫn CLI
```

## Ghi Chú Về MPI

Nếu chạy:

```sh
mpirun -np 4 ...
```

mà chương trình in kết quả 4 lần và báo `processors: 1`, nguyên nhân thường
là `mpirun` không cùng bản MPI với `mpicc` đã dùng để build. Khi đó hãy dùng:

```sh
/opt/homebrew/bin/mpirun -np 4 ...
```

thay vì `mpirun` của Anaconda.
