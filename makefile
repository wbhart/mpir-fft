FFT_LIBS=-L/home/wbhart/mpir-2.1.1/.libs
FFT_INC=-I/home/wbhart/mpir-2.1.1
GMP_LIBS=-L/home/wbhart/gmp-5.0.2/.libs
GMP_INC=-I/home/wbhart/gmp-5.0.2
FFT_FLAGS=-O2 -g

all: mul_fft.c mul_fft.c
	gcc $(FFT_FLAGS) mul_fft.c -o mul_fft $(FFT_INC) $(FFT_LIBS) -lmpir

time_gmp: time_gmp.c
	gcc $(FFT_FLAGS) time_gmp.c -o time_gmp $(GMP_INC) $(GMP_LIBS) -static -lgmp
