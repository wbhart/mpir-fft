FFT_LIBS=-L/home/wbhart/mpir-2.4.0/.libs
FFT_INC=-I/home/wbhart/mpir-2.4.0
FFT_FLAGS=-O2 -g

all: mul_fft.c mul_fft.c
	gcc $(FFT_FLAGS) mul_fft.c -o mul_fft $(FFT_INC) $(FFT_LIBS) -lmpir

