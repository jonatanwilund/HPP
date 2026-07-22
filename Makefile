CC = gcc
CFLAGS = -Wall -O3 -march=native -ffast-math -funroll-loops -fopenmp

all: strassen_algorithm strassen_algorithm_old

strassen_algorithm: strassen_algorithm.o
	$(CC) $(CFLAGS) -o strassen_algorithm strassen_algorithm.o

strassen_algorithm.o: strassen_algorithm.c
	$(CC) $(CFLAGS) -c strassen_algorithm.c

strassen_algorithm_old: strassen_algorithm_old.o
	$(CC) $(CFLAGS) -o strassen_algorithm_old strassen_algorithm_old.o

strassen_algorithm_old.o: strassen_algorithm_old.c
	$(CC) $(CFLAGS) -c strassen_algorithm_old.c

clean:
	rm -f strassen_algorithm strassen_algorithm_old *.o