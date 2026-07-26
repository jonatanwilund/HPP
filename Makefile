UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S), Darwin)
    CC = gcc-15
    CFLAGS = -Wall -Wextra -O3 -march=native -fopenmp -flax-vector-conversions
    LDFLAGS = -framework Accelerate
else
    CC = gcc
    CFLAGS = -Wall -Wextra -O3 -march=native -fopenmp
    LDFLAGS = -lopenblas -lm
endif

.PHONY: all clean

all: strassen_algorithm strassen_algorithm_old

strassen_algorithm: main.o strassen.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

main.o: main.c strassen.h
	$(CC) $(CFLAGS) -c main.c

strassen.o: strassen.c strassen.h
	$(CC) $(CFLAGS) -c strassen.c

strassen_algorithm_old: strassen_algorithm_old.o
	$(CC) $(CFLAGS) -o $@ $^

strassen_algorithm_old.o: strassen_algorithm_old.c
	$(CC) $(CFLAGS) -c strassen_algorithm_old.c

clean:
	rm -f strassen_algorithm strassen_algorithm_old *.o