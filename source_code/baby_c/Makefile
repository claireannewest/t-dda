CC=gcc
CFLAGS=-c -Wall -iquote/usr/local/lib/include

all: BABY

BABY: baby_code.o
	$(CC) baby_code.o -o baby_code -L/usr/local/bin -lgfortran -llapack -lblas -lm

baby_code.o: baby_code.c
	$(CC) $(CFLAGS) baby_code.c

clean:
	rm *.o baby_code
