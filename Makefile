CC = gcc
CFLAGS = -Iinclude -Wall -Wextra -O3
LDFLAGS = -lm

all: main

main: src/main.o src/physics.o
	$(CC) $(LDFLAGS) -o main src/main.o src/physics.o

src/main.o: src/main.c include/physics.h
	$(CC) $(CFLAGS) -c src/main.c -o src/main.o

src/physics.o: src/physics.c include/physics.h
	$(CC) $(CFLAGS) -c src/physics.c -o src/physics.o

clean:
	rm -f src/*.o main
