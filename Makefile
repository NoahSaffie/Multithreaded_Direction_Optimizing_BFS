CC = g++
CFLAGS = -g -Wall -Wextra -pedantic
OBJECTS = main.o

driver: $(OBJECTS)
	$(CC) $(CFLAGS) -o driver -fopenmp $(OBJECTS)
main.o: main.cpp graph.h graph.hpp wtime.h
	$(CC) $(CFLAGS) -fopenmp -c main.cpp -o main.o
clean:
	rm driver $(OBJECTS)