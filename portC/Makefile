CPC = g++
CC = gcc

CPFLAGS = -Wall -std=c++11
CFLAGS = -Wall

objects = main.o io.o fullproblem.o subproblem.o

gert: $(objects)
	$(CC) $(CFLAGS) -o gert $(objects)

main.o: main.c svm.h
	$(CC) $(CFLAGS) -c $<

io.o: io.c io.h svm.h
	$(CC) $(CFLAGS) -c $<

subproblem.o: subproblem.c subproblem.h
	$(CC) $(CFLAGS) -c $<

fullproblem.o: fullproblem.c fullproblem.h
	$(CC) $(CFLAGS) -c $<

.PHONY: clean
clean:
	rm -f gert $(objects)

test: gert
	./gert -f alt.txt
