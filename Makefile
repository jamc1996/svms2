CPC = g++
CC = g

CPFLAGS = -Wall -std=c++11
CFLAGS = -Wall

objects = main.o io.o

gert: $(objects)
	$(CPC) $(CPFLAGS) -o gert $(objects)

main.o: main.cpp svm.h
	$(CPC) $(CPFLAGS) -c $<

io.o: io.cpp io.h
	$(CPC) $(CPFLAGS) -c $<

.PHONY: clean
clean:
	rm -f gert $(objects)

test: gert
	./gert -f alt.txt
