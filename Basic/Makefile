CPC = g++
CC = g

CPFLAGS = -Wall -std=c++11
CFLAGS = -Wall

objects = main.o io.o subproblem.o fullproblem.o

gert: $(objects)
	$(CPC) $(CPFLAGS) -o gert $(objects)

main.o: main.cpp svm.h
	$(CPC) $(CPFLAGS) -c $<

io.o: io.cpp io.h svm.h
	$(CPC) $(CPFLAGS) -c $<

subproblem.o: subproblem.cpp subproblem.h
	$(CPC) $(CPFLAGS) -c $<

fullproblem.o: fullproblem.cpp fullproblem.h
	$(CPC) $(CPFLAGS) -c $<

.PHONY: clean
clean:
	rm -f gert $(objects)

test: gert
	./gert -f alt.txt
