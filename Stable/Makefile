CPC = g++
CC = gcc

CPFLAGS = -Wall -std=c++11
CFLAGS = -Wall -std=c99

objects = main.o io.o fullproblem.o subproblem.o

gert: $(objects)
	$(CC) $(CFLAGS) -o gert $(objects) -lm

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

test: testPHP testStates testBasketball

testPHP: gert
	./gert -f phpVDlhKL.csv -s phpModel.txt
	./gert -f phpVDlhKL.csv -k phpModel.txt

testStates: gert
	./gert -f alt.txt -s altModel.txt
	./gert -f alt.txt -k altModel.txt

testBasketball: gert
	./gert -f basketball.csv -s bball.txt
	./gert -f basketball.csv -k bball.txt

testKC1: gert
	./gert -f kc1.csv -s kc1Model.txt
	./gert -f kc1.csv -k kc1Model.txt

testFake: gert
	./gert -f bigData.txt -s bigDataModel.txt
	./gert -f bigData.txt -k bigDataModel.txt
