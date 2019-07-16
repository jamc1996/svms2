CPC = g++
CC = gcc

CPFLAGS = -Wall -std=c++11
CFLAGS = -Wall -std=c99

objects = main.o io.o fullproblem.o subproblem.o kernels.o

gert: $(objects)
	$(CC) $(CFLAGS) -o gert $(objects) -lm

main.o: main.c svm.h
	$(CC) $(CFLAGS) -c $<

io.o: io.c io.h svm.h
	$(CC) $(CFLAGS) -c $<

subproblem.o: subproblem.c subproblem.h svm.h
	$(CC) $(CFLAGS) -c $<

fullproblem.o: fullproblem.c fullproblem.h svm.h
	$(CC) $(CFLAGS) -c $<

kernels.o: kernels.c kernels.h svm.h
	$(CC) $(CFLAGS) -c $<

.PHONY: clean
clean:
	rm -f gert $(objects)

test: testPHP testStates testBasketball testCircle

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

testCircle: gert
	./gert -f circleData.txt -t 1 -d 2 -g 1 -s circleModel.txt
	./gert -f circleData.txt -k circleModel.txt -t 1 -d 2 -g 1
