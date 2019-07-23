CPC = g++
CC = gcc
MPCC = mpicc
CPFLAGS = -Wall -std=c++11
CFLAGS = -Wall -std=c99

objects = main.o io.o fullproblem.o subproblem.o kernels.o linked.o

gert: $(objects)
	$(MPCC) $(CFLAGS) -o gert $(objects) -lm

main.o: main.c svm.h
	$(MPCC) $(CFLAGS) -c $<

io.o: io.c io.h svm.h
	$(MPCC) $(CFLAGS) -c $<

subproblem.o: subproblem.c subproblem.h svm.h
	$(MPCC) $(CFLAGS) -c $<

fullproblem.o: fullproblem.c fullproblem.h svm.h
	$(MPCC) $(CFLAGS) -c $<

kernels.o: kernels.c kernels.h svm.h
	$(MPCC) $(CFLAGS) -c $<

linked.o: linked.c linked.h svm.h
	$(MPCC) $(CFLAGS) -c $<

.PHONY: clean
clean:
	rm -f gert $(objects)

test: testPHP testStates testBasketball testPoly

testPHP: gert
	./gert -f phpData.txt -s phpModel.txt
	./gert -f phpData.txt -t phpModel.txt

testStates: gert
	./gert -f alt.txt -s altModel.txt
	./gert -f alt.txt -t altModel.txt

testBasketball: gert
	./gert -f bballData.txt -s bballModel.txt
	./gert -f bballData.txt -t bballModel.txt

testKC1: gert
	./gert -f kc1.csv -s kc1Model.txt
	./gert -f kc1.csv -t kc1Model.txt

testFake: gert
	./gert -f bigData.txt -s bigDataModel.txt
	./gert -f bigData.txt -t bigDataModel.txt

testPoly: gert
	./gert -f circleData.txt -k 1 -d 2 -g 1 -s circleModel.txt
	./gert -f circleData.txt -t circleModel.txt -k 1 -d 2 -g 1

testExp: gert
	./gert -f circleData.txt -k 2 -d 2 -g 1 -s expModel.txt
	./gert -f circleData.txt -t expModel.txt -k 2 -d 2 -g 1

testMP: gert
	module load cports openmpi;\
	mpiexec -n 4 ./gert -f alt.txt -s multiTrain.txt
