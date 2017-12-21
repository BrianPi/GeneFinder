# Makes and runs the example code.

all: prog test

prog: testor

testor: testor.cpp *.h
	g++ -std=c++11 -o testor testor.cpp

test:
	./testor
