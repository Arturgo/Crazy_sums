VERSION=17
FLAG=-O2


all: build run
	

build:
	g++ -o main main.cpp -Wall -std=c++$(VERSION) $(FLAG) -lgmp

run:
	./main

timer: build runT

runT:
	time ./main