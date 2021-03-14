BIN := crazysums
#GPROF := -pg
OPT := -O3 -flto
OPT := -O2

all: build run

build: $(BIN)

-include $(BIN).d

$(BIN): main.cpp Makefile
	g++ -o "$@" $< -Wall -Wextra -std=c++17 $(OPT) -march=native -lgmp -MMD -g \
	    ${EXTRA} $(GPROF) -DHAS_COLOR \

run:
	./$(BIN)

timer: build runT

runT:
	time ./$(BIN)
	
clean:
	rm -f $(BIN)