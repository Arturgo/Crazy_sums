BIN := crazysums
#GPROF := -pg
OPT := -O3 -flto
#OPT := -O2

all: build run

build: $(BIN)

-include $(BIN).d

$(BIN): main.cpp Makefile
	g++ -o "$@" $< -Wall -Wextra -std=c++17 $(OPT) -march=native -lpthread -lstdc++fs -MMD -g \
	    ${EXTRA} $(GPROF) -DHAS_COLOR \

run:
	./$(BIN)
	@@echo Generating pdf from LaTeX logs...
	@(cd tex && pdflatex --interaction=batchmode -halt-on-error logs.tex \
	    | { grep -v "This is pdfTeX, " || true; } \
	    | { grep -v "restricted \\\\write18 enabled" || true; } \
	    | { grep -v "entering extended mode" || true; } \
	)

timer: build runT

runT:
	time ./$(BIN)

clean:
	rm -rf $(BIN) $(BIN).d tex
