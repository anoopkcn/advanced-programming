SRC = $(wildcard *.cc)
CXX = g++-8
CXXFLAGS = -g -std=c++17 -O3
EXE = $(SRC:.cc=.x)

.SUFFIXES:
.SUFFIXES: .cc .x

all: $(EXE)

%.x: %.cc
	$(CXX) $< -o $@ $(CXXFLAGS)


format: $(SRC)
	@clang-format -i $^ 2>/dev/null || echo "Please install clang-format to run this commands"

clean:
	rm -f $(EXE) *~ 
	rm -rf *.*.dSYM 

.PHONY: clean
