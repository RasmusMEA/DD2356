.PHONY all: build/main
	./build/main

build/finitevolume.o: src/finitevolume.cpp src/finitevolume.h
	$(CXX) -c $< -o $@

build/matrix.o: src/matrix.cpp src/matrix.h
	$(CXX) -c $< -o $@

build/main: src/main.cpp build/matrix.o build/finitevolume.o
	$(CXX) $^ -o $@

