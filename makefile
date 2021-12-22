

# Compilation flags
flags = -Wall -std=c++17 `root-config --libs` `root-config --cflags` -lgsl -lgslcblas -Os

# Header file
inc = include/functions.h


all: main error propagate


main: main.o objectDir
	g++ main.o $(flags) -o main
	mv main.o obj

error: error.o objectDir
	g++ error.o $(flags) -o error
	mv error.o obj

propagate: propagate.o objectDir
	g++ propagate.o $(flags) -o propagate
	mv propagate.o obj

objectDir:
	mkdir -p obj


main.o: src/main.cpp $(inc)
	g++ -c src/main.cpp $(inc)

error.o: src/error.cpp $(inc)
	g++ -c src/error.cpp $(inc)

propagate.o: src/propagate.cpp $(inc)
	g++ -c src/propagate.cpp $(inc)

clean:
	rm -rf obj
	rm error main propagate
	rm include/functions.h.gch



