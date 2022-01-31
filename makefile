

# Compilation flags
flags = -Wall -std=c++17 `root-config --libs` `root-config --cflags` -lgsl -lgslcblas -Os

# Header file
inc = include/functions.h


all: main error propagate create_data


main: main.o objectDir
	g++ main.o $(flags) -o main
	mv main.o obj

error: error.o objectDir
	g++ error.o $(flags) -o error
	mv error.o obj

propagate: propagate.o objectDir
	g++ propagate.o $(flags) -o propagate
	mv propagate.o obj

create_data: create_data.o objectDir
	g++ create_data.o $(flags) -o cdata
	mv create_data.o obj

objectDir:
	mkdir -p obj


main.o: src/main.cpp $(inc)
	g++ -c src/main.cpp $(inc)

error.o: src/error.cpp $(inc)
	g++ -c src/error.cpp $(inc)

propagate.o: src/propagate.cpp $(inc)
	g++ -c src/propagate.cpp $(inc)

create_data.o: src/create_data.cpp $(inc)
	g++ -c src/create_data.cpp $(inc)

clean:
	rm -rf obj
	rm error main propagate cdata
	rm dataLMDV.txt dataQ4.txt dataQ6.txt
	rm include/functions.h.gch



