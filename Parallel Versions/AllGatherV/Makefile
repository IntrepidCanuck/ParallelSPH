default: main

main.o: main.cpp
	mpicxx -c main.cpp

main: main.o
	mpicxx -o main main.o 

clean:
	rm *.o *.txt main
