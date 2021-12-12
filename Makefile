default: main

main.o: main.cpp
	g++-11 -c main.cpp

main: main.o
	g++-11 -o main main.o 

clean:
	rm *.o *.txt main
