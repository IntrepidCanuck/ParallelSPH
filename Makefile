default: simulation

simulation.o: simulation.cpp
	mpicxx -c simulation.cpp

simulation: simulation.o
	mpicxx -o simulation simulation.o 

clean:
	rm *.o *.txt *.png simulation
