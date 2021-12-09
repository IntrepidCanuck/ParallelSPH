#include <iostream>
#include <math.h>

// This is a program to perform a relatively simple SPH fluid particle simulation
// in a Parallelized manner. This problem simulates the dropping of a block of fluid
// 0.5x0.5 in size from the center of a 2-D domain unit width and unit height

using namespace std;

int main(int argc, char **argv){

// Initialize parameters of the problem which can be defined by the user
int N = 100; // Number of particles

// Do not allow the number of particles to be anything other than a perfect square
// for simplicity
if (ceil((double)sqrt(N)) != floor((double)sqrt(N))){
  printf("ERROR: The number of particles specified must be a perfect square.\n");
}

// Define arrays to store properties of each particle at each time step
double* density = new double[N];
double* pressure = new double[N];
double* x = new double[N]; // x and y location
double* y = new double[N];
double* u = new double[N]; // u and v velocity components
double* v = new double[N];

// Initialize the particles at the center of the domain
int linearN = sqrt(N);
double dx = 0.5/(double)linearN; // Calculate the spacing between each particle
// Loop through every particle and calculate the starting position of each particle within
// the centered 0.5x0.5 box
for (int i = 0;i < N; ++i){
  x[i] = ((i%linearN)/(double)(linearN-1)*0.5)+0.25;
  y[i] = ((i/linearN)/(double)(linearN-1)*0.5)+0.25;
}



printf("%f",dx);


}
