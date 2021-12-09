#include <iostream>
#include <fstream>
#include <math.h>

// This is a program to perform a relatively simple SPH fluid particle simulation
// in a Parallelized manner. This problem simulates the dropping of a block of fluid
// 0.5x0.5 in size from the center of a 2-D domain unit width and unit height

using namespace std;

double kernelCalc(double x1, double y1, double x2, double y2, double rad){
  double dist = (sqrt(pow(x1-x2,2)+pow(y1-y2,2)))/rad;

  if (dist < 1.0){
    return (4.0/3.14159/pow(rad,2))*pow(1-pow(dist,2),3);
  } else {
    return 0.0;
  }
}

double FirstOrderKernelCalcX(double x1, double y1, double x2, double y2, double rad){
  double dist = (sqrt(pow(x1-x2,2)+pow(y1-y2,2)))/rad;

  if (dist < 1.0){
    return (-30/3.14159/pow(rad,3))*(x1-x2)*(pow(1-dist,2)/dist);
  } else {
    return 0.0;
  }
}

double FirstOrderKernelCalcY(double x1, double y1, double x2, double y2, double rad){
  double dist = (sqrt(pow(x1-x2,2)+pow(y1-y2,2)))/rad;

  if (dist < 1.0){
    return (-30/3.14159/pow(rad,3))*(y1-y2)*(pow(1-dist,2)/dist);
  } else {
    return 0.0;
  }
}

double SecondOrderKernelCalcX(double x1, double y1, double x2, double y2, double u1, double u2, double rad){
  double dist = (sqrt(pow(x1-x2,2)+pow(y1-y2,2)))/rad;

  if (dist < 1.0){
    return (12/pow(rad,4))*(u1 - u2)*(1-dist);
  } else {
    return 0.0;
  }
}

double SecondOrderKernelCalcY(double x1, double y1, double x2, double y2, double v1, double v2, double rad){
  double dist = (sqrt(pow(x1-x2,2)+pow(y1-y2,2)))/rad;

  if (dist < 1.0){
    return (12/pow(rad,4))*(v1-v2)*(1-dist);
  } else {
    return 0.0;
  }
}

int main(int argc, char **argv){

// Open up a file to which we can output data
string file_name = "jacobi_" + to_string(rank) + ".txt";
ofstream fp(file_name);

// Initialize parameters of the problem which can be defined by the user
int N = 900; // Number of particles
int mass = 1; // Initial guess for mass (will be updated later)
double radInfluence = 1; // Radius of Influence of each particle
double restingDensity = 1000; // Resting density of the fluid
double T = 1; // Total time of simulation
double dt = 0.01; // Time step
double pressureConstant = 2000; // Pressure constant
double viscConstant = 1.0; // Viscosity Constant
double gravConstant = 9.8; // Gravitational Constant
double restitutionConst = 0.5; // Restitution Coeffecient

// Do not allow the number of particles to be anything other than a perfect square
// for simplicity
if (ceil((double)sqrt(N)) != floor((double)sqrt(N))){
  printf("ERROR: The number of particles specified must be a perfect square.\n");
}

// Define arrays to store properties of each particle at each time step
double* density = new double[N];
double* pressure = new double[N];
double* pressureX = new double[N]; // Pressure forces in the x and y direction
double* pressureY = new double[N];
double* viscosityX = new double[N]; // Viscosity forces in the x and y direction
double* viscosityY = new double[N];
double* gravY = new double[N]; // Gravity force
double* x = new double[N]; // x and y location
double* y = new double[N];
double* u = new double[N]; // u and v velocity components
double* v = new double[N];

// Initialize the particles at the center of the domain
int linearN = sqrt(N);
double dx = 0.5/(double)linearN; // Calculate the spacing between each particle
printf("dx = %f\n",dx);
// Loop through every particle and calculate the starting position of each particle within
// the centered 0.5x0.5 box
for (int i = 0;i < N; ++i){
  x[i] = ((i%linearN)/(double)(linearN-1)*0.5)+0.25;
  y[i] = ((i/linearN)/(double)(linearN-1)*0.5)+0.25;
}

// Initialize the velocities as stationary
for (int i = 0; i < N; ++i){
  u[i] = 0.0;
  v[i] = 0.0;
}

// We first need to determine the mass which each particle represents
double densityContribution;
for (int i = 0; i < N; ++i){
  densityContribution = 0.0;
  for (int j = 0; j < N; ++j){
    densityContribution += mass*kernelCalc(x[i],y[i],x[j],y[j],radInfluence);
  }
  density[i] = densityContribution;
}

// Now recalculate the mass, assuming mass is spread evenly among all particles
double densitySum = 0;
for (int i = 0; i < N; ++i){
  densitySum += density[i];
}
mass = ((double)N*restingDensity)/densitySum;

// Counter for the current timestep of the simulation
double currentTime = 0.0;

while (currentTime <= T){

  printf("%f \n",currentTime);

  // Perform the density calculation
  for (int i = 0; i < N; ++i){
    densityContribution = 0.0;
    for (int j = 0; j < N; ++j){
      densityContribution += mass*kernelCalc(x[i],y[i],x[j],y[j],radInfluence);
    }
    density[i] = densityContribution;
  }

  // Perform the Pressure calculation
  for (int i = 0; i < N; ++i){
    pressure[i] = pressureConstant*(density[i] - restingDensity);
  }

  // Perform the Pressure force Calculation
  double pressureContributionX;
  double pressureContributionY;
  for (int i = 0; i < N; ++i){
    pressureContributionX = 0.0;
    pressureContributionY = 0.0;
    for (int j = 0; j < N; ++j){
      pressureContributionX += (mass/density[j])*(pressure[i]+pressure[j])*(1/2)*FirstOrderKernelCalcX(x[i],y[i],x[j],y[j],radInfluence);
      pressureContributionY += (mass/density[j])*(pressure[i]+pressure[j])*(1/2)*FirstOrderKernelCalcX(x[i],y[i],x[j],y[j],radInfluence);
    }
    pressureX[i] = -1.0*pressureContributionX;
    pressureY[i] = -1.0*pressureContributionY;
  }

  // Perform the viscous force calculations
  double viscousContributionX;
  double viscousContributionY;
  for (int i = 0; i < N; ++i){
    viscousContributionX = 0.0;
    viscousContributionY = 0.0;
    for (int j = 0; j < N; ++j){
      viscousContributionX += (mass/density[j])*SecondOrderKernelCalcX(x[i],y[i],x[j],y[j],u[i],u[j],radInfluence);
      viscousContributionY += (mass/density[j])*SecondOrderKernelCalcY(x[i],y[i],x[j],y[j],v[i],v[j],radInfluence);
    }
    viscosityX[i] = -1.0*viscConstant*viscousContributionX;
    viscosityY[i] = -1.0*viscConstant*viscousContributionY;
  }

  // Perform the gravity force calculations
  for (int i = 0; i < N; ++i){
    gravY[i] = -1.0*density[i]*gravConstant;
  }

  // Update the velocity
  for (int i = 0; i < N; ++i){
    u[i] += (pressureX[i]+viscosityX[i])*dt/density[i];
    v[i] += (pressureY[i]+viscosityY[i]+gravY[i])*dt/density[i];
  }

  // Update the position of each particle
  for (int i = 0; i < N; ++i){
    x[i] += u[i]*dt;
    y[i] += v[i];
  }

  // Check for the boundary Conditions
  for (int i = 0; i < N; ++i){
    // Check bottom boundary condition
    if (y[i] <= radInfluence){
      y[i] = radInfluence;
      v[i] = -1.0*restitutionConst*v[i];
    }

    // Check top boundary Conditions
    if  (y[i] >= 1-radInfluence){
      y[i] = 1 - radInfluence;
      v[i] = -1.0*restitutionConst*v[i];
    }

    // Check left boundary conditions
    if (x[i] <= radInfluence){
      x[i] = radInfluence;
      u[i] = -1.0*restitutionConst*u[i];
    }

    // Check right boundary conditions
    if (x[i] >= 1-radInfluence){
      x[i] = 1-radInfluence;
      u[i] = -1.0*restitutionConst*u[i];
    }
  }


// Increment the counter
currentTime += dt;





}







}
