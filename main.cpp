#include <iostream>
#include <fstream>
#include <math.h>
#include <ctime>

// This is a program to perform a relatively simple SPH fluid particle simulation
// in a Parallelized manner. This problem simulates the dropping of a block of fluid
// 0.5x0.5 in size from the center of a 2-D domain unit width and unit height

using namespace std;

double kernelCalc(int i, int j, double x1, double y1, double x2, double y2, double rad){

  if (i == j){
    return (4.0/(3.14159*pow(rad,2)));
  }
  double dist = (sqrt(pow(x1-x2,2)+pow(y1-y2,2)))/rad;

  if (dist < 1.0){
    return (4.0/(3.14159*pow(rad,2)))*pow(1-pow(dist,2),3);
  } else {
    return 0.0;
  }
}

double FirstOrderKernelCalcX(double x1, double y1, double x2, double y2, double rad){
  double dist = (sqrt(pow(x1-x2,2)+pow(y1-y2,2)))/rad;

  if (dist < 1.0){
    return (-10/pow(rad,3))*(x1-x2)*(pow(1-dist,2)/dist);
  } else {
    return 0.0;
  }
}

double FirstOrderKernelCalcY(double x1, double y1, double x2, double y2, double rad){
  double dist = (sqrt(pow(x1-x2,2)+pow(y1-y2,2)))/rad;

  if (dist < 1.0){
    return (-10/pow(rad,3))*(y1-y2)*(pow(1-dist,2)/dist);
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

// Seed the random number generator but only change seed once a day so we can compare runs
srand(int(time(0))/86400);

// Open up a file to which we can output data
string file_name = "results.txt";
ofstream fp(file_name);

// Initialize parameters of the problem which can be defined by the user
int N = 225; // Number of particles
double mass = 1; // Initial guess for mass (will be updated later)
double radInfluence = .01; // Radius of Influence of each particle
double restingDensity = 1000; // Resting density of the fluid
double T = 1; // Total time of simulation
double dt = 0.0001; // Time step
double pressureConstant = 2000; // Pressure constant
double viscConstant = 1.0; // Viscosity Constant
double gravConstant = 9.8; // Gravitational Constant
double restitutionConst = 0.5; // Restitution Coeffecient

// Do not allow the number of particles to be anything other than a perfect square
// for simplicity
if (ceil((double)sqrt(N)) != floor((double)sqrt(N))){
  printf("ERROR: The number of particles specified must be a perfect square.\n");
  return 1;
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

double pressureContributionX;
double pressureContributionY;
double viscousContributionX;
double viscousContributionY;

// Initialize the particles at the center of the domain
int linearN = sqrt(N);
double dx = 0.5/(double)linearN; // Calculate the spacing between each particle
printf("dx = %f\n",dx);
// Loop through every particle and calculate the starting position of each particle within
// the centered 0.5x0.5 box
for (int i = 0;i < N; ++i){
  // Insert some randomness to ensure that the particles don't just bounce on top of each other
  double randFactor = ((rand()/double(RAND_MAX)) - 0.5)*(0.1*radInfluence);

  x[i] = (((i%linearN)/(double)(linearN-1)*0.5)+0.25)+randFactor;
  y[i] = (((i/linearN)/(double)(linearN-1)*0.5)+0.25)+randFactor;
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
    densityContribution += mass*kernelCalc(i,j,x[i],y[i],x[j],y[j],radInfluence);
  }
  density[i] = densityContribution;

  //printf("Density of Particle %i = %f\n",i,density[i]);
}

// Now recalculate the mass, assuming mass is spread evenly among all particles
double densitySum = 0;
for (int i = 0; i < N; ++i){
  densitySum += density[i];
}
mass = ((double)N*restingDensity)/densitySum;

//printf("Mass before Loop = %f\n",mass);

// Counter for the current timestep of the simulation
double currentTime = 0.0;

while (currentTime <= T){

  /*if (currentTime%(T/20) == 0.0){
    printf("Simulation %f \% complete \n",currentTime/(T/20));
  }*/

  printf("New TimeStep at %f\n",currentTime);
  //printf("Mass in Loop = %f\n",mass);

  // DEBUG USE only
  /*printf("BEG OF LOOP X\n");
  for (int i = 0; i < N; ++i){
    printf("%f ",x[i]);
  }
  printf("\n");*/


  // Perform the density calculation
  for (int i = 0; i < N; ++i){
    densityContribution = 0.0;
    for (int j = 0; j < N; ++j){
      densityContribution = densityContribution + mass*kernelCalc(i,j,x[i],y[i],x[j],y[j],radInfluence);
      //printf("Density Contribution = %f with mass = %f and %i,%i\n",densityContribution,mass,i,j);
      //printf("     xi=%f,yi=%f,xj=%f,yj=%f\n",x[i],y[i],x[j],y[j]);
      //printf("     Returned kernel val = %f\n",kernelCalc(i,j,x[i],y[i],x[j],y[j],radInfluence));
    }
    density[i] = densityContribution;
    //printf("Density of Particle %i = %f\n",i,density[i]);
  }

  // Perform the Pressure calculation
  for (int i = 0; i < N; ++i){
    pressure[i] = pressureConstant*(density[i] - restingDensity);
  }

  // DEBUG USE only
  /*for (int i = 0; i < N; ++i){
    printf("%f ",x[i]);
  }*/

  // Perform the Pressure force Calculation
  for (int i = 0; i < N; ++i){
    pressureContributionX = 0.0;
    pressureContributionY = 0.0;
    for (int j = 0; j < N; ++j){
      if (i != j){
        pressureContributionX += (mass/density[j])*(pressure[i]+pressure[j])*(1.0/2.0)*FirstOrderKernelCalcX(x[i],y[i],x[j],y[j],radInfluence);
        pressureContributionY += (mass/density[j])*(pressure[i]+pressure[j])*(1.0/2.0)*FirstOrderKernelCalcY(x[i],y[i],x[j],y[j],radInfluence);

        /*printf("densityj = %f, pressurei = %f, pressurej = %f, kernel = %f, at %i,%i\n",density[j],pressure[i],pressure[j],FirstOrderKernelCalcX(x[i],y[i],x[j],y[j],radInfluence),i,j);
        printf("Using x1=%f,y1=%f,x2=%f,y2=%f\n",x[i],y[i],x[j],y[j]);

        if (isnan(pressureContributionX)){
          printf("ISNAN Pressure Contribution found at %i,%i \n",i,j);
          return 1;
        }*/
      }
    }
    pressureX[i] = -1.0*pressureContributionX;
    pressureY[i] = -1.0*pressureContributionY;

    /*if (pressureY[i] != 0.0){
      printf("Pressure is not zero\n");
    }*/

    //printf("PressY = %f\n",pressureY[i]);
  }

  // Perform the viscous force calculations
  for (int i = 0; i < N; ++i){
    viscousContributionX = 0.0;
    viscousContributionY = 0.0;
    for (int j = 0; j < N; ++j){
      if (i != j){
        viscousContributionX += (mass/density[j])*SecondOrderKernelCalcX(x[i],y[i],x[j],y[j],u[i],u[j],radInfluence);
        viscousContributionY += (mass/density[j])*SecondOrderKernelCalcY(x[i],y[i],x[j],y[j],v[i],v[j],radInfluence);
      }
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

    //printf("u = %f\n",u[i]);
    //printf("v = %f\n",v[i]);
  }

  // Update the position of each particle
  for (int i = 0; i < N; ++i){
    x[i] += u[i]*dt;
    y[i] += v[i]*dt;
    /*printf("x = %f\n",x[i]);
    printf("y = %f\n",y[i]);*/
  }

  // DEBUG USE only
  /*printf("END OF LOOP X\n");
  for (int i = 0; i < N; ++i){
    printf("%f ",x[i]);
  }
  printf("\n");*/

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

  // DEBUG USE only
  /*printf("END OF BC X\n");
  for (int i = 0; i < N; ++i){
    printf("%f ",x[i]);
  }
  printf("\n");*/

// Increment the counter
currentTime += dt;

// Output positions to file
for (int i = 0; i < N; ++i){
  fp << x[i] << " " << y[i] << " ";
}

// Start a new line in the output for the next time step
fp << "\n";



}




// Close the file
fp.close();


}
