clear
close all
clc

% Import the results data
imported = importdata("results.txt");
xVals = imported(:,1:2:size(imported,2));
yVals = imported(:,2:2:size(imported,2));

% Calculate velocities as difference between subsequent location values
for i = 1:length(xVals)-1
    uVals(i,:) = (xVals(i+1,:) -xVals(i,:))./0.0001;
    vVals(i,:) = (yVals(i+1,:) -yVals(i,:))./0.0001;
end

% Calculate sum of velocity squared
kinEnergy = sum(uVals.^2 + vVals.^2,2);

% Open figure and plot kinetic energy
figure(1)
plot(kinEnergy,'linewidth',2)

% Make plot look nice
title('Kinetic Energy Analog of SPH Simulation System');
ylabel('Mass Independent Squared Velocity (s^{-2})')
xlabel("Time Steps")
set(gca,'fontsize',14)
grid on
grid minor

% Save figure
saveas(figure(1),'kinEnergy.png')
