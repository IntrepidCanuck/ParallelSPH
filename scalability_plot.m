clear
close all
clc

% Store the computation times found for each number of processes
compTime = [594 318 247 201 224 261 262 272];

% Open a figure and make the bar plot
figure(1)
bar(compTime)

% Make the plot look nice
xlabel("Number of Processes")
ylabel("Computation Time (s)")
title("Computation Time vs Number of Processes")
set(gca,'fontsize',14)
grid on
grid minor

% Save the figure
saveas(figure(1),'scalability.png')