clear
close all
clc

compTime = [594 318 247 201 224 261 262 272];

figure(1)
bar(compTime)

xlabel("Number of Processes")
ylabel("Computation Time (s)")
title("Computation Time vs Number of Processes")
set(gca,'fontsize',14)
grid on
grid minor

saveas(figure(1),'scalability.png')