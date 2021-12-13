clear
clc
close all

imported = importdata("results.txt");
xVals = imported(:,1:2:size(imported,2));
yVals = imported(:,2:2:size(imported,2));


figure(1)
hold on
for i = 1:100:length(imported)
    for j = 1:size(xVals,2)
        xlim([0 1]);
        ylim([0 1]);
        hold on
        plot(xVals(i,j), yVals(i,j),'ob');
    end
    pause(.1)
    clf('reset')
end