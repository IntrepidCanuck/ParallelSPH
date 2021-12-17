clear
clc
close all

dispMovie = 0;

imported = importdata("results.txt");
xVals = imported(:,1:2:size(imported,2));
yVals = imported(:,2:2:size(imported,2));

plotsToSave = [1 2800 4000 6400 10000];

if (dispMovie == 1)
    figure(1)
    hold on
    for i = 1:400:length(imported)
        for j = 1:size(xVals,2)
            xlim([0 1]);
            ylim([0 1]);
            hold on
            plot(xVals(i,j), yVals(i,j),'ob');
        end
        pause(.1)
        clf('reset')
        disp(i)
    end
end

figurecounter = 2;

for i = 1:length(plotsToSave)
    figure(figurecounter)
    for j = 1:size(xVals,2)
        xlim([0 1]);
        ylim([0 1]);
        hold on
        plot(xVals(plotsToSave(i),j), yVals(plotsToSave(i),j),'ob');
    end

    grid on
    grid minor
    title(strcat("Simulation Results at t = ",string(plotsToSave(i)*0.0001)))
    set(gca,'fontsize',14)

    saveas(figure(figurecounter),strcat("Results_",string(plotsToSave(i)),".png"))

    figurecounter = figurecounter+1;
end

