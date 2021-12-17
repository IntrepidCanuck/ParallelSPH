clear
clc
close all

% Set to 1 if you want to plot gif of results
dispMovie = 1;

% Import and sort data from simulation
imported = importdata("results.txt");
xVals = imported(:,1:2:size(imported,2));
yVals = imported(:,2:2:size(imported,2));

% Store the indexes of the timesteps you wish to plot at save
plotsToSave = [1 2800 4000 6400 10000];

% If the gif is desired, plot every 100 time steps to animate fluid motion
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

% Counter for the number of figure objects
figurecounter = 2;

% Plot the timesteps we want to save
for i = 1:length(plotsToSave)
    figure(figurecounter)
    for j = 1:size(xVals,2)
        xlim([0 1]);
        ylim([0 1]);
        hold on
        plot(xVals(plotsToSave(i),j), yVals(plotsToSave(i),j),'ob');
    end

    % Make the plot look nice
    grid on
    grid minor
    title(strcat("Simulation Results at t = ",string(plotsToSave(i)*0.0001)))
    set(gca,'fontsize',14)

    % Save the figure
    saveas(figure(figurecounter),strcat("Results_",string(plotsToSave(i)),".png"))

    % Increment the figure counter
    figurecounter = figurecounter+1;
end

