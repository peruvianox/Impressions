%% CreateCoPNorms

% this script will input a set of pre-processed impressions data and create
% a normals database of those results. 

%% select data to load
folder = uigetdir; % select folder
addpath(genpath(folder));
Contents = dir(folder); 
g = 1; 
for i = 1:length(Contents) % loop through main folder
    if length(Contents(i).name) > 3 % if folder has over 3 characters
        SubFolder = strcat(Contents(i).folder,'\', Contents(i).name); % combine name into string
        SubContents = dir(SubFolder); 
        for j = 1:length(SubContents)
            if strfind(SubContents(j).name, 'Results.mat')
                clc; disp({'Loading',SubContents(j).name}); 
                Results = load(SubContents(j).name);
                Data(g).Norm = Results.Norm;
                 g = g+1; % counter for data output
                close all % close all figures that pop up
                break % break inner for loop because data is already loaded for that folder
            end
        end
    end
end
clearvars i j g Results

%% Combine all the COPs into a 3D matrix
for i = 1:length(Data)
    kL = 1;
    kR = 1;
    for j = 1:length(Data(i).Norm)
        if isempty(Data(i).Norm(j).LCoPs) == 0
            Data(i).All_LCoPs(:,:,kL) =  Data(i).Norm(j).LCoPs; 
            Data(i).All_LSums(:,:,kL) =  Data(i).Norm(j).LSum;
            kL = kL + 1;
        end
         if isempty(Data(i).Norm(j).RCoPs) == 0
            Data(i).All_RCoPs(:,:,kR) =  Data(i).Norm(j).RCoPs; 
            Data(i).All_RSums(:,:,kL) =  Data(i).Norm(j).RSum;
            kR = kR + 1; 
         end
    end
end
clearvars kL kR i j SubFolder SubContents folder

%% Calculate averages and SDs for each participant
for i = 1:length(Data)
    Data(i).LMeanCoP = mean( Data(i).All_LCoPs, 3); 
    Data(i).LStdCoP = std( Data(i).All_LCoPs, 0, 3); 
    Data(i).RMeanCoP = mean( Data(i).All_RCoPs, 3);
    Data(i).RStdCoP = std( Data(i).All_RCoPs, 0, 3); 
    
    All_LMeanCoP(:,:,i) = Data(i).LMeanCoP; 
    All_LStdCoP(:,:,i) = Data(i).LStdCoP; 
    All_RMeanCoP(:,:,i) = Data(i).RMeanCoP;
    All_RStdCoP(:,:,i) = Data(i).RStdCoP; 
end

%% calculate group averages and group average of the SD
for i = 1:length(Data)
    LMeanCoP = mean(All_LMeanCoP, 3); 
    LStdCoP = std(All_LStdCoP, 0, 3); 
    RMeanCoP = mean(All_RMeanCoP, 3);
    RStdCoP = std(All_RStdCoP, 0, 3); 
end

%% plot CoP norms
figure;
subplot(121); hold on; 
plot(LMeanCoP(:,1), LMeanCoP(:,2), '-k'); 
plot(LMeanCoP(:,1) + LStdCoP(:,1), LMeanCoP(:,2) + LStdCoP(:,2), '--k'); 
plot(LMeanCoP(:,1) - LStdCoP(:,1), LMeanCoP(:,2) - LStdCoP(:,2), '--k'); 
xlim([20 80]); 
xlabel('% of foot width'); 
ylabel('% of foot length'); 
title('Left CoP Norms'); 

subplot(122); hold on; 
plot(RMeanCoP(:,1), RMeanCoP(:,2), '-k'); 
plot(RMeanCoP(:,1) + RStdCoP(:,1), RMeanCoP(:,2) + RStdCoP(:,2), '--k'); 
plot(RMeanCoP(:,1) - RStdCoP(:,1), RMeanCoP(:,2) - RStdCoP(:,2), '--k'); 
xlim([20 80]); 
xlabel('% of foot width'); 
ylabel('% of foot length'); 
title('Right CoP Norms'); 

%% save 
clearvars ans Contents i
save CoPNorms.mat LMeanCoP LStdCoP RMeanCoP RStdCoP 
    