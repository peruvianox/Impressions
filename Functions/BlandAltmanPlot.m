function [CompAB, CompCD]  = BlandAltmanPlot(DataA, DataB, DataC, DataD, Title)

Data_A.Load = DataA; 
Data_A.Clean = Data_A.Load; 
Data_A.Nan = isnan(Data_A.Load); 
 
Data_B.Load = DataB; 
Data_B.Clean = Data_B.Load; 
Data_B.Nan = isnan(Data_B.Load); 

Data.Clean = or(Data_A.Nan, Data_B.Nan);  

Data_A.Clean(Data.Clean) = []; 
Data_B.Clean(Data.Clean) = []; 

if nargin == 4 || DataD ~= 0 % if more than 2 inputs
    Data_C.Load = DataC;
    Data_C.Clean = Data_C.Load;
    Data_C.Nan = isnan(Data_C.Load);
    
    Data_D.Load = DataD;
    Data_D.Clean = Data_D.Load;
    Data_D.Nan = isnan(Data_D.Load);
    
    Data2.Clean = or(Data_C.Nan, Data_D.Nan);
    
    Data_C.Clean(Data2.Clean) = [];
    Data_D.Clean(Data2.Clean) = [];
end

%% Compute Means, Diffs, and CIs
Data.Mean = mean([Data_A.Clean, Data_B.Clean],2); 
Data.Diff = Data_A.Clean - Data_B.Clean;
Data.Std = std(Data.Mean); 
Data.StdDiff = std(Data.Diff); 
Data.CI = Data.StdDiff .* 1.96; 
if nargin == 4   || DataD ~= 0
    Data2.Mean = mean([Data_C.Clean, Data_D.Clean],2);
    Data2.Diff = Data_C.Clean - Data_D.Clean;
    Data2.Std = std(Data2.Mean);
    Data2.StdDiff = std(Data2.Diff);
    Data2.CI = Data2.StdDiff .* 1.96;
end


%% plot bland altman 
BAplot = figure; 
plot(Data.Mean,Data.Diff,'b.'); 
hold on;
if nargin == 4 % if more than 2 inputs
    plot(Data2.Mean,Data2.Diff,'r.'); 
end

xlabel('Average Measure');
ylabel('Difference Between Measures'); 
if exist('Title') ~=1
    Title = inputdlg('What would you like to title the plot?');
end
title(Title); 
% plot trendlines and CIs
[Data.p,Data.s] = polyfit(Data.Mean, Data.Diff,1); 
xLow = min(Data.Mean); 
xHigh = max(Data.Mean); 
Data.x = linspace(xLow,xHigh,100);
[Data.y,Data.delta] = polyval(Data.p,Data.x,Data.s);
plot(Data.x,Data.y,'b-.', 'LineWidth',2);
% TrendErr = questdlg('Would you like to plot trendline errors?');
% if strcmp(TrendErr, 'Yes') == 1
%     plot(Data.x,Data.y+Data.delta,'b.');
%     plot(Data.x,Data.y-Data.delta,'b.');
% end
Data.MeanDiff = mean(Data.Diff); 
Data.MeanDiff_x = ones(length(Data.x))*Data.MeanDiff;
plot(Data.x,Data.MeanDiff_x,'b-', 'LineWidth',2); 
plot(Data.x,Data.MeanDiff_x + Data.CI, 'b:'); 
plot(Data.x,Data.MeanDiff_x - Data.CI, 'b:'); 
% text(mean([xLow,xHigh]), Data.MeanDiff_x(1) + Data.CI(1) - 1, strcat('y = ',num2str(Data.p(1)), 'x + ',num2str(Data.p(2)))); 

if nargin == 4 % if more than 2 inputs
    % plot trendlines and CIs
    [Data2.p,Data2.s] = polyfit(Data2.Mean, Data2.Diff,1);
    xLow = min(Data2.Mean);
    xHigh = max(Data2.Mean);
    Data2.x = linspace(xLow,xHigh,100);
    [Data2.y,Data2.delta] = polyval(Data2.p,Data2.x,Data2.s);
    plot(Data2.x,Data2.y,'r-.', 'LineWidth',2);
%     if strcmp(TrendErr, 'Yes') == 1
%         plot(Data2.x,Data2.y+Data2.delta,'r.');
%         plot(Data2.x,Data2.y-Data2.delta,'r.');
%     end
    Data2.MeanDiff = mean(Data2.Diff);
    Data2.MeanDiff_x = ones(length(Data2.x))*Data2.MeanDiff;
    plot(Data2.x,Data2.MeanDiff_x,'r-', 'LineWidth',2);
    plot(Data2.x,Data2.MeanDiff_x + Data2.CI, 'r:');
    plot(Data2.x,Data2.MeanDiff_x - Data2.CI, 'r:');
end

CompAB = Data;

if nargin == 4
    CompCD = Data2;
else
    CompCD = [];
end

end
