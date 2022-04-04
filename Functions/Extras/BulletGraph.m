function [] = BulletGraph(Type, Poor, Ok, Actual, Target, Limits, Color, Num, SubPlot, Pair)

% BulletGraph
% creates a tiered set of bar plots that displays how a measure
% is performing relative to a target. For more information on Bullet Graphs, see:
% http://www.perceptualedge.com/articles/misc/Bullet_Graph_Design_Spec.pdf
%
% Definitions
% Type - can be 'V' for Vertical, 'H' for Hoizontal, 'ReverseV' for Reverse
% Vertical, or 'ReverseH' for Reverse Horizontal.
% Poor - the divide between Poor and Ok levels
% Ok - the divide between Ok and Good levels
% Actual - the measured or current value
% Target - the desired or optimal value
% Limits - the Limits of the graph
% Color - color of Actual value in ColorSpec, defaults to 'k'
% Num - the number of bars desired on the plot
% SubPlot - whether this graph is destined to be as a part of a subplot. 

% Ricky Pimentel
% Center for Gait and Movement Analysis, Children's Hospital Colorado,
% Aurora, CO, USA
% August 2017

if strcmp(SubPlot,'Yes') ~= 1
    figure('Position',[100 100 200 500]);
end


%% If vertical graph desired
if strcmp(Type,'V') == 1
    
    hold on;
    
    if exist('Color') ~= 1 % if Color isn't specified, default to black
        Color = [0 0 0];
    end
    
    % bar background colors
    bar(1,Ok,2,'FaceColor','k', 'EdgeColor', [0 0 0],'LineStyle','none');
    alpha(0.8);
    bar(1,Poor,2,'FaceColor','k', 'EdgeColor', [0 0 0],'LineStyle','none');
    alpha(0.4);
    
    if Num == 1 % if a single measure is desired
        bar(1,Actual,0.25,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
    elseif Num == 2  % for 2 measures
        if strcmp(Pair, 'Yes') == 1  
        bar(0.9,Actual(1),0.2,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        bar(1.1,Actual(2),0.2,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        else
        bar(0.85,Actual(1),0.2,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        bar(1.15,Actual(2),0.2,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        end
    elseif Num == 3  % for 3 measures
        bar(0.75,Actual(1),0.2,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        bar(1,Actual(2),0.2,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        bar(1.25,Actual(3),0.2,'FaceColor',Color(3,:), 'EdgeColor', [0 0 0],'LineStyle','none');
    elseif Num == 4 % for 4 measures
         if strcmp(Pair, 'Yes') == 1
             bar(0.75,Actual(1),0.15,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
             bar(0.9,Actual(2),0.15,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
             bar(1.1,Actual(3),0.15,'FaceColor',Color(3,:), 'EdgeColor', [0 0 0],'LineStyle','none');
             bar(1.25,Actual(4),0.15,'FaceColor',Color(4,:), 'EdgeColor', [0 0 0],'LineStyle','none');
         else             
             bar(0.7,Actual(1),0.15,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
             bar(0.9,Actual(2),0.15,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
             bar(1.1,Actual(3),0.15,'FaceColor',Color(3,:), 'EdgeColor', [0 0 0],'LineStyle','none');
             bar(1.3,Actual(4),0.15,'FaceColor',Color(4,:), 'EdgeColor', [0 0 0],'LineStyle','none');
         end
    elseif Num == 6 % for 6 measures
        if strcmp(Pair, 'Yes') == 1
            bar(0.6,Actual(1),0.125,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(0.725,Actual(2),0.125,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(0.9375,Actual(3),0.125,'FaceColor',Color(3,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(1.0625,Actual(4),0.125,'FaceColor',Color(4,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(1.275,Actual(5),0.125,'FaceColor',Color(5,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(1.4,Actual(6),0.125,'FaceColor',Color(6,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        else
            bar(0.6875,Actual(1),0.125,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(0.8125,Actual(2),0.125,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(0.9375,Actual(3),0.125,'FaceColor',Color(3,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(1.0625,Actual(4),0.125,'FaceColor',Color(4,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(1.1875,Actual(5),0.125,'FaceColor',Color(5,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(1.3125,Actual(6),0.125,'FaceColor',Color(6,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        end
    elseif Num == 8 % for 8 measures
        if strcmp(Pair, 'Yes') == 1
            bar(0.49,Actual(1),0.12,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(0.61,Actual(2),0.12,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(0.79,Actual(3),0.12,'FaceColor',Color(3,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(0.91,Actual(4),0.12,'FaceColor',Color(4,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(1.09,Actual(5),0.12,'FaceColor',Color(5,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(1.21,Actual(6),0.12,'FaceColor',Color(6,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(1.39,Actual(7),0.12,'FaceColor',Color(7,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(1.51,Actual(8),0.12,'FaceColor',Color(8,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        else
            bar(0.5625,Actual(1),0.125,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(0.6875,Actual(2),0.125,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(0.8125,Actual(3),0.125,'FaceColor',Color(3,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(0.9375,Actual(4),0.125,'FaceColor',Color(4,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(1.0625,Actual(5),0.125,'FaceColor',Color(5,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(1.1875,Actual(6),0.125,'FaceColor',Color(6,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(1.3125,Actual(7),0.125,'FaceColor',Color(7,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(1.4375,Actual(8),0.125,'FaceColor',Color(8,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        end
    end
    
    % set target value and limits
    h = hline(Target, 'k-');
    set(h,'LineWidth',2);
    xlim([0.4 1.6]);
    ylim(Limits);
    
    % Set Axes Properties
    ax = gca;
    ax.XTick = [];
    
    %% If horizontal graph desired
elseif strcmp(Type,'H') == 1
    %  figure('Position',[100 100 500 200]);
    hold on;
    
    if exist('Color') ~= 1 % if Color isn't specified, default to black
        Color = 'k';
    end
    
    % bar background colors
    barh(1,Ok,2,'FaceColor','k', 'EdgeColor', [0 0 0],'LineStyle','none');
    alpha(0.8);
    barh(1,Poor,2,'FaceColor','k', 'EdgeColor', [0 0 0],'LineStyle','none');
    alpha(0.4);
    
    if Num == 1 % if a single measure is desired
        barh(1,Actual,0.25,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
    elseif Num == 2  % for 2 measures
        if strcmp(Pair, 'Yes') == 1  
        barh(0.9,Actual(1),0.2,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        barh(1.1,Actual(2),0.2,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        else
        barh(0.85,Actual(1),0.2,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        barh(1.15,Actual(2),0.2,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        end
    elseif Num == 3  % for 3 measures
        barh(0.75,Actual(1),0.2,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        barh(1,Actual(2),0.2,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        barh(1.25,Actual(3),0.2,'FaceColor',Color(3,:), 'EdgeColor', [0 0 0],'LineStyle','none');
    elseif Num == 4 % for 4 measures
         if strcmp(Pair, 'Yes') == 1
             barh(0.75,Actual(1),0.15,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
             barh(0.9,Actual(2),0.15,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
             barh(1.1,Actual(3),0.15,'FaceColor',Color(3,:), 'EdgeColor', [0 0 0],'LineStyle','none');
             barh(1.25,Actual(4),0.15,'FaceColor',Color(4,:), 'EdgeColor', [0 0 0],'LineStyle','none');
         else             
             barh(0.7,Actual(1),0.15,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
             barh(0.9,Actual(2),0.15,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
             barh(1.1,Actual(3),0.15,'FaceColor',Color(3,:), 'EdgeColor', [0 0 0],'LineStyle','none');
             barh(1.3,Actual(4),0.15,'FaceColor',Color(4,:), 'EdgeColor', [0 0 0],'LineStyle','none');
         end
    elseif Num == 6 % for 6 measures
        if strcmp(Pair, 'Yes') == 1
            barh(0.6,Actual(1),0.125,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(0.725,Actual(2),0.125,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(0.9375,Actual(3),0.125,'FaceColor',Color(3,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(1.0625,Actual(4),0.125,'FaceColor',Color(4,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(1.275,Actual(5),0.125,'FaceColor',Color(5,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(1.4,Actual(6),0.125,'FaceColor',Color(6,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        else
            barh(0.6875,Actual(1),0.125,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(0.8125,Actual(2),0.125,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(0.9375,Actual(3),0.125,'FaceColor',Color(3,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(1.0625,Actual(4),0.125,'FaceColor',Color(4,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(1.1875,Actual(5),0.125,'FaceColor',Color(5,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(1.3125,Actual(6),0.125,'FaceColor',Color(6,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        end
    elseif Num == 8 % for 8 measures
        if strcmp(Pair, 'Yes') == 1
            barh(0.49,Actual(1),0.12,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(0.61,Actual(2),0.12,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(0.79,Actual(3),0.12,'FaceColor',Color(3,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(0.91,Actual(4),0.12,'FaceColor',Color(4,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(1.09,Actual(5),0.12,'FaceColor',Color(5,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(1.21,Actual(6),0.12,'FaceColor',Color(6,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(1.39,Actual(7),0.12,'FaceColor',Color(7,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(1.51,Actual(8),0.12,'FaceColor',Color(8,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        else
            barh(0.5625,Actual(1),0.125,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(0.6875,Actual(2),0.125,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(0.8125,Actual(3),0.125,'FaceColor',Color(3,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(0.9375,Actual(4),0.125,'FaceColor',Color(4,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(1.0625,Actual(5),0.125,'FaceColor',Color(5,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(1.1875,Actual(6),0.125,'FaceColor',Color(6,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(1.3125,Actual(7),0.125,'FaceColor',Color(7,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(1.4375,Actual(8),0.125,'FaceColor',Color(8,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        end
    end
    %barh(1,Actual,0.25,'FaceColor',Color, 'EdgeColor', [0 0 0],'LineStyle','none');
    
    h = vline(Target, 'k-');
    set(h,'LineWidth',2);
    
    xlim(Limits);
    ylim([0.4 1.6]);
    
    % Set Axes Properties
    ax = gca;
    ax.YTick = [];
    
    %% if Reverse Vertical graph desired
elseif strcmp(Type,'ReverseV') == 1
    
    % figure('Position',[100 100 200 500]);
    hold on;
    
    if exist('Color') ~= 1 % if Color isn't specified, default to black
        Color = 'k';
    end
    
    % bar background colors
    bar(1,Limits(2),2,'FaceColor','k', 'EdgeColor', [0 0 0],'LineStyle','none');
    bar(1,Poor,2,'FaceColor','w', 'EdgeColor', [0 0 0],'LineStyle','none');
    alpha(0.5);
    bar(1,Ok,2,'FaceColor','w', 'EdgeColor', [0 0 0],'LineStyle','none');
    
     if Num == 1 % if a single measure is desired
        bar(1,Actual,0.25,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
    elseif Num == 2  % for 2 measures
        if strcmp(Pair, 'Yes') == 1  
        bar(0.9,Actual(1),0.2,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        bar(1.1,Actual(2),0.2,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        else
        bar(0.85,Actual(1),0.2,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        bar(1.15,Actual(2),0.2,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        end
    elseif Num == 3  % for 3 measures
        bar(0.75,Actual(1),0.2,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        bar(1,Actual(2),0.2,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        bar(1.25,Actual(3),0.2,'FaceColor',Color(3,:), 'EdgeColor', [0 0 0],'LineStyle','none');
    elseif Num == 4 % for 4 measures
         if strcmp(Pair, 'Yes') == 1
             bar(0.75,Actual(1),0.15,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
             bar(0.9,Actual(2),0.15,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
             bar(1.1,Actual(3),0.15,'FaceColor',Color(3,:), 'EdgeColor', [0 0 0],'LineStyle','none');
             bar(1.25,Actual(4),0.15,'FaceColor',Color(4,:), 'EdgeColor', [0 0 0],'LineStyle','none');
         else             
             bar(0.7,Actual(1),0.15,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
             bar(0.9,Actual(2),0.15,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
             bar(1.1,Actual(3),0.15,'FaceColor',Color(3,:), 'EdgeColor', [0 0 0],'LineStyle','none');
             bar(1.3,Actual(4),0.15,'FaceColor',Color(4,:), 'EdgeColor', [0 0 0],'LineStyle','none');
         end
    elseif Num == 6 % for 6 measures
        if strcmp(Pair, 'Yes') == 1
            bar(0.6,Actual(1),0.125,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(0.725,Actual(2),0.125,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(0.9375,Actual(3),0.125,'FaceColor',Color(3,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(1.0625,Actual(4),0.125,'FaceColor',Color(4,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(1.275,Actual(5),0.125,'FaceColor',Color(5,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(1.4,Actual(6),0.125,'FaceColor',Color(6,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        else
            bar(0.6875,Actual(1),0.125,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(0.8125,Actual(2),0.125,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(0.9375,Actual(3),0.125,'FaceColor',Color(3,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(1.0625,Actual(4),0.125,'FaceColor',Color(4,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(1.1875,Actual(5),0.125,'FaceColor',Color(5,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(1.3125,Actual(6),0.125,'FaceColor',Color(6,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        end
    elseif Num == 8 % for 8 measures
        if strcmp(Pair, 'Yes') == 1
            bar(0.49,Actual(1),0.12,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(0.61,Actual(2),0.12,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(0.79,Actual(3),0.12,'FaceColor',Color(3,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(0.91,Actual(4),0.12,'FaceColor',Color(4,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(1.09,Actual(5),0.12,'FaceColor',Color(5,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(1.21,Actual(6),0.12,'FaceColor',Color(6,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(1.39,Actual(7),0.12,'FaceColor',Color(7,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(1.51,Actual(8),0.12,'FaceColor',Color(8,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        else
            bar(0.5625,Actual(1),0.125,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(0.6875,Actual(2),0.125,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(0.8125,Actual(3),0.125,'FaceColor',Color(3,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(0.9375,Actual(4),0.125,'FaceColor',Color(4,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(1.0625,Actual(5),0.125,'FaceColor',Color(5,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(1.1875,Actual(6),0.125,'FaceColor',Color(6,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(1.3125,Actual(7),0.125,'FaceColor',Color(7,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            bar(1.4375,Actual(8),0.125,'FaceColor',Color(8,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        end
    end
    % bar(1,Actual,0.25,'FaceColor',Color, 'EdgeColor', [0 0 0],'LineStyle','none');
    
    h = hline(Target, 'k-');
    set(h,'LineWidth',2);
    
    ylim(Limits);
    xlim([0.4 1.6]);
    
    % Set Axes Properties
    ax = gca;
    ax.XTick = [];
    
    %% if Reverse Horizontal graph desired
elseif strcmp(Type,'ReverseH') == 1
    
    % figure('Position',[100 100 500 200]);
    hold on;
    
    if exist('Color') ~= 1 % if Color isn't specified, default to black
        Color = 'k';
    end
    
    % bar background colors
    barh(1,Limits(2),2,'FaceColor','k', 'EdgeColor', [0 0 0],'LineStyle','none');
    %  alpha(0.4);
    barh(1,Poor,2,'FaceColor','w', 'EdgeColor', [0 0 0],'LineStyle','none');
    % barh(1,Poor,1,'FaceColor','k', 'EdgeColor', [0 0 0],'LineStyle','none');
    alpha(0.5);
    barh(1,Ok,2,'FaceColor','w', 'EdgeColor', [0 0 0],'LineStyle','none');
    
        if Num == 1 % if a single measure is desired
        barh(1,Actual,0.25,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
    elseif Num == 2  % for 2 measures
        if strcmp(Pair, 'Yes') == 1  
        barh(0.9,Actual(1),0.2,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        barh(1.1,Actual(2),0.2,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        else
        barh(0.85,Actual(1),0.2,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        barh(1.15,Actual(2),0.2,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        end
    elseif Num == 3  % for 3 measures
        barh(0.75,Actual(1),0.2,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        barh(1,Actual(2),0.2,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        barh(1.25,Actual(3),0.2,'FaceColor',Color(3,:), 'EdgeColor', [0 0 0],'LineStyle','none');
    elseif Num == 4 % for 4 measures
         if strcmp(Pair, 'Yes') == 1
             barh(0.75,Actual(1),0.15,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
             barh(0.9,Actual(2),0.15,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
             barh(1.1,Actual(3),0.15,'FaceColor',Color(3,:), 'EdgeColor', [0 0 0],'LineStyle','none');
             barh(1.25,Actual(4),0.15,'FaceColor',Color(4,:), 'EdgeColor', [0 0 0],'LineStyle','none');
         else             
             barh(0.7,Actual(1),0.15,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
             barh(0.9,Actual(2),0.15,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
             barh(1.1,Actual(3),0.15,'FaceColor',Color(3,:), 'EdgeColor', [0 0 0],'LineStyle','none');
             barh(1.3,Actual(4),0.15,'FaceColor',Color(4,:), 'EdgeColor', [0 0 0],'LineStyle','none');
         end
    elseif Num == 6 % for 6 measures
        if strcmp(Pair, 'Yes') == 1
            barh(0.6,Actual(1),0.125,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(0.725,Actual(2),0.125,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(0.9375,Actual(3),0.125,'FaceColor',Color(3,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(1.0625,Actual(4),0.125,'FaceColor',Color(4,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(1.275,Actual(5),0.125,'FaceColor',Color(5,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(1.4,Actual(6),0.125,'FaceColor',Color(6,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        else
            barh(0.6875,Actual(1),0.125,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(0.8125,Actual(2),0.125,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(0.9375,Actual(3),0.125,'FaceColor',Color(3,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(1.0625,Actual(4),0.125,'FaceColor',Color(4,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(1.1875,Actual(5),0.125,'FaceColor',Color(5,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(1.3125,Actual(6),0.125,'FaceColor',Color(6,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        end
    elseif Num == 8 % for 8 measures
        if strcmp(Pair, 'Yes') == 1
            barh(0.49,Actual(1),0.12,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(0.61,Actual(2),0.12,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(0.79,Actual(3),0.12,'FaceColor',Color(3,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(0.91,Actual(4),0.12,'FaceColor',Color(4,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(1.09,Actual(5),0.12,'FaceColor',Color(5,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(1.21,Actual(6),0.12,'FaceColor',Color(6,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(1.39,Actual(7),0.12,'FaceColor',Color(7,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(1.51,Actual(8),0.12,'FaceColor',Color(8,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        else
            barh(0.5625,Actual(1),0.125,'FaceColor',Color(1,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(0.6875,Actual(2),0.125,'FaceColor',Color(2,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(0.8125,Actual(3),0.125,'FaceColor',Color(3,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(0.9375,Actual(4),0.125,'FaceColor',Color(4,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(1.0625,Actual(5),0.125,'FaceColor',Color(5,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(1.1875,Actual(6),0.125,'FaceColor',Color(6,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(1.3125,Actual(7),0.125,'FaceColor',Color(7,:), 'EdgeColor', [0 0 0],'LineStyle','none');
            barh(1.4375,Actual(8),0.125,'FaceColor',Color(8,:), 'EdgeColor', [0 0 0],'LineStyle','none');
        end
    end
   % barh(1,Actual,0.25,'FaceColor',Color, 'EdgeColor', [0 0 0],'LineStyle','none');
    
    h = vline(Target, 'k-');
    set(h,'LineWidth',2);
    
    xlim(Limits);
    ylim([0.4 1.6]);
    
    % Set Axes Properties
    ax = gca;
    ax.YTick = [];
    
end

