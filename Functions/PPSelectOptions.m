function [SubjDemo, PPSettings, PPPlots] = PPSelectOptions

% % -------------------------------------------------------------------------
% Selects Options from Selector_File_Template
% -------------------------------------------------------------------------
% Author: Ricky Pimentel
% Date: December 21, 2016
%
% Description:	This function will read in the options from the
% selector_file_template and define the user inputs for kinemataics and
% kinetics data processing and plotting.
%
% Input:      Selector_File_Template.xlsx 
%
% Outputs: (all structural data)
%           KinematicsOptions           Options for Kinematics Processing

% % DEBUG LOOP %%%%%%%%%
% close all
% clear
% clc
% %%% END DEBUG LOOP %%%%%

%% Load Selector_File_Template
[~,~,RAW] = xlsread('PP_Selector_File.xlsx','Options');

%% Find Regions for each output
% Define start of regions
Ind = strcmp(RAW(:,1),'Subject Inputs') == 1;
SubjStart = find(Ind);

Ind = strcmp(RAW(:,1),'Processing Inputs') == 1;
ProcessingStart = find(Ind);

% Define end of regions
SubjEnd = ProcessingStart - 1;
Ind = strcmp(RAW(:,1),'END') == 1;
ProcessingEnd = find(Ind) - 1; 

%% Define Subject Inputs
for i = 1: SubjEnd - SubjStart
    SubjDemo(i).Input = RAW(SubjStart + i,1);
    if i == 1
        SubjDemo(i).Choice = RAW{SubjStart + i,2};
    else
        SubjDemo(i).Choice = double(cell2mat(RAW(SubjStart + i,2)));
    end
end

%% Define PP Option Range
for i = 1: ProcessingEnd - ProcessingStart
    PPOptions(i).Input = RAW(ProcessingStart + i,1);
    PPOptions(i).Choice = RAW(ProcessingStart + i,2);
end

%% Save PP Settings in Structure
for i = 1:length(PPOptions)
    if strcmp(PPOptions(i).Input, 'Initial Check')
        PPSettings.InitialCheck = PPOptions(i).Choice{1};
    end
    if strcmp(PPOptions(i).Input, 'PP Mat Type')
        PPSettings.PPMatType = PPOptions(i).Choice{1};
    end
    if strcmp(PPOptions(i).Input, 'Load Type')
        PPSettings.LoadNew = PPOptions(i).Choice{1};
    end
    if strcmp(PPOptions(i).Input, 'Auto Select All')
        PPSettings.AutoSelectAll = PPOptions(i).Choice{1};
    end
    if strcmp(PPOptions(i).Input, 'Mask Choice')
        PPSettings.MaskChoice = PPOptions(i).Choice{1};
    end
%     if strcmp(PPOptions(i).Input, 'Static Trial Analysis')
%         PPSettings.StaticAnalysis = PPOptions(i).Choice{1};
%     end
    if strcmp(PPOptions(i).Input, 'C3D Input')
        PPSettings.C3DInput = PPOptions(i).Choice{1};
    end
    if strcmp(PPOptions(i).Input, 'Movie Output')
        PPSettings.MovieOutput = PPOptions(i).Choice{1};
    end
    if strcmp(PPOptions(i).Input, 'Movie Types')
        PPSettings.MovieTypes = PPOptions(i).Choice{1};
    end
%     if strcmp(PPOptions(i).Input, 'GIF Speed')
%         PPSettings.GIFSpeed = PPOptions(i).Choice{1};
%     end
%     if strcmp(PPOptions(i).Input, 'Toe Walker')
%         PPSettings.ToeWalker = PPOptions(i).Choice{1};
%     end
    if strcmp(PPOptions(i).Input, 'Export Report')
        PPSettings.ExportReport = PPOptions(i).Choice{1};
    end
%     if strcmp(PPOptions(i).Input, 'Add in Offset')
%         PPSettings.AddOffset = PPOptions(i).Choice{1};
%     end

%% Plotting options
    if strcmp(PPOptions(i).Input, 'Plot Temporal Spatial')
        PPPlots.TempSpat = PPOptions(i).Choice{1};
    end
    if strcmp(PPOptions(i).Input, 'Plot Close Up')
        PPPlots.CloseUp= PPOptions(i).Choice{1};
    end
%     if strcmp(PPOptions(i).Input, 'COP Norms')
%         PPPlots.COPNorms = PPOptions(i).Choice{1};
%     end
    if strcmp(PPOptions(i).Input, 'Plot All Forces')
        PPPlots.AllForces = PPOptions(i).Choice{1};
    end
    if strcmp(PPOptions(i).Input, 'Plot All Pressures')
        PPPlots.AllPressures = PPOptions(i).Choice{1};
    end
    if strcmp(PPOptions(i).Input, 'Plot All Areas')
        PPPlots.AllAreas = PPOptions(i).Choice{1};
    end
    if strcmp(PPOptions(i).Input, 'Plot Average Forces')
        PPPlots.AvgForces = PPOptions(i).Choice{1};
    end
    if strcmp(PPOptions(i).Input, 'Plot Average Pressures')
        PPPlots.AvgPressures = PPOptions(i).Choice{1};
    end
    if strcmp(PPOptions(i).Input, 'Plot Validation') 
        PPPlots.Validation = PPOptions(i).Choice{1};
    end
end

%% Add defaults 
% to change defaults, must 
if isfield(PPSettings, 'InitialCheck') == 0
    PPSettings.InitialCheck = 'Yes';
end
if isfield(PPSettings, 'PPMatType') == 0
    PPSettings.PPMatType = 'Novel';
end
if isfield(PPSettings, 'LoadNew') == 0
    PPSettings.LoadNew = 'Raw';
end
if isfield(PPSettings, 'AutoSelectAll') == 0
    PPSettings.AutoSelectAll = 'Yes';
end
if isfield(PPSettings, 'MaskChoice') == 0
    PPSettings.MaskChoice = 'General';
end
if isfield(PPSettings, 'C3DInput') == 0
    PPSettings.C3DInput = 'No';
end
if isfield(PPSettings, 'MovieOutput') == 0
    PPSettings.MovieOutput = 'No';
end
if isfield(PPSettings, 'MovieTypes') == 0
    PPSettings.MovieTypes = 'Both';
end
% if isfield(PPSettings, 'GIFSpeed') == 0
%     PPSettings.GIFSpeed = 20;
% end
% if isfield(PPSettings, 'ToeWalker') == 0
%     PPSettings.ToeWalker = 'No';
% end
if isfield(PPSettings, 'ExportReport') == 0
    PPSettings.ExportReport = 'Yes';
end

% Plotting options
if isfield(PPSettings, 'TempSpat') == 0
    PPPlots.TempSpat =  'Yes';
end
if isfield(PPSettings, 'CloseUp') == 0
    PPPlots.CloseUp= 'Yes';
end
% if isfield(PPSettings, 'COP Norms') == 0
%     PPPlots.COPNorms = PPOptions(i).Choice{1};
% end
if isfield(PPPlots, 'AllForces') == 0
    PPPlots.AllForces = 'Yes';
end
if isfield(PPPlots, 'AllPressures') == 0
    PPPlots.AllPressures =  'No';
end
if isfield(PPPlots, 'AllAreas') == 0
    PPPlots.AllAreas =  'No';
end
if isfield(PPPlots, 'AvgForces') == 0
    PPPlots.AvgForces =  'No';
end
if isfield(PPPlots, 'AvgPressures') == 0
    PPPlots.AvgPressures =  'No';
end
if isfield(PPPlots, 'Validation') == 0
    PPPlots.AvgPressures =  'No';
end
end


