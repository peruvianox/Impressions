function [MAT] = plantar_Excel2struct (pathname, filenameRO, filenameCOF)
% Nicole Look
% Start: July 24, 2014
% Update: July 31, 2014
%%%---------------------------------------------------------------------%%%
% Directions for using plantar_Excel2struct
%
% Export "Entire Plate Roll-Off" and "Entire Centre of Force" data from RSscan.
% Save files in Excel using the binary extension .xlsb.
%
% Inputs:
%   1. Excel spreadsheet of "Entire Plate Roll-Off" as .xlsb
%   2. Excel spreadsheet of "Entire Centre of Force" as .xlsb
%
% Outputs:
%   1. Condenses Excel spreadsheet into a .mat file that can be used 
%   by mushy_foot.m
%
%%%---------------------------------------------------------------------%%%

%%
%%%---------------------------------------------------------------------%%%
    %%DEBUG LOOP%%
%     close all
%     clear all 
%     clc
    %%END DEBUG LOOP%%
    
    % Pre-allocate all big matrices
%     ndata = zeros(100000,63);
%     TXT = zeros(100000,1);
%     alldata = zeros(100000,63);

%%%---------------------------------------------------------------------%%%

%%
%%%---------------------------------------------------------------------%%%
    % Import plantar pressure data from Excel spreadsheets. Must be saved with
    % the Excel binary extension. Name of folder with data is used as
    % identifier for subject when output variables are exported to Excel
    % spreadsheet
 [filename, pathname] = uigetfile('*.xlsb', 'Select .xlsb of entire plate roll-off');
    [ndata, TXT, alldata] = xlsread([pathname filename], 1);
    locate_slash = find(pathname == '\');
    subject_name = pathname(locate_slash(end-1)+1:locate_slash(end)-1);
    study_date = filename(1:(end-29));

    % Create 3-dimensional matrix in which each frame in time contains the x- and y- 
    % coordinates. The z-coordinate is composed of each frame. In general,
    % there are usually 248 frames at 7.94 ms per frame. This is not always
    % true, hence the user input

    % Frames are separated by two rows of NaN
    not_number = isnan(ndata(:,1));
    not_number = double(not_number);
    halt = find(not_number == 1);
    
    % Display total number of frames for user 
    last_frame = TXT(length(TXT));
    tot_frames = input('Please enter number of frames based on previous output: ');
    tot_time = input('Please enter total time (ms) based on previous output: ');
    end_frame = halt(1)-1;
    
    time_matrix = zeros(end_frame,63,tot_frames);
    n_data = zeros(tot_frames,2);

    time_matrix(:,:,1) = ndata(1:end_frame,:);
    [c, r] = size(time_matrix);
    increment = c + 2;

    start = end_frame + 3;
    end_frame = end_frame + increment;

    for frame = 2:tot_frames
        time_matrix(:,:,frame) = ndata(start:end_frame,:);
        start = end_frame + 3;
        end_frame = end_frame+increment;
    end
%%%---------------------------------------------------------------------%%%


%%%---------------------------------------------------------------------%%%
% For now, just save the individual Excel spreadsheet as it's own .mat file

        AnalyzePPImages.alldata = alldata;
        AnalyzePPImages.time_matrix = time_matrix;
        AnalyzePPImages.subject_name = subject_name;
        AnalyzePPImages.study_date = study_date; 
        AnalyzePPImages.total_time = tot_time;
        
        prompt_msg = 'Please enter name of mat file: '; 
        prompt = {prompt_msg};
        dlg_title = 'Export Data';
        num_lines = 1;
        def = {''};
        mat_filename = strcat(char(cell2mat( inputdlg(prompt,dlg_title,num_lines,def) )), '.mat');
        save(mat_filename, 'AnalyzePPImages');
        
        
%%%---------------------------------------------------------------------%%%        
 

%%        
%%%---------------------------------------------------------------------%%%        
        
% Import center of pressure data from Excel spreadsheets. Must be saved with
% the Excel binary extension. Name of folder helps to determine if center
% of pressure matches the plantar pressure imported from above
[file_name, path_name] = uigetfile('*.xlsb', 'Select .xlsb of entire centre of force');
[n_data, TXT, alldata] = xlsread([pathname filename], 1);
check_study_date = file_name(1:(end-30)); 

while study_date ~= check_study_date
    [file_name, path_name] = uigetfile('*.xlsb', 'Select a different .xlsb of entire centre of force, ensuring it matches the entire plate roll-off.');
    [n_data, TXT, alldata] = xlsread([path_name file_name], 1);
    check_study_date = file_name(1:(end-30)); 
end

AnalyzePPImages.center_of_pressure = n_data; 

MAT = AnalyzePPImages.center_of_pressure;

% save(mat_filename, 'AnalyzePPImages');

%%%---------------------------------------------------------------------%%% 

end

