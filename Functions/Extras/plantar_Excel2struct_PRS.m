function[AnalyzePPImages] =  plantar_Excel2struct_PRS
% Paul Silva
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
    % Pre-allocate all big matrices
    ndata = sparse(100000,63);
    TXT = sparse(100000,1);
%%
%%%---------------------------------------------------------------------%%%
    % Import plantar pressure data from Excel spreadsheets. Must be saved with
    % the Excel binary extension. Name of folder with data is used as
    % identifier for subject when output variables are exported to Excel
    % spreadsheet
    uiwait(msgbox('Please select .xlsb of entire plate roll-off',...
        'Plate Roll-Off'));
    [filename, pathname] = uigetfile('*.xlsb');
    if isequal(filename,0)
       msgbox('User pressed cancel')
       return
    end
    [ndata, TXT, ~] = xlsread([pathname filename], 1);
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
   
    last_frame = cell2mat(TXT(length(TXT)));
    pstart = strfind(last_frame,'e ');
    pstart = pstart + 2;
    pend = strfind(last_frame,' (');
    pend = pend - 1;
    tot_frames = last_frame(pstart:pend);
    tot_frames = str2num(tot_frames);
    pstart = strfind(last_frame,'(');
    pend = strfind(last_frame,' ms');
    tot_time = last_frame(pstart+1:pend-1);
    tot_time = str2num(tot_time);
    end_frame = halt(1)-1;
    
    
% % % %     last_frame = msgbox(TXT(length(TXT)))
% % % %     tot_frames = inputdlg('Please enter number of frames based on previous output: ');
% % % %     tot_frames = str2num(tot_frames{1});
% % % %     tot_time = inputdlg('Please enter total time (ms) based on previous output: ');
% % % %     tot_time = str2num(tot_time{1});
% % % %     end_frame = halt(1)-1;
    
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
% For now, just save the individual Excel spreadsheet as it's own .mat file

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
%%        
%%%---------------------------------------------------------------------%%%        
        
% Import center of pressure data from Excel spreadsheets. Must be saved with
% the Excel binary extension. Name of folder helps to determine if center
% of pressure matches the plantar pressure imported from above
uiwait(msgbox( 'Select .xlsb of entire centre of force',...
        'Center of Force'));
    if isequal(filename,0)
       msgbox('User pressed cancel')
       return
    end
[file_name, path_name] = uigetfile('*.xlsb');
[n_data, TXT, ~] = xlsread([path_name file_name], 1);
check_study_date = file_name(1:(end-30)); 

while study_date ~= check_study_date
    [file_name, path_name] = uigetfile('*.xlsb', 'Select a different .xlsb of entire centre of force, ensuring it matches the entire plate roll-off.');
    [n_data, TXT, ~] = xlsread([path_name file_name], 1);
    check_study_date = file_name(1:(end-30)); 
end

AnalyzePPImages.center_of_pressure = n_data; 

save(mat_filename, 'AnalyzePPImages'); 
% end

