function[AnalyzePPImages, TrialType, SaveName] =  plantar_load(FileNames, TrialType, SaveName)
% Authors: Ricky Pimentel and Paul Silva
% Center for Gait and Movement Analysis, Children's Hospital Colorado
%%%---------------------------------------------------------------------%%%
% plantar_load
%
% Export "Entire Plate Roll-Off" and "Entire Centre of Force" data from
% RSscan for each plantar pressure trial to be analyzed.
%
% Inputs:
%   1. Generic File containing "Entire Plate Roll-Off"
%   2. Generic File containing  "Entire Centre of Force"
% Note: can load multiple pairs of these inputs - will all be saved in a
% single .mat file
%
% Outputs:
%   1. Condenses generic file into a .mat file that can be used
%   by Impressions for Plantar pressure Analysis AND returns a variable
%   containing the PP data from the incoming files.
%
%%%---------------------------------------------------------------------%%%

%% Select type of input data if not already defined
if exist('TrialType') == 0
    Question = 'What type of data are you importing?';
    TrialType= questdlg(Question,'Input Type','RSScan','Novel','Novel');
end
%% load series of data
if strcmp(TrialType,'Novel') == 1
    if exist('FileNames') == 0 || FileNames == 0
        filename = uigetfile('*.lst*','Select Novel data .lst file', 'MultiSelect','On');
    end
    if iscell(filename)
        NumTrials = length(filename); 
    else
        NumTrials = 1; 
    end
    
    for z = 1:NumTrials
        %% Initialize variables.
        delimiter = '\t';
        startRow = 9;
        
        % Read columns of data as text:
        % For more information, see the TEXTSCAN documentation.
        formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

       % Open the text file.
       if NumTrials == 1
           fileID = fopen(filename,'r');
       else
           fileID = fopen(char(filename(z)),'r');
       end
        
       % Read columns of data according to the format.
       % This call is based on the structure of the file used to generate this
       % code. If an error occurs for a different file, try regenerating the code
       % from the Import Tool.
       textscan(fileID, '%[^\n\r]', startRow-1 , 'WhiteSpace', '', 'ReturnOnError', false);
       dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'ReturnOnError', false);
       
       % Close the text file.
        fclose(fileID);
        
        dataArray{1,1}(1) = []; % delete first cell to ensure similar length of all input types
        m = length(dataArray{1,1}); % number of rows
        n = length(dataArray); % number of columms
        PRO = zeros(m,n); % initialize data to be transcribed
        
        ind = find(str2double(dataArray{1,1})==1); % find indicies of all the 1's - where the matricies start
        nFrames = length(ind); % number of frames in the dataset
        
        tic;
        for i = 1:5 % slower method for first few columns to contain the text fields
            PRO(:,i) = str2double(dataArray{1,i}); % input dataArray cells into one large matrix
        end
        for i = 6:89 % faster method for the rest to only include numbers
            PRO(:,i) = str2doubles(dataArray{1,i}); % input dataArray cells into one large matrix
        end
        toc
        
        row = 288; % define matrix size
        col = 88;
        
        TimeMatrix =zeros(row,col,nFrames); % initialize TimeMatrix of data
        
        for i = 1:nFrames
            TimeMatrix(:,:,i) = PRO(ind(i)+1:ind(i)+row,2:89); % split large matrix into a 3D series of row,col,frames
        end
        
%         TimeMatrix = ans.time_matrix;
%         nFrames = 199;
%         % check to make sure import is correct
%         uiwait(msgbox(['Video of plantar pressures.']));
%         figure('Position',[100 100 200 600]);
%         for i = 1:nFrames
%             contour(TimeMatrix(:,:,i), 20);
%             f(i)= getframe;
%         end
%         uiwait(msgbox(['Video of plantar pressures.']));
%         close;
%         figure('Position',[100 100 200 600]);
%         movie(f,10,50)
%         
        %% Initialize variables for importing trial metrics
        delimiter = '\t';
        startRow = 1;
        endRow = 8;
        formatSpec = '%*s%s%[^\n\r]';
        
        %% Open the text file.
        if NumTrials == 1
            fileID = fopen(filename,'r');
        else
            fileID = fopen(char(filename(z)),'r');
        end
        
        %% Read columns of data according to the format.
        textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false);
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'ReturnOnError', false);
        
        %% Close the text file.
        fclose(fileID);
        
        %% save trial parameters
        dataString = [dataArray{1:end-1}]; % find data description string
        C = strsplit(dataString{2}, '\'); % split the string
        
%         subject_name = C(7); % location of the subject name
%         if strcmp(TrialType, 'RSScan') == 1
%             study_date = datestr(datenum(char(C(6)), 'yyyymmdd')); % location of the date
%         else
%             study_date = datestr(datenum(char(C(6)), 'yyyymmdd')); % location of the date
%         end
        
        %% Save Data into Structure format
        if NumTrials == 1
            AnalyzePPImages(z).file_name = filename;
        else
            AnalyzePPImages(z).file_name = filename(z);
        end
        AnalyzePPImages(z).time_matrix = TimeMatrix;
%         AnalyzePPImages(z).subject_name = subject_name;
%          AnalyzePPImages(z).study_date = study_date;
        
        %% compute center of pressure
%         center_of_pressure = zeros(nFrames,2);
%         for i = 1:nFrames
%             CC = bwconncomp(TimeMatrix(:,:,i));
%             CCC = labelmatrix(CC);
%             CCC(CCC~=0) = 1;
%             Cent = regionprops(CCC,'centroid');
%             if isempty(Cent) == 0
%                 center_of_pressure(i,1:2) = Cent.Centroid;
%             else
%                 center_of_pressure(i,1:2) = NaN;
%             end
%         end
%         
% % %         contour(TimeMatrix(:,:,nFrames), 100); 
% %         figure('Position',[100,100, 400, 700]); 
% %         axis equal; 
% %         hold on; 
% %         for i = 1:nFrames
% %         contour(TimeMatrix(:,:,i),30); 
% %         hold on; 
% %         plot(center_of_pressure(i,1),center_of_pressure(i,2) , 'r*'); 
% %         pause(0.0001)
% %         hold off; 
% %         end
%         
%         AnalyzePPImages(z).center_of_pressure = center_of_pressure;
        
    end
    
else 
    %% if RS Scan data
    if exist('FileNames') == 0 || FileNames == 0
        % Determine how many trials to import
        %Question = 'How many trials would you like to load?';
        
        [FileNames, ~] = uigetfile('*.*','Select "Entire Plate Roll Off" file','MultiSelect','On');
        if iscell(FileNames)
            NumTrialsImport = length(FileNames);
        else
            NumTrialsImport = 1;
        end
    else
        if iscell(FileNames)
            NumTrialsImport = length(FileNames);
        else
            NumTrialsImport = 1;
        end
    end
    
    for z = 1:NumTrialsImport
%         %% Import Entire Center of Force File
%         [file_name, ~] = uigetfile('*.*','Select "Entire Center of Force" file');
%         
%         CoFFileID = fopen(file_name);
%         CoF = textscan(CoFFileID,'%s %s %*s %*s %f %f');
%         
%         center_of_pressure = [CoF{1,3} CoF{1,4}];
%         indcop = find(isnan(center_of_pressure(:,1)),50);
%         center_of_pressure(indcop,:) = [];
        
        %% Import Entire Plate Roll Off File
%         if exist('FileNames') == 0 % if no filenames input
%             [filename, pathname] = uigetfile('*.*','Select "Entire Plate Roll Off" file');
%         else 
            if NumTrialsImport == 1
                filename = FileNames;
                [pathname, ~, ~] = fileparts(filename);
            else
                filename = FileNames{z};
                 [pathname, ~, ~] = fileparts(filename);
            end
%         end
        % ROFileID = fopen(filename);
        delimiter = '\t';
        startRow = 2;
        
        % Read columns of data as strings:
        % For more information, see the TEXTSCAN documentation.
        formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';
        
        % Open the text file.
        ROFileID = fopen(filename,'r');
        
        % Read columns of data according to format string.
        % This call is based on the structure of the file used to generate this
        % code. If an error occurs for a different file, try regenerating the code
        % from the Import Tool.
        textscan(ROFileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false);
        dataArray = textscan(ROFileID, formatSpec, 'Delimiter', delimiter, 'ReturnOnError', false);
        
        % Close the text file.
        fclose(ROFileID);
        clearvars col fileID CoFFileID ROFileID startRow indcop formatSpec
        
        %% translate input data into a double-type matrix
        % this section takes the longest - approx 2 minutes!
        tic;
        m = length(dataArray{1,1});
        if m > 63992
            for i = 1:63
                dataArray{1,i}(63993:end) = [];
            end
            m = 63992;
        end
        n = length(dataArray);
        if n == 64
            dataArray{1,64} = [];
            n = 63;
        end
        PRO = zeros(m,n);
        
        for i = 1:5 % slower method for first few columns to contain the text fields
            PRO(:,i) = str2double(dataArray{1,i});
        end
        for i = 6:63 % faster method for the rest to only include numbers
            PRO(:,i) = str2doubles(dataArray{1,i}); % input dataArray cells into one large matrix
        end
        
        ind = find(isnan(PRO(:,1)), 248);
        toc
        
            
        %%  Define Subject Name and Study Date
%         locate_slash = find(pathname == '\');
%         subject_name = pathname(locate_slash(end-1)+1:locate_slash(end)-1);
%         study_date = filename(1:(end-29));
        
        %% Create 3-dimensional time-series matrix in which
        % each frame in time contains the x- and y-coordinates.
        % The z-coordinate is composed of each frame.
        nFrames = length(ind);
        row = 256;
        col = 63;
        
        TimeMatrix = zeros(row,col,nFrames);
        TimeMatrix(:,:,1) = PRO(1:256,:);
        % if there is a double row of NaN spacing between timeseries
        if ind(1) == ind(2) - 1
            ind = find(isnan(PRO(:,1)));
            ind(mod(ind,2) == 0) = [];
            nFrames = length(ind);
            for i = 2:nFrames
                TimeMatrix(:,:,i) = PRO(ind(i-1)+2:ind(i)-1,:);
            end
        else % if single row of NaN spacing between timeseries
            for i = 2:nFrames
                TimeMatrix(:,:,i) = PRO(ind(i-1)+1:ind(i)-1,:);
            end
        end
        
          %% compute center of pressure
        center_of_pressure = zeros(nFrames,2);
        for i = 1:nFrames
            CC = bwconncomp(TimeMatrix(:,:,i));
            CCC = labelmatrix(CC);
            CCC(CCC~=0) = 1;
            Cent = regionprops(CCC,'centroid');
            if isempty(Cent) == 1
                center_of_pressure(i,1:2) = [0,0];
            else
                center_of_pressure(i,1:2) = Cent.Centroid;
            end
        end
        
%         SUMTM = sum(TimeMatrix(:,:,:),3); 
%         PlantarData = contour(SUMTM); hold on;
%         RSScan = plot(Alex_CoP_Test.center_of_pressure(:,1), Alex_CoP_Test.center_of_pressure(:,2),'b.'); 
%         SelfCompute = plot(center_of_pressure(:,1), center_of_pressure(:,2),'r.'); 
%         legend({'PlantarData','RSScan','SelfCompute'}); 
        
        %% Save Data into Structure format
        AnalyzePPImages(z).time_matrix = TimeMatrix;
%         AnalyzePPImages(z).subject_name = subject_name;
%         AnalyzePPImages(z).study_date = study_date;
        AnalyzePPImages(z).center_of_pressure = center_of_pressure;
    end
end

%% save aggregate data as mat file
if exist('SaveName') == 0
    prompt_msg = 'Please enter name of mat file to save to reduce import next time';
    prompt = {prompt_msg};
    dlg_title = 'Export Data';
    num_lines = 1;
    def = {''};
    SaveName = strcat(char(cell2mat( inputdlg(prompt,dlg_title,num_lines,def) )), '_RawData.mat');
end

save(SaveName, 'AnalyzePPImages');
end

