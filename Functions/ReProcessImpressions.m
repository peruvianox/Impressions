%% IMPRESSIONS
% Plantar Pressure Analysis
% Version 1.1       February 6, 2019

% Impressions will input raw plantar pressure data from a pressure plate and calculate a
% variety of metrics related to walking. Impressions requires version
% Matlab version 2017b or newer and the image processing toolbox.

%                                       Inputs:
% Time-Series Matrix of raw pressure data (aka Entire Plate Roll Off)
% This data file should contain a 2-D (length by width) matrix of pressure or force
% varables arranged in a consecutive time series. A quality recording will contain
% 2-5 steps along the length of the mat and at least 4 trials for consistency.


%                                       Outputs:
% Plantar Pressure Report
% An excel spreadsheet containing the most relevant figures, metrics, and
% ratios from the plantar pressure analysis.
%
% Saved .mat file
% This file contains all the saved data from the matlab workspace so that
% the data may be saved in an analysis-friendly version for future
% reference.

% GIF movie files for each


% Supported PP Mats: RS Scan 2m, Novel emedXL
% These files can initially be loaded as generic files stright from RS Scan (slowest),
%  or Novel .lst files from manual conversion, or as
% pre-processed .mat files (fastest).

% Author: Ricky Pimentel
% Center for Gait and Movement Analysis, Children's Hospital Colorado

clear;
close all;
clc;
warning off;
dbstop if error;
addpath(genpath('Functions'));
addpath(genpath('Data'));

[filename, path] = uigetfile('*.mat','Select .mat file(s) for PP Analysis','MultiSelect','On');
addpath(genpath(path));
load(filename); 
close all

%% define input options
% clc; disp('Obtaining input selections');
% [SubjDemo, PPSettings, PPPlots] = PPSelectOptions;
% if isnan(SubjDemo(5).Choice) % if foot length not defined, divide height by 6
%     SubjDemo(5).Choice = SubjDemo(3).Choice / 6;
% end
% 
% % define folder for saving if not already input
if isempty(SubjDemo(1).Choice)
    folder = inputdlg('Input the folder namer for data saving');
else
    folder = SubjDemo(1).Choice;
end
mkdir(folder);
% 
% %% select and load PP data
% clc; disp('Selection and loading of PP data');
% if strcmp(PPSettings.LoadNew,'Raw') == 1
%     % FYI: this part takes at least 12 seconds per trial to load. Elapsed times per trial will
%     % be shown in the command window for reference.
%     [AnalyzePPImages, PPSettings.PPMatType, filename] = plantar_load(0, PPSettings.PPMatType);
% elseif strcmp(PPSettings.LoadNew,'Mat')
%     % If loading mat files, select
%     [filename, ~] = uigetfile('*.mat','Select .mat file(s) for PP Analysis','MultiSelect','On');
%     if iscell(filename)
%         NumTrials = length(filename);
%         for q = 1:NumTrials % load Loop for multiple trials
%             load(filename{q});
%         end
%         for q = 1:NumTrials
%             center_of_pressure = zeros(size(AnalyzePPImages(q).time_matrix, 3),2);
%             for i = 1:size(AnalyzePPImages(q).time_matrix, 3)
%                 CC = bwconncomp(PPTrials(q).time_matrix(:,:,i));
%                 CCC = labelmatrix(CC);
%                 CCC(CCC~=0) = 1;
%                 Cent = regionprops(CCC,'centroid');
%                 if isempty(Cent) == 1
%                     AnalyzePPImages(q).center_of_pressure(i,1:2) = [0,0];
%                 else
%                     AnalyzePPImages(q).center_of_pressure(i,1:2) = Cent.Centroid;
%                 end
%             end
%         end
%     else % if multiple trials saved on 1 mat file
%         load(filename);
%         NumTrials = length(AnalyzePPImages);
%         for q = 1:NumTrials
%             center_of_pressure = zeros(size(AnalyzePPImages(q).time_matrix, 3),2);
%             for i = 1:size(AnalyzePPImages(q).time_matrix, 3)
%                 CC = bwconncomp(AnalyzePPImages(q).time_matrix(:,:,i));
%                 CCC = labelmatrix(CC);
%                 CCC(CCC~=0) = 1;
%                 Cent = regionprops(CCC,'centroid');
%                 if isempty(Cent) == 1
%                     AnalyzePPImages(q).center_of_pressure(i,1:2) = [0,0];
%                 else
%                     AnalyzePPImages(q).center_of_pressure(i,1:2) = Cent.Centroid;
%                 end
%             end
%         end
%     end
% end
% 
% %% Plantar Pressure Plate Dimensions
% [m,n, ~] = size(AnalyzePPImages(1).time_matrix);
% if m == 256 && n == 63
%     PPSettings.PPMatType = 'RSScan';
%     PPDim.cmLength = 196;
%     PPDim.cmWidth = 34;
%     %     PPMeasure.Device = 'RSScan';
%     PPDim.Length = m;
%     PPDim.Width = n;
% elseif m == 288 && n == 88
%     PPSettings.PPMatType = 'Novel';
%     PPDim.cmLength = 144;
%     PPDim.cmWidth = 44;
%     %     PPMeasure.Device = 'Novel Emed XL';
%     PPDim.Length = 288;
%     PPDim.Width = 88;
% end
% Adj.Length = PPDim.cmLength/PPDim.Length;
% Adj.Width  = PPDim.cmWidth/PPDim.Width;
% Adj.Area = Adj.Width * Adj.Length; % Area adjustment from cells to cm^2
% Adj.Pressure = Adj.Area; % Pressure Adjustment = N* Adj.Area, UNITS = Newtons
% 
% clearvars CC CCC Cent center_of_pressure i m n q
% 
% %% Initial PP quality check
% % Sets the units of your root object (screen) to pixels
% set(0,'units','pixels');
% Pix_SS = get(0,'screensize');
% Pix_SS(1:2) = [];
% DisplayDim = [40, 40, Pix_SS(1)*.66, Pix_SS(2)*.8];
% if strcmp(PPSettings.InitialCheck, 'Yes')
%     Test = figure;
%     set(Test, 'Position', DisplayDim);
%     for i = 1:length(AnalyzePPImages)
%         subplot(1, length(AnalyzePPImages), i);
%         A  = sum(AnalyzePPImages(i).time_matrix, 3);
%         contour(A,500,'LineWidth', 1);
%         ax = gca;
%         ax.XTick = [];
%         ax.YTick = [];
%         title(strcat('Trial ', num2str(i)));
%     end
%     saveas(Test, strcat(folder,'/', folder, '_Raw.png'));
%     Proceed = questdlg('Do all trials contain at least 2 consecutive steps on the mat along its length?');
%     if strcmp(Proceed, 'No')
%         NumDynTrials = length(AnalyzePPImages); % determine # of trials
%         if NumDynTrials == 1
%             ListStr = {'1'};
%         elseif NumDynTrials == 2
%             ListStr = {'1','2'};
%         elseif NumDynTrials == 3
%             ListStr = {'1','2','3'};
%         elseif NumDynTrials == 4
%             ListStr = {'1','2','3','4'};
%         elseif NumDynTrials == 5
%             ListStr = {'1','2','3','4','5'};
%         elseif NumDynTrials == 6
%             ListStr = {'1','2','3','4','5','6'};
%         elseif NumDynTrials == 7
%             ListStr = {'1','2','3','4','5','6','7'};
%         elseif NumDynTrials == 8
%             ListStr = {'1','2','3','4','5','6','7','8'};
%         elseif NumDynTrials == 9
%             ListStr = {'1','2','3','4','5','6','7','8','9'};
%         elseif NumDynTrials == 10
%             ListStr = {'1','2','3','4','5','6','7','8','9','10'};
%         end
%         
%         % option to edit trials
%         EditTrials = questdlg('Would you like to edit any trials?');
%         if strcmp(EditTrials, 'Yes')
%             ToEdit = listdlg('PromptString','Select Trials to edit','ListString', ListStr);
%             
%             for i = ToEdit
%                 uiwait(msgbox(['Select area(s) to ignore in trial ',num2str(i)]));
%                 j = 1;
%                 EditMore = 'Yes';
%                 while strcmp(EditMore, 'Yes') == 1
%                     subplot(1, length(AnalyzePPImages), i);
%                     h = imrect; % select areas to delete
%                     Delete(j).Coord = round(getPosition(h)); % save area coordinates
%                     if Delete(j).Coord(1) < 1 % QA for coordinates to ensure they are on the mat and indexible
%                         Delete(j).Coord(1) = 1;
%                     end
%                     if Delete(j).Coord(2) < 1
%                         Delete(j).Coord(2) = 1;
%                     end
%                     AddX = Delete(j).Coord(1)+Delete(j).Coord(3);
%                     if AddX > PPDim.Width
%                         AddX = PPDim.Width;
%                     end
%                     AddY = Delete(j).Coord(2)+Delete(j).Coord(4);
%                     if AddY> PPDim.Length
%                         AddY = PPDim.Length;
%                     end
%                     AnalyzePPImages(i).time_matrix(Delete(j).Coord(2):AddY, Delete(j).Coord(1):AddX,:) = 0; % "zero" areas to delete
%                     j = j + 1;
%                     subplot(1, length(AnalyzePPImages), i); % re-plot with areas deleted
%                     A = sum(AnalyzePPImages(i).time_matrix, 3);
%                     contour(A,500,'LineWidth', 1);
%                     ax = gca;
%                     ax.XTick = [];
%                     ax.YTick = [];
%                      title(strcat('Trial ', num2str(i)));
%                     Question = 'Would you like to remove more sections?';
%                     EditMore = questdlg(Question,'Remove areas','Yes','No','No');
%                     clearvars h
%                 end
%             end
%         end
%         
%         % option to delete trials
%         DeleteTrials = questdlg('Would you like to delete any rogue trials?');
%         if strcmp(DeleteTrials,'Yes')
%             ToDel = listdlg('PromptString','Select Trials to delete','ListString', ListStr);
%             AnalyzePPImages(ToDel) = [];
%             NumTrials = length(AnalyzePPImages);
%         end
%         clf;
%         for i = 1:length(AnalyzePPImages)
%             subplot(1, length(AnalyzePPImages), i);
%             A  = sum(AnalyzePPImages(i).time_matrix, 3);
%             contour(A, 500,'LineWidth', 1);
%             ax = gca;
%             ax.XTick = [];
%             ax.YTick = [];
%             title(strcat('Trial ', num2str(i)));
%         end
%         % if no changes, exit script and try again
%         if strcmp(EditTrials, 'No') && strcmp(DeleteTrials, 'No')
%             uiwait(msgbox('Please select a different set of trials and try again'));
%             return
%         end
%         saveas(Test, strcat(folder,'/', folder, '_Edited.png'));
%     end
%     close all;
% end
% % Assign edited data to a new variable, PPTrials
% for q = 1:length(AnalyzePPImages)
%     [~, ~, p] = size( AnalyzePPImages(q).time_matrix);
%     for i = 1:p
%         if sum(sum(AnalyzePPImages(q).time_matrix(:,:,i))) > 0
%             PPTrials(q).TM(:,:,i) = AnalyzePPImages(q).time_matrix(:,:,i);
%             PPTrials(q).CoP(i,:) = AnalyzePPImages(q).center_of_pressure(i,:);
%         end
%     end
%     PPTrials(q).SumTM = sum(PPTrials(q).TM,3);
% end
% clearvars ListStr ToDel NumDynTrials DeleteTrials EditTrials Delete EditMore
% NumTrials = length(PPTrials);
% 
% %% Subject Qualities
% if exist('SubjDemo','var')
%     Subject.Age = SubjDemo(2).Choice;
%     Subject.Height = SubjDemo(3).Choice;
%     Subject.Mass = SubjDemo(4).Choice;
%     Subject.FootLength = SubjDemo(5).Choice;
% end
% % if subject demographics aren't imput, prompt for their values
% if isempty(Subject.Age) || isempty(SubjDemo(2).Choice)
%     Subject.Age = inputdlg('What is the age of the subject in years?');
% end
% if isempty(Subject.Mass) || isempty(SubjDemo(3).Choice)
%     Subject.Mass = inputdlg('What is the mass of the subject in kg?');
% end
% if isempty(Subject.Height) || isempty(SubjDemo(4).Choice)
%     Subject.Height = inputdlg('What is the height of the subject in cm?');
% end
% if isempty(Subject.FootLength) || isempty(SubjDemo(5).Choice)
%     Subject.FootLength = str2double(inputdlg('What is the average footlength of the subject in cm?'));
% end
% 
% %% determine which trials need to be flipped
% for q = 1:NumTrials % loop through trials to determine if they need to be flipped
%     % find start of PPs
%     p = size(PPTrials(q).TM,3);
%     for i = 1:p
%         if sum(sum(PPTrials(q).TM(:,:,i))) > 0
%             PPTrials(q).Start = i;
%             break
%         end
%     end
%     for i = 1:p
%         if sum(sum(PPTrials(q).TM(:,:,p-i+1))) > 0
%             PPTrials(q).End = p-i+1;
%             break
%         end
%     end
%     clearvars p i
%     
%     PPTrials(q).Type = 'Dynamic';
%     if mean(PPTrials(q).CoP(PPTrials(q).Start:PPTrials(q).Start+10,2)) < mean(PPTrials(q).CoP(PPTrials(q).End-10:PPTrials(q).End,2)) % if first cop value is less than the last
%         PPTrials(q).Flip = 0;
%         PPTrials(q).SumTM = fliplr(PPTrials(q).SumTM);
%         PPTrials(q).TM = fliplr(PPTrials(q).TM);
%         PPTrials(q).CoP(:,1) =  PPDim.Width - PPTrials(q).CoP(:,1);
%         PPTrials(q).CoP(:,2) = PPTrials(q).CoP(:,2);
%     else
%         PPTrials(q).Flip = 1;
%         PPTrials(q).SumTM = flipud(PPTrials(q).SumTM);
%         PPTrials(q).TM = flipud(PPTrials(q).TM);
%         PPTrials(q).CoP(:,1) = PPTrials(q).CoP(:,1);
%         PPTrials(q).CoP(:,2) = PPDim.Length - PPTrials(q).CoP(:,2);
%     end
% end
% 
% %% Static Trial Analysis
% for i = 1:NumTrials
%     LogDyn(i)  = strcmp(PPTrials(i).Type, 'Dynamic');
% end
% % if strcmp(PPSettings.StaticAnalysis, 'Yes') % split trials into static and dynamic
% %     StaticPPTrials = PPTrials(~LogDyn);
% % end
% DynamicPPTrials = PPTrials(LogDyn);
% NumDynTrials = length(DynamicPPTrials);
% 
% %% Ensure Quality foot pressures
% FLest =  1.1 * Subject.FootLength / Adj.Length;
% Regions = cell(1,NumDynTrials);
% % Identify foot strikes for Spatio Temporal Measurements
% for i = 1:NumDynTrials
%     [Regions{i}, ToeDrag] = SpatioTemporalPP(DynamicPPTrials(i).TM, FLest, PPSettings, PPDim, PPTrials(i).CoP);
% end
% 
% %% Check L/R classifications
% clc; disp('Left/Right Classifications');
% LRTest = figure;
% set(LRTest, 'Position', DisplayDim);
% for i = 1:length(DynamicPPTrials)
%     subplot(1, length(DynamicPPTrials), i);
%     contour(DynamicPPTrials(i).SumTM, 100);
%     title(strcat('Trial ', num2str(i)));
%     for j = 1:length(Regions{i})
%         if strcmp(Regions{i}(j).Side, 'Left')
%             rectangle('Position', [Regions{i}(j).BoundingBox], 'EdgeColor','r');
%         else
%             rectangle('Position', [Regions{i}(j).BoundingBox], 'EdgeColor','g');
%         end
%     end
%     set(gca, 'XTick', []);
%     set(gca, 'YTick', []);
% end
% ReClassify = questdlg('LEFT and RIGHT feet are defined using red (LEFT) and green (RIGHT) boxes. Would you like to re-classify any of the LEFTS or RIGHTS?');
% while strcmp(ReClassify, 'Yes')
%     if NumDynTrials == 1
%         ListStr = {'1'};
%     elseif NumDynTrials == 2
%         ListStr = {'1','2'};
%     elseif NumDynTrials == 3
%         ListStr = {'1','2','3'};
%     elseif NumDynTrials == 4
%         ListStr = {'1','2','3','4'};
%     elseif NumDynTrials == 5
%         ListStr = {'1','2','3','4','5'};
%     elseif NumDynTrials == 6
%         ListStr = {'1','2','3','4','5','6'};
%     elseif NumDynTrials == 7
%         ListStr = {'1','2','3','4','5','6','7'};
%     elseif NumDynTrials == 8
%         ListStr = {'1','2','3','4','5','6','7','8'};
%     elseif NumDynTrials == 9
%         ListStr = {'1','2','3','4','5','6','7','8','9'};
%     elseif NumDynTrials == 10
%         ListStr = {'1','2','3','4','5','6','7','8','9','10'};
%     end
%     Trials2ReClass = listdlg('PromptString','Select trial(s) to reclassify sides','ListString', ListStr);
%     % edit foot progression angles
%     for i = 1:length(Trials2ReClass)
%         subplot(1, NumDynTrials, Trials2ReClass(i));
% %         title(strcat('Trial ',  Trials2ReClass(i)));
%         % LEFT side
%         LeftAdj = questdlg(['For Trial ' num2str(Trials2ReClass(i)) ,'Would you like to define any LEFT strikes?']);
%         if strcmp(LeftAdj,'Yes') == 1
%             LeftDefine = str2double(cell2mat(inputdlg('How many LEFT feet would you like to define?')));
%             for j = 1:LeftDefine
%                 h = imrect;
%                 Coords(Trials2ReClass(i)).Left(1:4) = getPosition(h);
%                 for k = 1:length(Regions{Trials2ReClass(i)})
%                     if Regions{Trials2ReClass(i)}(k).Centroid(1) > Coords(Trials2ReClass(i)).Left(1) && Regions{Trials2ReClass(i)}(k).Centroid(1) < Coords(Trials2ReClass(i)).Left(1) + Coords(Trials2ReClass(i)).Left(3)
%                         if Regions{Trials2ReClass(i)}(k).Centroid(2) > Coords(Trials2ReClass(i)).Left(2) && Regions{Trials2ReClass(i)}(k).Centroid(2) < Coords(Trials2ReClass(i)).Left(2) + Coords(Trials2ReClass(i)).Left(4)
%                             Regions{Trials2ReClass(i)}(k).Side = 'Left';
%                         end
%                     end
%                 end
%             end
%         end
%         
%         % RIGHT side
%         RightAdj = questdlg(['For Trial ' num2str(Trials2ReClass(i)) ,'Would you like to define any RIGHT strikes?']);
%         if strcmp(RightAdj,'Yes') == 1
%             RightDefine = str2double(cell2mat(inputdlg('How many RIGHT feet would you like to define?')));
%             for j = 1:RightDefine
%                 h = imrect;
%                 Coords(Trials2ReClass(i)).Right(1:4) = getPosition(h);
%                 for k = 1:length(Regions{Trials2ReClass(i)})
%                     if Regions{Trials2ReClass(i)}(k).Centroid(1) > Coords(Trials2ReClass(i)).Right(1) && Regions{Trials2ReClass(i)}(k).Centroid(1) < Coords(Trials2ReClass(i)).Right(1) + Coords(Trials2ReClass(i)).Right(3)
%                         if Regions{Trials2ReClass(i)}(k).Centroid(2) > Coords(Trials2ReClass(i)).Right(2) && Regions{Trials2ReClass(i)}(k).Centroid(2) < Coords(Trials2ReClass(i)).Right(2) + Coords(Trials2ReClass(i)).Right(4)
%                             Regions{Trials2ReClass(i)}(k).Side = 'Right';
%                         end
%                     end
%                 end
%             end
%         end
%     end
%     clf;
%     for i = 1:length(DynamicPPTrials)
%         subplot(1, length(DynamicPPTrials), i);
%         contour(DynamicPPTrials(i).SumTM, 100);
%         title(strcat('Trial ', num2str(i)));
%         set(gca, 'XTick', []);
%         set(gca, 'YTick', []);
%         for j = 1:length(Regions{i})
%             if strcmp(Regions{i}(j).Side, 'Left')
%                 rectangle('Position', [Regions{i}(j).BoundingBox], 'EdgeColor','r');
%             else
%                 rectangle('Position', [Regions{i}(j).BoundingBox], 'EdgeColor','g');
%             end
%         end
%     end
%     ReClassify = questdlg('Would you like to re-classify right and left feet again?');
% end
% clearvars Trials2ReClass RightAdj RightDefine LeftAdj LeftDefine Coords
% close;
% 
% %% Load C3D Data if desired
if strcmp(PPSettings.C3DInput, 'Yes') == 1
    clc; disp('Loading C3D data');
    [C3Ddata,TrimSpat, DynamicPPTrials] = PPInputC3D(PPDim, DynamicPPTrials, Subject, Regions);
end
% 
% %% Center of Pressure Line for each foot;
% % create subset of foot strikes for each new CoP
% for i = 1:NumDynTrials
%     for j = 1:length(Regions{1,i})
%         x = ceil(Regions{1,i}(j).BoundingBox(1));
%         y = ceil(Regions{1,i}(j).BoundingBox(2));
%         xp = x+Regions{1,i}(j).BoundingBox(3);
%         yp = y+Regions{1,i}(j).BoundingBox(4);
%         if xp > PPDim.Width
%             xp = PPDim.Width;
%         end
%         if yp > PPDim.Length
%             yp = PPDim.Length;
%         end
%         s =  Regions{1,i}(j).Strike;
%         f =  Regions{1,i}(j).Off;
%         if s == 0
%             break
%         end
%         Regions{1,i}(j).Step = DynamicPPTrials(i).TM(y:yp, x:xp, s:f);
%     end
% end
% 
% % calculate CoP
% for i = 1:NumDynTrials
%     for j = 1:length(Regions{1,i})
%         [h,w, L] = size(Regions{1,i}(j).Step(:,:,:)); % determine size of matrix
%         CoP = zeros(L,2);
%         for k = 1:L
%             % manual CoP calculation
%             Cx = sum(sum(Regions{1,i}(j).Step(:,:,k)).*(1:w))/sum(sum(Regions{1,i}(j).Step(:,:,k))); % x location centroid
%             Cy = sum(sum(Regions{1,i}(j).Step(:,:,k),2)'.*(1:h))/sum(sum(Regions{1,i}(j).Step(:,:,k),2)'); % y location centroid
%             CoP(k,:) = [Cx,Cy];
%             clearvars Cx Cy
%         end
%         Regions{1,i}(j).StepCoP = CoP; % define CoP for each step
%         Regions{1,i}(j).StepCoP66 = CoP(1:round(0.66*L), :);
%         clearvars CoP
%         Regions{1,i}(j).StepSum(:,:) = sum(Regions{1,i}(j).Step(:,:,:),3);  % create sum image
%     end
% end
% clearvars w x xp y yp q h i j k L f Cx Cy Reg Test LogDyn
% 
% %% Identify pressure peaks
% clc; disp('Identifying CoP and pressure regional peaks');
% for i = 1:NumDynTrials
%     for j = 1:length(Regions{1,i})
%         % identify longitudinal peaks in pressure
%         [h,~] = size(Regions{1,i}(j).StepSum);
%         for k = 1:h
%             Regions{1,i}(j).LongSum(k) = sum(Regions{1,i}(j).StepSum(k,:));
%         end
%         % set thresholds for peaks
%         if strcmp(PPSettings.PPMatType, 'Novel')
%             PkThresh = 10000;
%         elseif strcmp(PPSettings.PPMatType, 'RSScan')
%             PkThresh = 25;
%         end
%         % find said peaks
%         [Regions{1,i}(j).LongPks, Regions{1,i}(j).LongLocs] = findpeaks(Regions{1,i}(j).LongSum, 'MinPeakProminence',PkThresh);
%         % if no peaks, signify to exclude IP and HC methods
%         if isempty(Regions{1,i}(j).LongPks)
%             Block.IPs = 1;
%             Block.HCs = 1;
%         else
%             Regions{1,i}(j).SurfPks = FastPeakFind(Regions{1,i}(j).StepSum);
%             % classify peaks for labelling
%             l = 1;
%             for k = 1:length(Regions{1,i}(j).SurfPks)
%                 if mod(k,2) == 1 % if odd
%                     Mat(l,1) = Regions{1,i}(j).SurfPks(k);
%                 else
%                     Mat(l,2) = Regions{1,i}(j).SurfPks(k);
%                     l = l+1;
%                 end
%             end
%             Regions{1,i}(j).SurfPks = [];
%             Regions{1,i}(j).SurfPks = Mat;
%             
%             % Determine initial arch index for improved prog line detection
%             if length(Regions{1,i}(j).LongLocs) > 1
%                 Block.IPs(i,j) = 0;
%                 % between 1st and 2nd peaks
%                 Regions{1,i}(j).MeanArchPress = nanmean(Regions{1,i}(j).LongSum(Regions{1,i}(j).LongLocs(1): Regions{1,i}(j).LongLocs(2)));
%                 Regions{1,i}(j).ArchRatio = Regions{1,i}(j).MeanArchPress / mean(Regions{1,i}(j).LongPks(1:2));
%             else
%                 Block.IPs(i,j) = 1;
%             end
%             
%             % Identify heel peak
%             [~,ind] = min(Regions{1,i}(j).SurfPks(:,2));
%             if ind >= 15
%                 Block.HCs(i,j) = 1;
%             else
%                 Regions{1,i}(j).HeelPt = Regions{1,i}(j).SurfPks(ind,:);
%                 Block.HCs(i,j) = 0;
%             end
%             % uncomment to plot peaks
%             %                 figure;
%             %                 subplot(121);
%             %                 plot(Regions{1,i}(j).LongSum); hold on;
%             % %                 plot([Regions{1,i}(j).LongLocs(1), Regions{1,i}(j).LongLocs(2)], [Regions{1,i}(j).MeanArchPress, Regions{1,i}(j).MeanArchPress],'k');
%             % %                 plot(Regions{1,i}(j).LongLocs(:), Regions{1,i}(j).LongPks(:),  'om');
%             %
%             %                 subplot(122);
%             %                 contour(Regions{1,i}(j).StepSum, 25)
%             % %                 surf(Regions{1,i}(j).StepSum);
%             %                 hold on;
%             %                 plot(Regions{1,i}(j).SurfPks(:,1), Regions{1,i}(j).SurfPks(:,2), 'ok');
%         end
%         clearvars Mat k l h w ind
%     end
% end
% % determine legitimacy of IP method
% if sum(sum(Block.IPs)) > 0
%     Block.IP = 1;
% else
%     Block.IP = 0;
% end
% % determine legitimacy of HC method
% if sum(sum(Block.HCs)) > 0
%     Block.HC = 1;
% else
%     Block.HC = 0;
% end
% 
% %% use CoP to determine foot progression angle
% for i = 1:NumDynTrials
%     % foot prog for full CoP trajectory
%     for j = 1:length(Regions{1,i})
%         % translate centroid to smaller window
%         Regions{1,i}(j).TranCent = [Regions{1,i}(j).Centroid(1) - Regions{1,i}(j).BoundingBox(1), Regions{1,i}(j).Centroid(2) - Regions{1,i}(j).BoundingBox(2)];
%         
%         % CoP trendlining
%         indX = isnan(Regions{1,i}(j).StepCoP(:,2));
%         indY = isnan(Regions{1,i}(j).StepCoP(:,1));
%         IND = indX + indY;
%         Regions{1,i}(j).StepCoP(IND>0, :) = [];
%         x = Regions{1,i}(j).StepCoP(:,2); % x values in
%         y = Regions{1,i}(j).StepCoP(:,1);  % y values in
%         p = polyfit(x,y,1); % first order polynomial
%         xFit = linspace(0,size(Regions{1,i}(j).Step, 1),100); % fit x
%         yFit = polyval(p, xFit); % evaluate y for every x
%         Regions{1,i}(j).CoP.Cent = [nanmean(y) + Regions{1,i}(j).BoundingBox(1), nanmean(x) + Regions{1,i}(j).BoundingBox(2)];
%         % uncomment to plot individual foot strikes in loop
%         %         figure;
%         %         contour(Regions{1,i}(j).StepSum, 20); hold on;
%         %         plot(Regions{1,i}(j).StepCoP(:,1), Regions{1,i}(j).StepCoP(:,2),'k.');
%         %         plot(yFit, xFit, 'r*')
%         %         axis equal;
%         if strcmp(Regions{1,i}(j).Side, 'Left')
%             Regions{1,i}(j).CoP.FootProg =  -180/pi*atan2(yFit(end) - yFit(1), xFit(end) - xFit(1));
%         else
%             Regions{1,i}(j).CoP.FootProg =  180/pi*atan2(yFit(end) - yFit(1), xFit(end) - xFit(1));
%         end
%         clearvars x y p xFit yFit
%     end
%     
%     % foot prog for 66% of CoP trajectory
%     for j = 1:length(Regions{1,i})
%         % CoP trendlining
%         indX = isnan(Regions{1,i}(j).StepCoP66(:,2));
%         indY = isnan(Regions{1,i}(j).StepCoP66(:,1));
%         IND = indX + indY;
%         Regions{1,i}(j).StepCoP66(IND>0, :) = [];
%         x = Regions{1,i}(j).StepCoP66(:,2); % x values in
%         y = Regions{1,i}(j).StepCoP66(:,1);  % y values in
%         p = polyfit(x,y,1); % first order polynomial
%         xFit = linspace(0,size(Regions{1,i}(j).Step, 1),100); % fit x
%         yFit = polyval(p, xFit); % evaluate y for every x
%         Regions{1,i}(j).CoP66.Cent = [nanmean(y)  + Regions{1,i}(j).BoundingBox(1), nanmean(x) + Regions{1,i}(j).BoundingBox(2)];
%         % uncomment to plot individual foot strikes in loop
%         %                 figure;
%         %                 contour(Regions{1,i}(j).StepSum, 20); hold on;
%         %                 plot(Regions{1,i}(j).StepCoP66(:,1), Regions{1,i}(j).StepCoP66(:,2),'k.');
%         %                 plot(yFit, xFit, 'r*')
%         %                 axis equal;
%         if strcmp(Regions{1,i}(j).Side, 'Left')
%             Regions{1,i}(j).CoP66.FootProg =  -180/pi*atan2(yFit(end) - yFit(1), xFit(end) - xFit(1));
%         else
%             Regions{1,i}(j).CoP66.FootProg =  180/pi*atan2(yFit(end) - yFit(1), xFit(end) - xFit(1));
%         end
%         clearvars x y p xFit yFit
%     end
%     
%     % foot prog Inter-peak CoP trajectory
%     if Block.IP == 0
%         for j = 1:length(Regions{1,i})
%             % CoP trendlining
%             S = Regions{1,i}(j).LongLocs(1);
%             F = Regions{1,i}(j).LongLocs(2);
%             x = Regions{1,i}(j).StepCoP(S:F,2); % x values in
%             y = Regions{1,i}(j).StepCoP(S:F,1);  % y values in
%             p = polyfit(x,y,1); % first order polynomial
%             xFit = linspace(0,size(Regions{1,i}(j).Step, 1),100); % fit x
%             yFit = polyval(p, xFit); % evaluate y for every x
%             Regions{1,i}(j).IP.Cent = [nanmean(y) + Regions{1,i}(j).BoundingBox(1), nanmean(x)  + Regions{1,i}(j).BoundingBox(2)];
%             % uncomment to plot individual foot strikes in loop
%             %                 figure;
%             %                 subplot(121);
%             %                 contour(Regions{1,i}(j).StepSum, 20); hold on;
%             %                  plot(Regions{1,i}(j).StepCoP(S:F,1), Regions{1,i}(j).StepCoP(S:F,2),'ko');
%             %                 plot(yFit, xFit, 'r.')
%             %                 axis equal;
%             %                 subplot(122);
%             %                 plot(Regions{1,i}(j).LongSum); hold on;
%             %                 plot(Regions{1,i}(j).LongLocs, Regions{1,i}(j).LongPks, 'om');
%             if strcmp(Regions{1,i}(j).Side, 'Left')
%                 Regions{1,i}(j).IP.FootProg =  -180/pi*atan2(yFit(end) - yFit(1), xFit(end) - xFit(1));
%             else
%                 Regions{1,i}(j).IP.FootProg =  180/pi*atan2(yFit(end) - yFit(1), xFit(end) - xFit(1));
%             end
%             clearvars x y p xFit yFit
%         end
%     end
%     
%     % foot prog heel peak to centroid
%     if Block.HC == 0
%         for j = 1:length(Regions{1,i})
%             % distances from heel point to centroid
%             dx = Regions{1,i}(j).TranCent(1) - Regions{1,i}(j).HeelPt(1);
%             dy = Regions{1,i}(j).TranCent(2) - Regions{1,i}(j).HeelPt(2);
%             % uncomment to plot individual footstrikes in loop
%             %         figure;
%             %         contour(Regions{1,i}(j).StepSum, 20); hold on;
%             %         plot(Regions{1,i}(j).StepCoP(:,1), Regions{1,i}(j).StepCoP(:,2),'k.');
%             %         plot(Regions{1,i}(j).TranCent(1), Regions{1,i}(j).TranCent(2), 'ok');
%             %         plot(Regions{1,i}(j).HeelPt(1), Regions{1,i}(j).HeelPt(2), 'ok');
%             %         axis equal;
%             if strcmp(Regions{1,i}(j).Side, 'Left')
%                 Regions{1,i}(j).HC.FootProg =  -180/pi*atan2(dx, dy);
%             else
%                 Regions{1,i}(j).HC.FootProg =  180/pi*atan2(dx, dy);
%             end
%             clearvars x y p xFit yFit dx dy dx2 dy2
%         end
%     end
% end
% 
% %% Identify and correct poor pressure trials
% clc; disp('Checking for invalid pressures');
% clearvars ToDel
% w = 1; 
% % identify poor pressure recordings
% % if the estimated or measured foot length is greater than 125% of the actual
% for i = 1:NumDynTrials
%     L = 0;
%     R = 0;
%     for j = 1:length(Regions{1,i})
%         if FLest > Regions{1,i}(j).MajorAxisLength*1.25
%             Regions{1,i}(j).PoorPressure = 1;
%             if strcmp(Regions{1,i}(j).Side,'Left')
%                 L = L+1;
%                 DynamicPPTrials(i).PoorPressure.Left(L) = 1;
%             else
%                 R = R+1;
%                 DynamicPPTrials(i).PoorPressure.Right(R) = 1;
%             end
%         else
%             Regions{1,i}(j).PoorPressure = 0;
%             if strcmp(Regions{1,i}(j).Side,'Left')
%                 L = L+1;
%                 DynamicPPTrials(i).PoorPressure.Left(L) = 0;
%             else
%                 R = R+1;
%                 DynamicPPTrials(i).PoorPressure.Right(R) = 0;
%             end
%         end
%     end
% end
% % create logical variable that identifies which trials have poor pressures
% indPoorPressureTrial = zeros(1,NumDynTrials);
% for i = 1:NumDynTrials
%     if sum([Regions{1,i}.PoorPressure]) > 0 %sum(PoorPressure) > 0
%         indPoorPressureTrial(i) = 1;
%     else
%         indPoorPressureTrial(i) = 0;
%     end
%     clearvars PoorPressure
% end
% if sum(indPoorPressureTrial) > 0
%     % block automated methods that use peak definition if pressures are invalid
%     Block.IP = 1;
%     Block.HC = 1;
%     uiwait(msgbox(['Impressions has detected ' num2str(sum(indPoorPressureTrial)) 'trial(s) with incomplete pressures.']));
%     % plot the contours of the poor pressure trials
%     PoorPress = find(indPoorPressureTrial);
%     figure('Position',  [DisplayDim(1) DisplayDim(2) 250*length(PoorPress) DisplayDim(4)]);
%     j = 1;
%     for i = 1:length(PoorPress)
%         subplot(1,sum(indPoorPressureTrial),i);
%         contour(DynamicPPTrials(PoorPress(i)).SumTM,100);
%         title(['Trial' num2str(PoorPress(i))]);
%         set(gca, 'XTick', []);
%         set(gca, 'YTick', []);
%         for k = 1: length(Regions{1,PoorPress(i)})
%             if Regions{1,PoorPress(i)}(k).PoorPressure == 1
%                 rectangle('Position', Regions{1,PoorPress(i)}(k).BoundingBox, 'EdgeColor','r');
%             end
%         end
%         axis equal;
%         j = j+1;
%     end
%     % instructions for adding/deleting foot pressures from
%     uiwait(msgbox(['This figure displays the trial(s) with invalid pressures.  Identify the steps you wish to edit or delete.  Ensure final selections are continuous.'], 'Instructions'));
%     % mention or determine how many pressures are detected
%     for i = 1:sum(indPoorPressureTrial)
%         % choose to add pressure, delete pressure, or delete trial
%         DelStrikes = questdlg(['For Trial ' num2str(PoorPress(i)) ', would you like to keep or  delete pressure(s)?'], 'Keep/Edit/Delete','Keep Pressure(s)','Delete Pressure(s)','Delete Whole Trial','Keep Presure(s)');
%        
%         if strcmp(DelStrikes,'Delete Pressure(s)') == 1 % if pressures are deleted
%             Trial = str2double(cell2mat(inputdlg(['How many presures would you like to delete in Trial ?' num2str(PoorPress(i)) '?'])));
%             for j = 1:Trial
%                 uiwait(msgbox('Using your cursor, draw a box around the step you would like to delete.', 'Instructions'));
%                 subplot(1,length(PoorPress),i);
%                 h = imrect;
%                 pos(j).Coord = getPosition(h);
%             end
%             for j = 1:Trial % loop through # to delete
%                 for  k = 1:length(Regions{1,PoorPress(i)}) % loop through valid foot stikes to find the one to delete
%                     if Regions{1,PoorPress(i)}(k).Centroid(1) > pos(j).Coord(1) &&  Regions{1,PoorPress(i)}(k).Centroid(1) < pos(j).Coord(1) + pos(j).Coord(3)
%                         if Regions{1,PoorPress(i)}(k).Centroid(2) > pos(j).Coord(2) &&  Regions{1,PoorPress(i)}(k).Centroid(2) < pos(j).Coord(2) + pos(j).Coord(4)
%                             Press2Delete(j) = k;
%                         end
%                     end
%                 end
%             end
%             
%             % delete poor pressure(s) from structures
%             for j = 1:length(Press2Delete)
%                 if Press2Delete(j) < 3
%                     if strcmp(Regions{1,PoorPress(i)}(Press2Delete(j)).Side, 'Left')
%                         DynamicPPTrials(PoorPress(i)).PoorPressure.Left(1) = [];
%                     else
%                         DynamicPPTrials(PoorPress(i)).PoorPressure.Right(1) = [];
%                     end
%                 elseif Press2Delete(j) > 2 && Press2Delete < 5
%                     if strcmp(Regions{1,PoorPress(i)}(Press2Delete(j)).Side, 'Left')
%                         DynamicPPTrials(PoorPress(i)).PoorPressure.Left(2) = [];
%                     else
%                         DynamicPPTrials(PoorPress(i)).PoorPressure.Right(2) = [];
%                     end
%                 else
%                     if strcmp(Regions{1,PoorPress(i)}(Press2Delete(j)).Side, 'Left')
%                         DynamicPPTrials(PoorPress(i)).PoorPressure.Left(3) = [];
%                     else
%                         DynamicPPTrials(PoorPress(i)).PoorPressure.Right(3) = [];
%                     end
%                 end
%             end
%             Regions{1,PoorPress(i)}(Press2Delete) = [];
%             
%         elseif strcmp(DelStrikes, 'Keep Pressure(s)') % if pressures are kept
%             Trial = sum([Regions{1,PoorPress(i)}.PoorPressure]);
%             Step = 1;
%             for j = 1:length(Regions{1,PoorPress(i)})
%                 
%                 if Regions{1,PoorPress(i)}(j).PoorPressure == 0
%                     continue
%                 end
%                 uiwait(msgbox(['Please manually draw a line of progression for the ' num2str(Trial) ' invalid pressure(s) in Trial ' num2str(PoorPress(i)) '. Draw this line as if the whole foot was touching the ground.'], 'Instructions'));
%                 
%                 if j > 2
%                     Step = 2;
%                 elseif j > 4
%                     Step = 3;
%                 end
%                 
%                 Locs = ginput(2); % choose locations on figure
%                 line([Locs(1,1), Locs(2,1)], [Locs(1,2),Locs(2,2)],'Color', 'k'); % draw line
%                 
%                 for k = 1:length(Regions{1,PoorPress(i)}) % loop through valid foot stikes to find the one reassign
%                     if Regions{1,PoorPress(i)}(k).Centroid(1) > min(Locs(:,1))-5 &&  Regions{1,PoorPress(i)}(k).Centroid(1) < max(Locs(:,1))+5
%                         if Regions{1,PoorPress(i)}(k).Centroid(2) > min(Locs(:,2))-10 &&  Regions{1,PoorPress(i)}(k).Centroid(2) < max(Locs(:,2))+10
%                             Regions{1,PoorPress(i)}(k).PoorPressHT = Locs;
%                             Regions{1,PoorPress(i)}(k).Centroid = [mean(Locs(:,1)), mean(Locs(:,2))];
%                             if strcmp(Regions{1,PoorPress(i)}(j).Side, 'Left')
%                                 Regions{1,PoorPress(i)}(k).Orientation = 90 - (-180/pi*atan2(Locs(2,1) - Locs(1,1),Locs(2,2) - Locs(1,2)));
%                                 DynamicPPTrials(PoorPress(i)).Lht.Manual{Step} = Locs;
%                             else
%                                 Regions{1,PoorPress(i)}(k).Orientation = 90 - (180/pi*atan2(Locs(2,1) - Locs(1,1),Locs(2,2) - Locs(1,2)));
%                                 DynamicPPTrials(PoorPress(i)).Rht.Manual{Step} = Locs;
%                             end
%                         end
%                     end
%                 end
%                 clearvars Locs
%             end
%             
%         elseif strcmp(DelStrikes,'Delete Whole Trial') == 1 % delete whole trial
%             ToDel(w) = PoorPress(i);
%             w = w +1; 
%         end
%     end
%     close;
% end
% % if w>1
% %             Regions{PoorPress(i)} = [];
% % DynamicPPTrials
% %             Selection(PoorPress(i)) = [];
% %             NumDynTrials = length(Selection); 
% % end
% % clearvars NumStrikeEst pos h i j k button question Question q Press2Delete Trial AddDelStrikes Proceed
% 
% 
% %% plot a series of N trials for visualization.
% clc; disp('Displaying trials and initial masking');
% clearvars i j p Ti ToDel
% % determine screen size in pixels
% % Sets the units of your root object (screen) to pixels
% set(0,'units','pixels');
% % Obtains this pixel information
Pix_SS = get(0,'screensize');
Pix_SS(1:2) = [];
DisplayDim = [40, 40, NumDynTrials*Pix_SS(1)*0.1, Pix_SS(2)*.8];
% 
% % plot all the dynamic trials
% Cont = figure;
% set(Cont, 'Position', DisplayDim);
% for i = 1:NumDynTrials
%     subplot(1,NumDynTrials, i);
%     contour(DynamicPPTrials(i).SumTM,100);
%     axis equal;
%     title(['Trial ' num2str(i)], 'FontSize',8);
%     for j = 1:length(Regions{1,i})
%           % add boxes around feet for R/L classification
%         if  strcmp(Regions{1,i}(j).Side,'Left') == 1 % if left foot strike
%             rectangle('Position', Regions{1,i}(j).BoundingBox, 'EdgeColor','r'); % draw red box
%         else % if right foot strike
%             rectangle('Position', Regions{1,i}(j).BoundingBox, 'EdgeColor','g'); % draw green box
%         end
%         
%         % automatically generate foot progression angle using general image processing line
%         ydist = Subject.FootLength ./ Adj.Length .* sind(Regions{1,i}(j).Orientation);
%         xdist = Subject.FootLength ./ Adj.Length .* cosd(Regions{1,i}(j).Orientation);
%         Pt1 = [Regions{1,i}(j).Centroid(1) + (xdist/2), Regions{1,i}(j).Centroid(2) - (ydist/2)];
%         Pt2 = [Regions{1,i}(j).Centroid(1) - (xdist/2), Regions{1,i}(j).Centroid(2) + (ydist/2)];
%         Regions{1,i}(j).HT.General = [Pt1(1) Pt2(1); Pt1(2) Pt2(2)];
%         Regions{1,i}(j).FP_LineGeneral = line([Pt1(1) Pt2(1)],[Pt1(2) Pt2(2)] , 'Color', rgb('Gold'));
%         
%         % if its a poor pressure that was manually labeled, automatically generate foot progression angle
%         if Regions{1,i}(j).PoorPressure == 1
%             Regions{1,i}(j).HT.Manual = Regions{1,i}(j).PoorPressHT;
%             Regions{1,i}(j).FP_LineManual = line(Regions{1,i}(j).PoorPressHT(:,1),Regions{1,i}(j).PoorPressHT(:,2) , 'Color', rgb('HotPink'));
%         end
%         
%         % automatically generate foot progression angle using CoP trajectory line
%         ydist = Regions{1,i}(j).MajorAxisLength .* cosd(Regions{1,i}(j).CoP.FootProg);
%         if  strcmp(Regions{1,i}(j).Side,'Left') == 1 % if left foot strike
%             xdist = Regions{1,i}(j).MajorAxisLength .* sind(Regions{1,i}(j).CoP.FootProg);
%         else
%             xdist = Regions{1,i}(j).MajorAxisLength .* -sind(Regions{1,i}(j).CoP.FootProg);
%         end
%         Pt1 = [Regions{1,i}(j).CoP.Cent(1) + (xdist/2), Regions{1,i}(j).CoP.Cent(2) - (ydist/2)];
%         Pt2 = [Regions{1,i}(j).CoP.Cent(1) - (xdist/2), Regions{1,i}(j).CoP.Cent(2) + (ydist/2)];
%         Regions{1,i}(j).HT.CoP = [Pt1(1) Pt2(1); Pt1(2) Pt2(2)];
%         Regions{1,i}(j).FP_LineCoP = line([Pt1(1) Pt2(1)],[Pt1(2) Pt2(2)] , 'Color', rgb('Purple'));
%         % plot CoP trajectory
%         hold on;
%         for k = 1:length(Regions{1,i}(j).StepCoP)
%             x = Regions{1,i}(j).BoundingBox(1) + Regions{1,i}(j).StepCoP(k,1);
%             y = Regions{1,i}(j).BoundingBox(2) + Regions{1,i}(j).StepCoP(k,2);
%             plot(x,y,'.k','MarkerSize', 8);
%         end
%         
%         % 66% CoP trajectory line
%         ydist = Regions{1,i}(j).MajorAxisLength .* cosd(Regions{1,i}(j).CoP66.FootProg);
%         if  strcmp(Regions{1,i}(j).Side,'Left') == 1 % if left foot strike
%             xdist = Regions{1,i}(j).MajorAxisLength .* sind(Regions{1,i}(j).CoP66.FootProg);
%         else
%             xdist = Regions{1,i}(j).MajorAxisLength .* -sind(Regions{1,i}(j).CoP66.FootProg);
%         end
%         Pt1 = [Regions{1,i}(j).CoP66.Cent(1) + (xdist/2), Regions{1,i}(j).CoP66.Cent(2) - (ydist/2)];
%         Pt2 = [Regions{1,i}(j).CoP66.Cent(1) - (xdist/2), Regions{1,i}(j).CoP66.Cent(2) + (ydist/2)];
%         hold on;
%         Regions{1,i}(j).HT.CoP66 = [Pt1(1) Pt2(1); Pt1(2) Pt2(2)];
%         Regions{1,i}(j).FP_Line66 = line([Pt1(1) Pt2(1)],[Pt1(2) Pt2(2)] , 'Color', rgb('Cyan'));
%         % plot 66% CoP trajectory
%         %         hold on;
%         %         for k = 1:length(Regions{1,i}(j).StepCoP66)
%         %             x = Regions{1,i}(j).BoundingBox(1) + Regions{1,i}(j).StepCoP66(k,1);
%         %             y = Regions{1,i}(j).BoundingBox(2) + Regions{1,i}(j).StepCoP66(k,2);
%         %             plot(x,y,'om','MarkerSize', 6);
%         %         end
%         
%         % InterPeak CoP trajectory line
%         if Block.IP == 0
%             ydist = Regions{1,i}(j).MajorAxisLength .* cosd(Regions{1,i}(j).IP.FootProg);
%             if  strcmp(Regions{1,i}(j).Side,'Left') == 1 % if left foot strike
%                 xdist = Regions{1,i}(j).MajorAxisLength .* sind(Regions{1,i}(j).IP.FootProg);
%             else
%                 xdist = Regions{1,i}(j).MajorAxisLength .* -sind(Regions{1,i}(j).IP.FootProg);
%             end
%             Pt1 = [Regions{1,i}(j).IP.Cent(1) + (xdist/2), Regions{1,i}(j).IP.Cent(2) - (ydist/2)];
%             Pt2 = [Regions{1,i}(j).IP.Cent(1) - (xdist/2), Regions{1,i}(j).IP.Cent(2) + (ydist/2)];
%             hold on;
%             Regions{1,i}(j).HT.IP = [Pt1(1) Pt2(1); Pt1(2) Pt2(2)];
%             Regions{1,i}(j).FP_LineIP = line([Pt1(1) Pt2(1)],[Pt1(2) Pt2(2)] , 'Color', rgb('ForestGreen'));
%         end
%         
%         % heel peak point to centroid line
%         if Block.HC == 0
%             ydist = Regions{1,i}(j).MajorAxisLength .* cosd(Regions{1,i}(j).HC.FootProg);
%             if  strcmp(Regions{1,i}(j).Side,'Left') == 1 % if left foot strike
%                 xdist = Regions{1,i}(j).MajorAxisLength .* sind(Regions{1,i}(j).HC.FootProg);
%             else
%                 xdist = Regions{1,i}(j).MajorAxisLength .* -sind(Regions{1,i}(j).HC.FootProg);
%             end
%             Pt1 = [Regions{1,i}(j).Centroid(1) + (xdist/2), Regions{1,i}(j).Centroid(2) - (ydist/2)];
%             Pt2 = [Regions{1,i}(j).Centroid(1) - (xdist/2), Regions{1,i}(j).Centroid(2) + (ydist/2)];
%             hold on;
%             Regions{1,i}(j).HT.HC = [Pt1(1) Pt2(1); Pt1(2) Pt2(2)];
%             Regions{1,i}(j).FP_LineHeelCent = line([Pt1(1) Pt2(1)],[Pt1(2) Pt2(2)] , 'Color', rgb('Red'));
%         end
%         
%         if strcmp(PPSettings.C3DInput, 'Yes')
%             % plot c3d data heel/arch and forefoot arch points
%             if strcmp(PPSettings.C3DInput,'Yes')
%                 if isempty(C3Ddata(i).LHA) == 0
% %                     plot(C3Ddata(i).LHA(:,2), C3Ddata(i).LHA(:,1), 'ok');
% %                     plot(C3Ddata(i).LFA(:,2), C3Ddata(i).LFA(:,1), 'ok');
%                     plot(C3Ddata(i).LHEE.base(:,2), C3Ddata(i).LHEE.base(:,1), '*k');
%                     plot(C3Ddata(i).LTOE.top(:,2), C3Ddata(i).LTOE.top(:,1), 'ok');
%                     
%                     plot(C3Ddata(i).LHAb(:,2), C3Ddata(i).LHAb(:,1), 'ok');
%                     plot(C3Ddata(i).LFAb(:,2), C3Ddata(i).LFAb(:,1), 'ok');
%                     
%                     plot(C3Ddata(i).LLCA.step(:,2), C3Ddata(i).LLCA.step(:,1), '*k');
%                     plot(C3Ddata(i).LMCA.step(:,2), C3Ddata(i).LMCA.step(:,1), '*k');
%                     plot(C3Ddata(i).LD1M.step(:,2), C3Ddata(i).LD1M.step(:,1), '*k');
%                     plot(C3Ddata(i).LD5M.step(:,2), C3Ddata(i).LD5M.step(:,1), '*k');
%                     
%                 end
%                 if isempty(C3Ddata(i).RHA) == 0
% %                     plot(C3Ddata(i).RHA(:,2), C3Ddata(i).RHA(:,1), 'ok');
% %                     plot(C3Ddata(i).RFA(:,2), C3Ddata(i).RFA(:,1), 'ok');
%                     plot(C3Ddata(i).RHEE.base(:,2), C3Ddata(i).RHEE.base(:,1), '*k');
%                     plot(C3Ddata(i).RTOE.top(:,2), C3Ddata(i).RTOE.top(:,1), 'ok');
%                     
%                     plot(C3Ddata(i).RHAb(:,2), C3Ddata(i).RHAb(:,1), 'ok');
%                     plot(C3Ddata(i).RFAb(:,2), C3Ddata(i).RFAb(:,1), 'ok');
%                     
%                     plot(C3Ddata(i).RLCA.step(:,2), C3Ddata(i).RLCA.step(:,1), '*k');
%                     plot(C3Ddata(i).RMCA.step(:,2), C3Ddata(i).RMCA.step(:,1), '*k');
%                     plot(C3Ddata(i).RD1M.step(:,2), C3Ddata(i).RD1M.step(:,1), '*k');
%                     plot(C3Ddata(i).RD5M.step(:,2), C3Ddata(i).RD5M.step(:,1), '*k');
%                 end
%             end
%         end
%     end
%     clearvars pCoP x y
%     set(gca, 'XTick', []);
%     set(gca, 'YTick', []);
% end
% clearvars ToDel xdist ydist Pt1 Pt2 Ti
% 
% %% Double check that all trails and pressures look correct
% if NumDynTrials == 1
%     ListStr = {'1'};
% elseif NumDynTrials == 2
%     ListStr = {'1','2'};
% elseif NumDynTrials == 3
%     ListStr = {'1','2','3'};
% elseif NumDynTrials == 4
%     ListStr = {'1','2','3','4'};
% elseif NumDynTrials == 5
%     ListStr = {'1','2','3','4','5'};
% elseif NumDynTrials == 6
%     ListStr = {'1','2','3','4','5','6'};
% elseif NumDynTrials == 7
%     ListStr = {'1','2','3','4','5','6','7'};
% elseif NumDynTrials == 8
%     ListStr = {'1','2','3','4','5','6','7','8'};
% elseif NumDynTrials == 9
%     ListStr = {'1','2','3','4','5','6','7','8','9'};
% elseif NumDynTrials == 10
%     ListStr = {'1','2','3','4','5','6','7','8','9','10'};
% end
% clearvars ToDel
% DblCheck = questdlg('Would you like to delete any of the pressure/trial(s)?');
% if strcmp(DblCheck,'Yes') == 1
%     DblEdit = listdlg( 'PromptString','Select Trials to edit','ListString',ListStr); % determine which trials to re-edit
%     % mention or determine how many pressures are detected
%     z = 1;
%     for i = 1:length(DblEdit)
%         % choose to add pressure, delete pressure, or delete trial
%         AddDelStrikes = questdlg(['For Trial ' num2str(DblEdit(i)) ', would you like to delete a pressure or delete the whole trial?'], 'Add/Delete/Cancel','Delete Pressure','Delete Trial','Add');
%         if strcmp(AddDelStrikes,'Delete Pressure') == 1
%             Num2Del = str2double(cell2mat(inputdlg(['How many presures would you like to delete in Trial ' num2str(DblEdit(i)) '?'])));
%             for j = 1:Num2Del
%                 uiwait(msgbox('Using your cursor, draw a box around the foot you would like to delete.', 'Instructions'));
%                 subplot(1,NumDynTrials,DblEdit(i));
%                 h = imrect;
%                 pos(j).Coord = getPosition(h);
%             end
%             for j = 1:Num2Del % loop through # to delete
%                 for  k = 1:length(Regions{1,DblEdit(i)}) % loop through valid foot stikes to find the one to delete
%                     if Regions{1,DblEdit(i)}(k).Centroid(1) > pos(j).Coord(1) &&  Regions{1,DblEdit(i)}(k).Centroid(1) < pos(j).Coord(1) + pos(j).Coord(3)
%                         if Regions{1,DblEdit(i)}(k).Centroid(2) > pos(j).Coord(2) &&  Regions{1,DblEdit(i)}(k).Centroid(2) < pos(j).Coord(2) + pos(j).Coord(4)
%                             Press2Delete(j) = k;
%                         end
%                     end
%                 end
%             end
%             Regions{1,DblEdit(i)}(Press2Delete) = [];
%         else % delete whole trial
%             Regions{DblEdit(i)} = [];
%             ToDel(z) = DblEdit(i);
%             z = z + 1;
%         end
%     end
%     if exist('ToDel','var')
%         Regions{[ToDel]} = [];
%         DynamicPPTrials([ToDel]) = [];
% %         Selection(ToDel) = [];
%         NumDynTrials = length(DynamicPPTrials);
%     end
% end
% 
% %% Select Trials for further analysis
% if NumDynTrials == 1
%     Selection = 1;
% end
% 
% if strcmp(PPSettings.AutoSelectAll,'Yes')
%     Selection = 1:NumDynTrials;
% else
%     % determine all trials for final selection
%     uiwait(msgbox(['Select trials for further analysis that have footstrikes cleanly on the pressure plate and quality foot pressures.'], 'Instructions'));
%     Selection = listdlg( 'PromptString','Select Trials to further analyze','ListString',ListStr);
% end
% 
% %% Determine Masking Type if "select after display" chosen
% if strcmp(PPSettings.MaskChoice, 'Select after Display')
%     MaskList = {'General (yellow)','CoP (purple)','66% CoP (cyan)','Inter-Peak (green)','Heel-Centroid (red)','Manual (pink)'};
%     A = listdlg( 'PromptString','Choose a masking method','ListString',MaskList, 'SelectionMode','single');
%     PPSettings.MaskChoice = MaskList{A};
% end
% 
% %% Determine numbers of steps on each side for each trial
% for i = 1:length(Selection)
%     % LEFT
%     k = 1;
%     for j = 1:length(Regions{1,Selection(i)})
%         FindLeft(j) = strcmp(Regions{1,Selection(i)}(j).Side, 'Left');
%         if FindLeft(j) == 1
%             DynamicPPTrials(Selection(i)).Lht.General{k} = [Regions{1,Selection(i)}(j).FP_LineGeneral.XData' Regions{1,Selection(i)}(j).FP_LineGeneral.YData'];
%             DynamicPPTrials(Selection(i)).Lht.CoP{k} = [Regions{1,Selection(i)}(j).FP_LineCoP.XData' Regions{1,Selection(i)}(j).FP_LineCoP.YData'];
%             DynamicPPTrials(Selection(i)).Lht.CoP66{k} = [Regions{1,Selection(i)}(j).FP_Line66.XData' Regions{1,Selection(i)}(j).FP_Line66.YData'];
%             if Block.IP == 0
%                 DynamicPPTrials(Selection(i)).Lht.IP{k} = [Regions{1,Selection(i)}(j).FP_LineIP.XData' Regions{1,Selection(i)}(j).FP_LineIP.YData'];
%             end
%             if Block.HC == 0
%                 DynamicPPTrials(Selection(i)).Lht.HeelCent{k} = [Regions{1,Selection(i)}(j).FP_LineHeelCent.XData' Regions{1,Selection(i)}(j).FP_LineHeelCent.YData'];
%             end
%             if Regions{1,Selection(i)}(j).PoorPressure == 1
%                 DynamicPPTrials(Selection(i)).Lht.Manual{k} = [Regions{1,Selection(i)}(j).FP_LineManual.XData' Regions{1,Selection(i)}(j).FP_LineManual.YData'];
%             end
%             k = k+1;
%         end
%     end
%     DynamicPPTrials(Selection(i)).NumLeft = sum(FindLeft);
%     
%     % RIGHT
%     k = 1;
%     for j = 1:length(Regions{1,Selection(i)})
%         FindRight(j) = strcmp(Regions{1,Selection(i)}(j).Side, 'Right');
%         if FindRight(j) == 1
%             DynamicPPTrials(Selection(i)).Rht.General{k} = [Regions{1,Selection(i)}(j).FP_LineGeneral.XData' Regions{1,Selection(i)}(j).FP_LineGeneral.YData'];
%             DynamicPPTrials(Selection(i)).Rht.CoP{k} = [Regions{1,Selection(i)}(j).FP_LineCoP.XData' Regions{1,Selection(i)}(j).FP_LineCoP.YData'];
%             DynamicPPTrials(Selection(i)).Rht.CoP66{k} = [Regions{1,Selection(i)}(j).FP_Line66.XData' Regions{1,Selection(i)}(j).FP_Line66.YData'];
%             if Block.IP == 0
%                 DynamicPPTrials(Selection(i)).Rht.IP{k} = [Regions{1,Selection(i)}(j).FP_LineIP.XData' Regions{1,Selection(i)}(j).FP_LineIP.YData'];
%             end
%             if Block.HC == 0
%                 DynamicPPTrials(Selection(i)).Rht.HeelCent{k} = [Regions{1,Selection(i)}(j).FP_LineHeelCent.XData' Regions{1,Selection(i)}(j).FP_LineHeelCent.YData'];
%             end
%             if Regions{1,Selection(i)}(j).PoorPressure == 1
%                 DynamicPPTrials(Selection(i)).Rht.Manual{k} = [Regions{1,Selection(i)}(j).FP_LineManual.XData' Regions{1,Selection(i)}(j).FP_LineManual.YData'];
%             end
%             k = k+1;
%         end
%     end
%     DynamicPPTrials(Selection(i)).NumRight = sum(FindRight);
%     clearvars FindRight FindLeft
% end
% 
% % redefine number of steps on each side after deletions
% % for i = 1:length(Selection)
% %     % LEFT
% %     k = 1;
% %     for j = 1:length(Regions{1,Selection(i)})
% %         FindLeft(j) = strcmp(Regions{1,Selection(i)}(j).Side, 'Left');
% %     end
% %     DynamicPPTrials(Selection(i)).NumLeft = sum(FindLeft);
% %     
% %     % RIGHT
% %     k = 1;
% %     for j = 1:length(Regions{1,Selection(i)})
% %         FindRight(j) = strcmp(Regions{1,Selection(i)}(j).Side, 'Right');
% %     end
% %     DynamicPPTrials(Selection(i)).NumRight = sum(FindRight);
% %     clearvars FindRight FindLeft
% % end
% 
% %% Manual Mask
% if strcmp(PPSettings.MaskChoice, 'Manual') || strcmp(PPSettings.MaskChoice, 'Manual (pink)') || strcmp(PPSettings.MaskChoice, 'Validation')
%     clc; disp('Manual masking');
%     Trials2EditFP = listdlg('PromptString','Select trial(s) to manually edit foot prog angle','ListString', ListStr(Selection));
%     close;
%     M = figure('Position', [DisplayDim(1) DisplayDim(2) (DisplayDim(3)/NumDynTrials)+100 DisplayDim(4)]);
%     % edit foot progression angles
%     for i = 1:length(Trials2EditFP)
%         contour(DynamicPPTrials(Trials2EditFP(i)).SumTM,100);
%         axis equal;
%         title(['Trial ' num2str(i)], 'FontSize',8);
%         ax = gca;
%         ax.XTick = [];
%         ax.YTick = [];
%         % plot FP angles from poor pressures
%         if DynamicPPTrials(Trials2EditFP(i)).NumLeft > 0
%             for j = 1:DynamicPPTrials(Trials2EditFP(i)).NumLeft
%                 if DynamicPPTrials(Trials2EditFP(i)).PoorPressure.Left(j) == 1
%                     line(DynamicPPTrials(Trials2EditFP(i)).Lht.Manual{(j)}(:,1),DynamicPPTrials(Trials2EditFP(i)).Lht.Manual{(j)}(:,2),'Color','k','LineWidth',1.5);
%                 end
%             end
%         end
%         if DynamicPPTrials(Trials2EditFP(i)).NumRight > 0
%             for j = 1:DynamicPPTrials(Trials2EditFP(i)).NumRight
%                 if DynamicPPTrials(Trials2EditFP(i)).PoorPressure.Right(j) == 1
%                     line(DynamicPPTrials(Trials2EditFP(i)).Rht.Manual{(j)}(:,1),DynamicPPTrials(Trials2EditFP(i)).Rht.Manual{(j)}(:,2),'Color','k','LineWidth',1.5);
%                 end
%             end
%         end
%         % manually define foot progression angles
%         % LEFT side
%         if DynamicPPTrials(Trials2EditFP(i)).NumLeft > 0
%             for j = 1:DynamicPPTrials(Trials2EditFP(i)).NumLeft
%                 if DynamicPPTrials(Trials2EditFP(i)).PoorPressure.Left(j) ~= 1
%                     uiwait(msgbox({['On trial ' num2str(Trials2EditFP(i)) ' use the cursor to select points just posterior to the heel and just anterior to the 2nd toe on LEFT step # ',num2str(j)]}));
%                     DynamicPPTrials(Trials2EditFP(i)).Lht.Manual{j} = ginput(2);
%                     line(DynamicPPTrials(Trials2EditFP(i)).Lht.Manual{(j)}(:,1),DynamicPPTrials(Trials2EditFP(i)).Lht.Manual{(j)}(:,2),'Color','k','LineWidth',1.5);
%                     question = 'Are you satisfied wtih the new foot progression angle?';
%                     ReEditFPAngles = questdlg(question,'ReDo Foot Progression Angles?', 'Yes','No','Yes');
%                     if strcmp(ReEditFPAngles, 'No') == 1
%                         uiwait(msgbox({['On trial ' num2str(Trials2EditFP(i)) ' use the cursor to select points just posterior to the heel and just anterior to the 2nd toe of the LEFT step # ',num2str(j)]}));
%                         DynamicPPTrials(Trials2EditFP(i)).Lht.Manual{j} = ginput(2);
%                         line(DynamicPPTrials(Trials2EditFP(i)).Lht.Manual{j}(:,1),DynamicPPTrials(Trials2EditFP(i)).Lht.Manual{j}(:,2),'Color','k','LineWidth',1.5);
%                     end
%                 end
%             end
%         end
%         % RIGHT side
%         if DynamicPPTrials(Trials2EditFP(i)).NumRight > 0
%             for j = 1:DynamicPPTrials(Trials2EditFP(i)).NumRight
%                 if DynamicPPTrials(Trials2EditFP(i)).PoorPressure.Right(j) ~= 1
%                     uiwait(msgbox({['On trial ' num2str(Trials2EditFP(i)) ' use the cursor to select points just posterior to the heel and just anterior to the 2nd toe of the RIGHT step # ',num2str(j)]}));
%                     DynamicPPTrials(Trials2EditFP(i)).Rht.Manual{j} = ginput(2);
%                     line(DynamicPPTrials(Trials2EditFP(i)).Rht.Manual{j}(:,1),DynamicPPTrials(Trials2EditFP(i)).Rht.Manual{j}(:,2),'Color','k','LineWidth',1.5);
%                     question = 'Are you satisfied wtih the new foot progression angle?';
%                     ReEditFPAngles = questdlg(question,'ReDo Foot Progression Angles?', 'Yes','No','Yes');
%                     if strcmp(ReEditFPAngles,'No') == 1
%                         uiwait(msgbox({['On trial ' num2str(Trials2EditFP(i)) ' use the cursor to select points just posterior to the heel and just anterior to the 2nd toe of the RIGHT step # ',num2str(j)]}));
%                         DynamicPPTrials(Trials2EditFP(i)).Rht.Manual{j} = ginput(2);
%                         line(DynamicPPTrials(Trials2EditFP(i)).Rht.Manual{j}(:,1),DynamicPPTrials(Trials2EditFP(i)).Rht.Manual{j}(:,2),'Color','k','LineWidth',1.5);
%                     end
%                 end
%             end
%         end
%         clf;
%     end
% end
% close;
% 
% %% ensure toe and heel points are oriented correctly
% clc; disp('Checking coordinates and calculating progression angles');
% % the top row should have a lower Y value than the bottom row, aka define the heel point before the toe point
% for i = 1:NumDynTrials % loop for length of trials
%     % LEFT
%     for j = 1:DynamicPPTrials(i).NumLeft % loop for # of Left feet
%         % manual points
%         if isempty(DynamicPPTrials(i).Lht.Manual{j}) == 0
%             if DynamicPPTrials(i).Lht.General{j}(1,2) > DynamicPPTrials(i).Lht.Manual{j}(2,2) % if top y value is greater than bottom y value -> swap
%                 Switch = DynamicPPTrials(i).Lht.Manual{j}(1,:);
%                 DynamicPPTrials(i).Lht.Manual{j}(1,:) = DynamicPPTrials(i).Lht.Manual{j}(2,:);
%                 DynamicPPTrials(i).Lht.Manual{j}(2,:) = Switch;
%             end
%         end
%         clearvars Switch
%         % original image processing FP points
%         if isempty(DynamicPPTrials(i).Lht.General) == 0
%             if DynamicPPTrials(i).Lht.General{j}(1,2) > DynamicPPTrials(i).Lht.General{j}(2,2) % if top y value is greater than bottom y value -> swap
%                 Switch = DynamicPPTrials(i).Lht.General{j}(1,:);
%                 DynamicPPTrials(i).Lht.General{j}(1,:) = DynamicPPTrials(i).Lht.General{j}(2,:);
%                 DynamicPPTrials(i).Lht.General{j}(2,:) = Switch;
%             end
%         end
%         clearvars Switch
%         % CoP FP points
%         if isempty(DynamicPPTrials(i).Lht.CoP) == 0
%             if DynamicPPTrials(i).Lht.CoP{j}(1,2) > DynamicPPTrials(i).Lht.CoP{j}(2,2) % if top y value is greater than bottom y value -> swap
%                 Switch = DynamicPPTrials(i).Lht.CoP{j}(1,:);
%                 DynamicPPTrials(i).Lht.CoP{j}(1,:) = DynamicPPTrials(i).Lht.CoP{j}(2,:);
%                 DynamicPPTrials(i).Lht.CoP{j}(2,:) = Switch;
%             end
%         end
%         clearvars Switch
%         % 66% CoP FP Points
%         if isempty(DynamicPPTrials(i).Lht.CoP66) == 0
%             if DynamicPPTrials(i).Lht.CoP66{j}(1,2) > DynamicPPTrials(i).Lht.CoP66{j}(2,2) % if top y value is greater than bottom y value -> swap
%                 Switch = DynamicPPTrials(i).Lht.CoP66{j}(1,:);
%                 DynamicPPTrials(i).Lht.CoP66{j}(1,:) = DynamicPPTrials(i).Lht.CoP66{j}(2,:);
%                 DynamicPPTrials(i).Lht.CoP66{j}(2,:) = Switch;
%             end
%         end
%         clearvars Switch
%         % inter-peak CoP Points
%         if Block.IP == 0
%             if isempty(DynamicPPTrials(i).Lht.IP) == 0
%                 if DynamicPPTrials(i).Lht.IP{j}(1,2) > DynamicPPTrials(i).Lht.IP{j}(2,2) % if top y value is greater than bottom y value -> swap
%                     Switch = DynamicPPTrials(i).Lht.IP{j}(1,:);
%                     DynamicPPTrials(i).Lht.IP{j}(1,:) = DynamicPPTrials(i).Lht.IP{j}(2,:);
%                     DynamicPPTrials(i).Lht.IP{j}(2,:) = Switch;
%                 end
%             end
%             clearvars Switch
%         end
%         % Heel peak to centroid FP Points
%         if Block.HC == 0
%             if isempty(DynamicPPTrials(i).Lht.HeelCent) == 0
%                 if DynamicPPTrials(i).Lht.HeelCent{j}(1,2) > DynamicPPTrials(i).Lht.HeelCent{j}(2,2) % if top y value is greater than bottom y value -> swap
%                     Switch = DynamicPPTrials(i).Lht.HeelCent{j}(1,:);
%                     DynamicPPTrials(i).Lht.HeelCent{j}(1,:) = DynamicPPTrials(i).Lht.HeelCent{j}(2,:);
%                     DynamicPPTrials(i).Lht.HeelCent{j}(2,:) = Switch;
%                 end
%             end
%         end
%         clearvars Switch
%     end
%     
%     % RIGHT
%     clearvars Switch
%     for j = 1:DynamicPPTrials(i).NumRight% loop for # of Right feet
%         % manual points
%         if isempty(DynamicPPTrials(i).Rht.Manual{j}) == 0
%             if DynamicPPTrials(i).Rht.Manual{j}(1,2) > DynamicPPTrials(i).Rht.Manual{j}(2,2) % if top y value is greater than bottom y value -> swap
%                 Switch = DynamicPPTrials(i).Rht.Manual{j}(1,:);
%                 DynamicPPTrials(i).Rht.Manual{j}(1,:) = DynamicPPTrials(i).Rht.Manual{j}(2,:);
%                 DynamicPPTrials(i).Rht.Manual{j}(2,:) = Switch;
%             end
%         end
%         if isempty(DynamicPPTrials(i).Rht.General) == 0
%             if DynamicPPTrials(i).Rht.General{j}(1,2) > DynamicPPTrials(i).Rht.General{j}(2,2) % if top y value is greater than bottom y value -> swap
%                 Switch = DynamicPPTrials(i).Rht.General{j}(1,:);
%                 DynamicPPTrials(i).Rht.General{j}(1,:) = DynamicPPTrials(i).Rht.General{j}(2,:);
%                 DynamicPPTrials(i).Rht.General{j}(2,:) = Switch;
%             end
%         end
%         clearvars Switch
%         % CoP FP points
%         if isempty(DynamicPPTrials(i).Rht.CoP) == 0
%             if DynamicPPTrials(i).Rht.CoP{j}(1,2) > DynamicPPTrials(i).Rht.CoP{j}(2,2) % if top y value is greater than bottom y value -> swap
%                 Switch = DynamicPPTrials(i).Rht.CoP{j}(1,:);
%                 DynamicPPTrials(i).Rht.CoP{j}(1,:) = DynamicPPTrials(i).Rht.CoP{j}(2,:);
%                 DynamicPPTrials(i).Rht.CoP{j}(2,:) = Switch;
%             end
%         end
%         clearvars Switch
%         % 66% CoP FP Points
%         if isempty(DynamicPPTrials(i).Rht.CoP66) == 0
%             if DynamicPPTrials(i).Rht.CoP66{j}(1,2) > DynamicPPTrials(i).Rht.CoP66{j}(2,2) % if top y value is greater than bottom y value -> swap
%                 Switch = DynamicPPTrials(i).Rht.CoP66{j}(1,:);
%                 DynamicPPTrials(i).Rht.CoP66{j}(1,:) = DynamicPPTrials(i).Rht.CoP66{j}(2,:);
%                 DynamicPPTrials(i).Rht.CoP66{j}(2,:) = Switch;
%             end
%         end
%         clearvars Switch
%         % inter-peak CoP Points
%         if Block.IP == 0
%             if isempty(DynamicPPTrials(i).Rht.IP) == 0
%                 if DynamicPPTrials(i).Rht.IP{j}(1,2) > DynamicPPTrials(i).Rht.IP{j}(2,2) % if top y value is greater than bottom y value -> swap
%                     Switch = DynamicPPTrials(i).Rht.IP{j}(1,:);
%                     DynamicPPTrials(i).Rht.IP{j}(1,:) = DynamicPPTrials(i).Rht.IP{j}(2,:);
%                     DynamicPPTrials(i).Rht.IP{j}(2,:) = Switch;
%                 end
%             end
%             clearvars Switch
%         end
%         % Heel peak to centroid FP Points
%         if Block.HC == 0
%             if isempty(DynamicPPTrials(i).Rht.HeelCent) == 0
%                 if DynamicPPTrials(i).Rht.HeelCent{j}(1,2) > DynamicPPTrials(i).Rht.HeelCent{j}(2,2) % if top y value is greater than bottom y value -> swap
%                     Switch = DynamicPPTrials(i).Rht.HeelCent{j}(1,:);
%                     DynamicPPTrials(i).Rht.HeelCent{j}(1,:) = DynamicPPTrials(i).Rht.HeelCent{j}(2,:);
%                     DynamicPPTrials(i).Rht.HeelCent{j}(2,:) = Switch;
%                 end
%             end
%         end
%         clearvars Switch
%     end
% end
% 
% %% calculate prog angles
% for j = 1:length(Selection)
%     for i = 1:DynamicPPTrials(Selection(j)).NumLeft % LEFT side
%          % manual points
%          if strcmp(PPSettings.MaskChoice, 'Manual') || strcmp(PPSettings.MaskChoice, 'Manual (pink)') || strcmp(PPSettings.MaskChoice, 'Validation')
%              DynamicPPTrials(Selection(j)).LProg.AngManual{i} = round(abs(180/pi*(atan((DynamicPPTrials(Selection(j)).Lht.Manual{i}(2,1)-DynamicPPTrials(Selection(j)).Lht.Manual{i}(1,1))/(DynamicPPTrials(Selection(j)).Lht.Manual{i}(2,2)-DynamicPPTrials(Selection(j)).Lht.Manual{i}(1,2))))));
%              DynamicPPTrials(Selection(j)).LProg.AngActualManual{i} = -180/pi*(atan((DynamicPPTrials(Selection(j)).Lht.Manual{i}(2,1)-DynamicPPTrials(Selection(j)).Lht.Manual{i}(1,1)) / (DynamicPPTrials(Selection(j)).Lht.Manual{i}(2,2)-DynamicPPTrials(Selection(j)).Lht.Manual{i}(1,2))));
%              DynamicPPTrials(Selection(j)).Lht.ManualV{i} = [(DynamicPPTrials(Selection(j)).Lht.Manual{i}(2,1)-DynamicPPTrials(Selection(j)).Lht.Manual{i}(1,1)) (DynamicPPTrials(Selection(j)).Lht.Manual{i}(2,2)-DynamicPPTrials(Selection(j)).Lht.Manual{i}(1,2))]  / norm([(DynamicPPTrials(Selection(j)).Lht.Manual{i}(2,1) - DynamicPPTrials(Selection(j)).Lht.Manual{i}(1,1)) (DynamicPPTrials(Selection(j)).Lht.Manual{i}(2,2) - DynamicPPTrials(Selection(j)).Lht.Manual{i}(1,2))]);
%              DynamicPPTrials(Selection(j)).Lht.ManualInV{i} = [-DynamicPPTrials(Selection(j)).Lht.ManualV{i}(2) DynamicPPTrials(Selection(j)).Lht.ManualV{i}(1)]; % calculate inverted slope of prog angles - for drawing perpendicular segmentation lines
%          end
%         % general image recognition points
%         DynamicPPTrials(Selection(j)).LProg.AngGeneral{i} = round(abs(180/pi*(atan((DynamicPPTrials(Selection(j)).Lht.General{i}(2,1)-DynamicPPTrials(Selection(j)).Lht.General{i}(1,1))/(DynamicPPTrials(Selection(j)).Lht.General{i}(2,2)-DynamicPPTrials(Selection(j)).Lht.General{i}(1,2))))));
%         DynamicPPTrials(Selection(j)).LProg.AngActualGeneral{i} = -180/pi*(atan((DynamicPPTrials(Selection(j)).Lht.General{i}(2,1)-DynamicPPTrials(Selection(j)).Lht.General{i}(1,1)) / (DynamicPPTrials(Selection(j)).Lht.General{i}(2,2)-DynamicPPTrials(Selection(j)).Lht.General{i}(1,2))));
%         DynamicPPTrials(Selection(j)).Lht.GeneralV{i} = [(DynamicPPTrials(Selection(j)).Lht.General{i}(2,1)-DynamicPPTrials(Selection(j)).Lht.General{i}(1,1)) (DynamicPPTrials(Selection(j)).Lht.General{i}(2,2)-DynamicPPTrials(Selection(j)).Lht.General{i}(1,2))]  / norm([(DynamicPPTrials(Selection(j)).Lht.General{i}(2,1) - DynamicPPTrials(Selection(j)).Lht.General{i}(1,1)) (DynamicPPTrials(Selection(j)).Lht.General{i}(2,2) - DynamicPPTrials(Selection(j)).Lht.General{i}(1,2))]);
%         DynamicPPTrials(Selection(j)).Lht.GeneralInV{i} = [-DynamicPPTrials(Selection(j)).Lht.GeneralV{i}(2) DynamicPPTrials(Selection(j)).Lht.GeneralV{i}(1)]; % calculate inverted slope of prog angles - for drawing perpendicular segmentation lines
%         % whole center of pressure points
%         DynamicPPTrials(Selection(j)).LProg.AngCoP{i} = round(abs(180/pi*(atan((DynamicPPTrials(Selection(j)).Lht.CoP{i}(2,1)-DynamicPPTrials(Selection(j)).Lht.CoP{i}(1,1))/(DynamicPPTrials(Selection(j)).Lht.CoP{i}(2,2)-DynamicPPTrials(Selection(j)).Lht.CoP{i}(1,2))))));
%         DynamicPPTrials(Selection(j)).LProg.AngActualCoP{i} = -180/pi*(atan((DynamicPPTrials(Selection(j)).Lht.CoP{i}(2,1)-DynamicPPTrials(Selection(j)).Lht.CoP{i}(1,1)) / (DynamicPPTrials(Selection(j)).Lht.CoP{i}(2,2)-DynamicPPTrials(Selection(j)).Lht.CoP{i}(1,2))));
%         DynamicPPTrials(Selection(j)).Lht.VCoP{i} = [(DynamicPPTrials(Selection(j)).Lht.CoP{i}(2,1)-DynamicPPTrials(Selection(j)).Lht.CoP{i}(1,1)) (DynamicPPTrials(Selection(j)).Lht.CoP{i}(2,2)-DynamicPPTrials(Selection(j)).Lht.CoP{i}(1,2))]  / norm([(DynamicPPTrials(Selection(j)).Lht.CoP{i}(2,1) - DynamicPPTrials(Selection(j)).Lht.CoP{i}(1,1)) (DynamicPPTrials(Selection(j)).Lht.CoP{i}(2,2) - DynamicPPTrials(Selection(j)).Lht.CoP{i}(1,2))]);
%         DynamicPPTrials(Selection(j)).Lht.InVCoP{i} = [-DynamicPPTrials(Selection(j)).Lht.VCoP{i}(2) DynamicPPTrials(Selection(j)).Lht.VCoP{i}(1)]; % calculate inverted slope of prog angles - for drawing perpendicular segmentation lines
%         % 66% CoP points
%         DynamicPPTrials(Selection(j)).LProg.Ang66{i} = round(abs(180/pi*(atan((DynamicPPTrials(Selection(j)).Lht.CoP66{i}(2,1)-DynamicPPTrials(Selection(j)).Lht.CoP66{i}(1,1))/(DynamicPPTrials(Selection(j)).Lht.CoP66{i}(2,2)-DynamicPPTrials(Selection(j)).Lht.CoP66{i}(1,2))))));
%         DynamicPPTrials(Selection(j)).LProg.AngActual66{i} = -180/pi*(atan((DynamicPPTrials(Selection(j)).Lht.CoP66{i}(2,1)-DynamicPPTrials(Selection(j)).Lht.CoP66{i}(1,1)) / (DynamicPPTrials(Selection(j)).Lht.CoP66{i}(2,2)-DynamicPPTrials(Selection(j)).Lht.CoP66{i}(1,2))));
%         DynamicPPTrials(Selection(j)).Lht.CoPV66{i} = [(DynamicPPTrials(Selection(j)).Lht.CoP66{i}(2,1)-DynamicPPTrials(Selection(j)).Lht.CoP66{i}(1,1)) (DynamicPPTrials(Selection(j)).Lht.CoP66{i}(2,2)-DynamicPPTrials(Selection(j)).Lht.CoP66{i}(1,2))]  / norm([(DynamicPPTrials(Selection(j)).Lht.CoP66{i}(2,1) - DynamicPPTrials(Selection(j)).Lht.CoP66{i}(1,1)) (DynamicPPTrials(Selection(j)).Lht.CoP66{i}(2,2) - DynamicPPTrials(Selection(j)).Lht.CoP66{i}(1,2))]);
%         DynamicPPTrials(Selection(j)).Lht.CoPInV66{i} = [-DynamicPPTrials(Selection(j)).Lht.CoPV66{i}(2) DynamicPPTrials(Selection(j)).Lht.CoPV66{i}(1)]; % calculate inverted slope of prog angles - for drawing perpendicular segmentation lines
%         % interpeak points
%         if Block.IP == 0
%             DynamicPPTrials(Selection(j)).LProg.AngIP{i} = round(abs(180/pi*(atan((DynamicPPTrials(Selection(j)).Lht.IP{i}(2,1)-DynamicPPTrials(Selection(j)).Lht.IP{i}(1,1))/(DynamicPPTrials(Selection(j)).Lht.IP{i}(2,2)-DynamicPPTrials(Selection(j)).Lht.IP{i}(1,2))))));
%             DynamicPPTrials(Selection(j)).LProg.AngActualIP{i} = -180/pi*(atan((DynamicPPTrials(Selection(j)).Lht.IP{i}(2,1)-DynamicPPTrials(Selection(j)).Lht.IP{i}(1,1)) / (DynamicPPTrials(Selection(j)).Lht.IP{i}(2,2)-DynamicPPTrials(Selection(j)).Lht.IP{i}(1,2))));
%             DynamicPPTrials(Selection(j)).Lht.VIP{i} = [(DynamicPPTrials(Selection(j)).Lht.IP{i}(2,1)-DynamicPPTrials(Selection(j)).Lht.IP{i}(1,1)) (DynamicPPTrials(Selection(j)).Lht.IP{i}(2,2)-DynamicPPTrials(Selection(j)).Lht.IP{i}(1,2))]  / norm([(DynamicPPTrials(Selection(j)).Lht.IP{i}(2,1) - DynamicPPTrials(Selection(j)).Lht.IP{i}(1,1)) (DynamicPPTrials(Selection(j)).Lht.IP{i}(2,2) - DynamicPPTrials(Selection(j)).Lht.IP{i}(1,2))]);
%             DynamicPPTrials(Selection(j)).Lht.InVIP{i} = [-DynamicPPTrials(Selection(j)).Lht.VIP{i}(2) DynamicPPTrials(Selection(j)).Lht.VIP{i}(1)]; % calculate inverted slope of prog angles - for drawing perpendicular segmentation lines
%         end
%         % Heel to centroid points
%         if Block.HC == 0
%             DynamicPPTrials(Selection(j)).LProg.AngHeelCent{i} = round(abs(180/pi*(atan((DynamicPPTrials(Selection(j)).Lht.HeelCent{i}(2,1)-DynamicPPTrials(Selection(j)).Lht.HeelCent{i}(1,1))/(DynamicPPTrials(Selection(j)).Lht.HeelCent{i}(2,2)-DynamicPPTrials(Selection(j)).Lht.HeelCent{i}(1,2))))));
%             DynamicPPTrials(Selection(j)).LProg.AngActualHeelCent{i} = -180/pi*(atan((DynamicPPTrials(Selection(j)).Lht.HeelCent{i}(2,1)-DynamicPPTrials(Selection(j)).Lht.HeelCent{i}(1,1)) / (DynamicPPTrials(Selection(j)).Lht.HeelCent{i}(2,2)-DynamicPPTrials(Selection(j)).Lht.HeelCent{i}(1,2))));
%             DynamicPPTrials(Selection(j)).Lht.VHeelCent{i} = [(DynamicPPTrials(Selection(j)).Lht.HeelCent{i}(2,1)-DynamicPPTrials(Selection(j)).Lht.HeelCent{i}(1,1)) (DynamicPPTrials(Selection(j)).Lht.HeelCent{i}(2,2)-DynamicPPTrials(Selection(j)).Lht.HeelCent{i}(1,2))]  / norm([(DynamicPPTrials(Selection(j)).Lht.HeelCent{i}(2,1) - DynamicPPTrials(Selection(j)).Lht.HeelCent{i}(1,1)) (DynamicPPTrials(Selection(j)).Lht.HeelCent{i}(2,2) - DynamicPPTrials(Selection(j)).Lht.HeelCent{i}(1,2))]);
%             DynamicPPTrials(Selection(j)).Lht.InVHeelCent{i} = [-DynamicPPTrials(Selection(j)).Lht.VHeelCent{i}(2) DynamicPPTrials(Selection(j)).Lht.VHeelCent{i}(1)]; % calculate inverted slope of prog angles - for drawing perpendicular segmentation lines
%         end
%     end
% end
% for j = 1:length(Selection)
%     for i = 1:DynamicPPTrials(Selection(j)).NumRight % RIGHT side
%         % manual points
%         if strcmp(PPSettings.MaskChoice, 'Manual') || strcmp(PPSettings.MaskChoice, 'Manual (pink)') || strcmp(PPSettings.MaskChoice, 'Validation')
%             DynamicPPTrials(Selection(j)).RProg.AngManual{i} = round(abs(180/pi*(atan((DynamicPPTrials(Selection(j)).Rht.Manual{i}(2,1)-DynamicPPTrials(Selection(j)).Rht.Manual{i}(1,1))/(DynamicPPTrials(Selection(j)).Rht.Manual{i}(2,2)-DynamicPPTrials(Selection(j)).Rht.Manual{i}(1,2))))));
%             DynamicPPTrials(Selection(j)).RProg.AngActualManual{i} = -180/pi*(atan((DynamicPPTrials(Selection(j)).Rht.Manual{i}(2,1)-DynamicPPTrials(Selection(j)).Rht.Manual{i}(1,1))/(DynamicPPTrials(Selection(j)).Rht.Manual{i}(2,2)-DynamicPPTrials(Selection(j)).Rht.Manual{i}(1,2))));
%             DynamicPPTrials(Selection(j)).Rht.ManualV{i} = [(DynamicPPTrials(Selection(j)).Rht.Manual{i}(2,1)-DynamicPPTrials(Selection(j)).Rht.Manual{i}(1,1)) (DynamicPPTrials(Selection(j)).Rht.Manual{i}(2,2)-DynamicPPTrials(Selection(j)).Rht.Manual{i}(1,2))]  / norm([(DynamicPPTrials(Selection(j)).Rht.Manual{i}(2,1) - DynamicPPTrials(Selection(j)).Rht.Manual{i}(1,1)) (DynamicPPTrials(Selection(j)).Rht.Manual{i}(2,2) - DynamicPPTrials(Selection(j)).Rht.Manual{i}(1,2))]);
%             DynamicPPTrials(Selection(j)).Rht.ManualInV{i} = [-DynamicPPTrials(Selection(j)).Rht.ManualV{i}(2) DynamicPPTrials(Selection(j)).Rht.ManualV{i}(1)];  % calculate inverted slope of prog angles - for drawing perpendicular segmentation lines
%         end
%         % general image recognition points
%         DynamicPPTrials(Selection(j)).RProg.AngGeneral{i} = round(abs(180/pi*(atan((DynamicPPTrials(Selection(j)).Rht.General{i}(2,1)-DynamicPPTrials(Selection(j)).Rht.General{i}(1,1))/(DynamicPPTrials(Selection(j)).Rht.General{i}(2,2)-DynamicPPTrials(Selection(j)).Rht.General{i}(1,2))))));
%         DynamicPPTrials(Selection(j)).RProg.AngActualGeneral{i} = -180/pi*(atan((DynamicPPTrials(Selection(j)).Rht.General{i}(2,1)-DynamicPPTrials(Selection(j)).Rht.General{i}(1,1))/(DynamicPPTrials(Selection(j)).Rht.General{i}(2,2)-DynamicPPTrials(Selection(j)).Rht.General{i}(1,2))));
%         DynamicPPTrials(Selection(j)).Rht.GeneralV{i} = [(DynamicPPTrials(Selection(j)).Rht.General{i}(2,1)-DynamicPPTrials(Selection(j)).Rht.General{i}(1,1)) (DynamicPPTrials(Selection(j)).Rht.General{i}(2,2)-DynamicPPTrials(Selection(j)).Rht.General{i}(1,2))]  / norm([(DynamicPPTrials(Selection(j)).Rht.General{i}(2,1) - DynamicPPTrials(Selection(j)).Rht.General{i}(1,1)) (DynamicPPTrials(Selection(j)).Rht.General{i}(2,2) - DynamicPPTrials(Selection(j)).Rht.General{i}(1,2))]);
%         DynamicPPTrials(Selection(j)).Rht.GeneralInV{i} = [-DynamicPPTrials(Selection(j)).Rht.GeneralV{i}(2) DynamicPPTrials(Selection(j)).Rht.GeneralV{i}(1)];  % calculate inverted slope of prog angles - for drawing perpendicular segmentation lines
%         % whole center of pressure points
%         DynamicPPTrials(Selection(j)).RProg.AngCoP{i} = round(abs(180/pi*(atan((DynamicPPTrials(Selection(j)).Rht.CoP{i}(2,1)-DynamicPPTrials(Selection(j)).Rht.CoP{i}(1,1))/(DynamicPPTrials(Selection(j)).Rht.CoP{i}(2,2)-DynamicPPTrials(Selection(j)).Rht.CoP{i}(1,2))))));
%         DynamicPPTrials(Selection(j)).RProg.AngActualCoP{i} = -180/pi*(atan((DynamicPPTrials(Selection(j)).Rht.CoP{i}(2,1)-DynamicPPTrials(Selection(j)).Rht.CoP{i}(1,1)) / (DynamicPPTrials(Selection(j)).Rht.CoP{i}(2,2)-DynamicPPTrials(Selection(j)).Rht.CoP{i}(1,2))));
%         DynamicPPTrials(Selection(j)).Rht.VCoP{i} = [(DynamicPPTrials(Selection(j)).Rht.CoP{i}(2,1)-DynamicPPTrials(Selection(j)).Rht.CoP{i}(1,1)) (DynamicPPTrials(Selection(j)).Rht.CoP{i}(2,2)-DynamicPPTrials(Selection(j)).Rht.CoP{i}(1,2))]  / norm([(DynamicPPTrials(Selection(j)).Rht.CoP{i}(2,1) - DynamicPPTrials(Selection(j)).Rht.CoP{i}(1,1)) (DynamicPPTrials(Selection(j)).Rht.CoP{i}(2,2) - DynamicPPTrials(Selection(j)).Rht.CoP{i}(1,2))]);
%         DynamicPPTrials(Selection(j)).Rht.InVCoP{i} = [-DynamicPPTrials(Selection(j)).Rht.VCoP{i}(2) DynamicPPTrials(Selection(j)).Rht.VCoP{i}(1)]; % calculate inverted slope of prog angles - for drawing perpendicular segmentation lines
%         % 66% CoP points
%         DynamicPPTrials(Selection(j)).RProg.Ang66{i} = round(abs(180/pi*(atan((DynamicPPTrials(Selection(j)).Rht.CoP66{i}(2,1)-DynamicPPTrials(Selection(j)).Rht.CoP66{i}(1,1))/(DynamicPPTrials(Selection(j)).Rht.CoP66{i}(2,2)-DynamicPPTrials(Selection(j)).Rht.CoP66{i}(1,2))))));
%         DynamicPPTrials(Selection(j)).RProg.AngActual66{i} = -180/pi*(atan((DynamicPPTrials(Selection(j)).Rht.CoP66{i}(2,1)-DynamicPPTrials(Selection(j)).Rht.CoP66{i}(1,1)) / (DynamicPPTrials(Selection(j)).Rht.CoP66{i}(2,2)-DynamicPPTrials(Selection(j)).Rht.CoP66{i}(1,2))));
%         DynamicPPTrials(Selection(j)).Rht.CoPV66{i} = [(DynamicPPTrials(Selection(j)).Rht.CoP66{i}(2,1)-DynamicPPTrials(Selection(j)).Rht.CoP66{i}(1,1)) (DynamicPPTrials(Selection(j)).Rht.CoP66{i}(2,2)-DynamicPPTrials(Selection(j)).Rht.CoP66{i}(1,2))]  / norm([(DynamicPPTrials(Selection(j)).Rht.CoP66{i}(2,1) - DynamicPPTrials(Selection(j)).Rht.CoP66{i}(1,1)) (DynamicPPTrials(Selection(j)).Rht.CoP66{i}(2,2) - DynamicPPTrials(Selection(j)).Rht.CoP66{i}(1,2))]);
%         DynamicPPTrials(Selection(j)).Rht.CoPInV66{i} = [-DynamicPPTrials(Selection(j)).Rht.CoPV66{i}(2) DynamicPPTrials(Selection(j)).Rht.CoPV66{i}(1)]; % calculate inverted slope of prog angles - for drawing perpendicular segmentation lines
%         % interpeak points
%         if Block.IP == 0
%             DynamicPPTrials(Selection(j)).RProg.AngIP{i} = round(abs(180/pi*(atan((DynamicPPTrials(Selection(j)).Rht.IP{i}(2,1)-DynamicPPTrials(Selection(j)).Rht.IP{i}(1,1))/(DynamicPPTrials(Selection(j)).Rht.IP{i}(2,2)-DynamicPPTrials(Selection(j)).Rht.IP{i}(1,2))))));
%             DynamicPPTrials(Selection(j)).RProg.AngActualIP{i} = -180/pi*(atan((DynamicPPTrials(Selection(j)).Rht.IP{i}(2,1)-DynamicPPTrials(Selection(j)).Rht.IP{i}(1,1)) / (DynamicPPTrials(Selection(j)).Rht.IP{i}(2,2)-DynamicPPTrials(Selection(j)).Rht.IP{i}(1,2))));
%             DynamicPPTrials(Selection(j)).Rht.VIP{i} = [(DynamicPPTrials(Selection(j)).Rht.IP{i}(2,1)-DynamicPPTrials(Selection(j)).Rht.IP{i}(1,1)) (DynamicPPTrials(Selection(j)).Rht.IP{i}(2,2)-DynamicPPTrials(Selection(j)).Rht.IP{i}(1,2))]  / norm([(DynamicPPTrials(Selection(j)).Rht.IP{i}(2,1) - DynamicPPTrials(Selection(j)).Rht.IP{i}(1,1)) (DynamicPPTrials(Selection(j)).Rht.IP{i}(2,2) - DynamicPPTrials(Selection(j)).Rht.IP{i}(1,2))]);
%             DynamicPPTrials(Selection(j)).Rht.InVIP{i} = [-DynamicPPTrials(Selection(j)).Rht.VIP{i}(2) DynamicPPTrials(Selection(j)).Rht.VIP{i}(1)]; % calculate inverted slope of prog angles - for drawing perpendicular segmentation lines
%         end
%         % Heel to centroid points
%         if Block.HC == 0
%             DynamicPPTrials(Selection(j)).RProg.AngHeelCent{i} = round(abs(180/pi*(atan((DynamicPPTrials(Selection(j)).Rht.HeelCent{i}(2,1)-DynamicPPTrials(Selection(j)).Rht.HeelCent{i}(1,1))/(DynamicPPTrials(Selection(j)).Rht.HeelCent{i}(2,2)-DynamicPPTrials(Selection(j)).Rht.HeelCent{i}(1,2))))));
%             DynamicPPTrials(Selection(j)).RProg.AngActualHeelCent{i} = -180/pi*(atan((DynamicPPTrials(Selection(j)).Rht.HeelCent{i}(2,1)-DynamicPPTrials(Selection(j)).Rht.HeelCent{i}(1,1)) / (DynamicPPTrials(Selection(j)).Rht.HeelCent{i}(2,2)-DynamicPPTrials(Selection(j)).Rht.HeelCent{i}(1,2))));
%             DynamicPPTrials(Selection(j)).Rht.VHeelCent{i} = [(DynamicPPTrials(Selection(j)).Rht.HeelCent{i}(2,1)-DynamicPPTrials(Selection(j)).Rht.HeelCent{i}(1,1)) (DynamicPPTrials(Selection(j)).Rht.HeelCent{i}(2,2)-DynamicPPTrials(Selection(j)).Rht.HeelCent{i}(1,2))]  / norm([(DynamicPPTrials(Selection(j)).Rht.HeelCent{i}(2,1) - DynamicPPTrials(Selection(j)).Rht.HeelCent{i}(1,1)) (DynamicPPTrials(Selection(j)).Rht.HeelCent{i}(2,2) - DynamicPPTrials(Selection(j)).Rht.HeelCent{i}(1,2))]);
%             DynamicPPTrials(Selection(j)).Rht.InVHeelCent{i} = [-DynamicPPTrials(Selection(j)).Rht.VHeelCent{i}(2) DynamicPPTrials(Selection(j)).Rht.VHeelCent{i}(1)]; % calculate inverted slope of prog angles - for drawing perpendicular segmentation lines
%         end
%     end
% end
% 
% %% Prepare for box searching
% for j = 1:length(Selection)
%     DynamicPPTrials(Selection(j)).MainMask = DynamicPPTrials(Selection(j)).SumTM;
%     DynamicPPTrials(Selection(j)).MaskLog = logical(DynamicPPTrials(Selection(j)).SumTM);
% end
% 
% % redefine heel and toe points to entire length of foot
% % extract bounding boxes for each foot strike
% for j = 1:length(Selection)
%     for i = 1:length(Regions{Selection(j)})
%         logs(i) = strcmp(Regions{1,Selection(j)}(i).Side, 'Left') == 1;
%     end
%     % LEFT foot timing
%     logind = find(logs);
%     for k = 1:sum(logs)
%         DynamicPPTrials(Selection(j)).BoundingBox(k).Left = Regions{1,Selection(j)}(logind(k)).BoundingBox;
%         DynamicPPTrials(Selection(j)).Centroid(k).Left = Regions{1,Selection(j)}(logind(k)).Centroid;
%         DynamicPPTrials(Selection(j)).TotalArea(k).Left = Regions{1,Selection(j)}(logind(k)).Area;
%     end
%     % RIGHT foot timing
%     logind = find(~logs);
%     for k = 1:length(logind)
%         DynamicPPTrials(Selection(j)).BoundingBox(k).Right = Regions{1,Selection(j)}(logind(k)).BoundingBox;
%         DynamicPPTrials(Selection(j)).Centroid(k).Right = Regions{1,Selection(j)}(logind(k)).Centroid;
%         DynamicPPTrials(Selection(j)).TotalArea(k).Right = Regions{1,Selection(j)}(logind(k)).Area;
%     end
%     clearvars logs logind
% end
% 
% % find intersection points of prog line with the horizontal front and back bounding boxes of each foot
% for q = 1:length(Selection)
%     % LEFT
%     for j = 1: DynamicPPTrials(Selection(q)).NumLeft
%         if DynamicPPTrials(Selection(q)).PoorPressure.Left(j) == 0
%             % manual points
%             if strcmp(PPSettings.MaskChoice, 'Manual') || strcmp(PPSettings.MaskChoice, 'Manual (pink)') || strcmp(PPSettings.MaskChoice, 'Validation')
%                 if abs(DynamicPPTrials(Selection(q)).LProg.AngManual{j}) < 45 % dont perform QA on points if prog angle is > 45
%                     [HTout] = FindIntersect(DynamicPPTrials(Selection(q)).Lht.Manual{j}, DynamicPPTrials(Selection(q)).BoundingBox(j).Left, DynamicPPTrials(Selection(q)).LProg.AngManual{j});
%                     DynamicPPTrials(Selection(q)).Lht.Manual{j} = HTout;
%                     clearvars HTout
%                 end
%             end
%             % general image processing points
%             if abs(DynamicPPTrials(Selection(q)).LProg.AngGeneral{j}) < 45 % dont perform QA on points if prog angle is > 45
%                 [HTout] = FindIntersect(DynamicPPTrials(Selection(q)).Lht.General{j}, DynamicPPTrials(Selection(q)).BoundingBox(j).Left, DynamicPPTrials(Selection(q)).LProg.AngGeneral{j});
%                 DynamicPPTrials(Selection(q)).Lht.General{j} = HTout;
%                 clearvars HTout
%             end
%             % CoP points
%             if abs(DynamicPPTrials(Selection(q)).LProg.AngCoP{j}) < 45 % dont perform QA on points if prog angle is > 45
%                 [HTout] = FindIntersect(DynamicPPTrials(Selection(q)).Lht.CoP{j}, DynamicPPTrials(Selection(q)).BoundingBox(j).Left, DynamicPPTrials(Selection(q)).LProg.AngCoP{j});
%                 DynamicPPTrials(Selection(q)).Lht.CoP{j} = HTout;
%                 clearvars HTout
%             end
%             % 66% CoP points
%             if abs(DynamicPPTrials(Selection(q)).LProg.Ang66{j}) < 45 % dont perform QA on points if prog angle is > 45
%                 [HTout] = FindIntersect(DynamicPPTrials(Selection(q)).Lht.CoP66{j}, DynamicPPTrials(Selection(q)).BoundingBox(j).Left, DynamicPPTrials(Selection(q)).LProg.Ang66{j});
%                 DynamicPPTrials(Selection(q)).Lht.CoP66{j} = HTout;
%                 clearvars HTout
%             end
%             % IP points
%             if Block.IP == 0
%                 if abs(DynamicPPTrials(Selection(q)).LProg.AngIP{j}) < 45 % dont perform QA on points if prog angle is > 45
%                     [HTout] = FindIntersect(DynamicPPTrials(Selection(q)).Lht.IP{j}, DynamicPPTrials(Selection(q)).BoundingBox(j).Left, DynamicPPTrials(Selection(q)).LProg.AngIP{j});
%                     DynamicPPTrials(Selection(q)).Lht.IP{j} = HTout;
%                     clearvars HTout
%                 end
%             end
%             % Heel to centroid points
%             if Block.HC == 0
%                 if abs(DynamicPPTrials(Selection(q)).LProg.AngHeelCent{j}) < 45 % dont perform QA on points if prog angle is > 45
%                     [HTout] = FindIntersect(DynamicPPTrials(Selection(q)).Lht.HeelCent{j}, DynamicPPTrials(Selection(q)).BoundingBox(j).Left, DynamicPPTrials(Selection(q)).LProg.AngHeelCent{j});
%                     DynamicPPTrials(Selection(q)).Lht.HeelCent{j} = HTout;
%                     clearvars HTout
%                 end
%             end
%         end
%     end
%     
%     %    RIGHT
%     for j = 1: DynamicPPTrials(Selection(q)).NumRight
%         if DynamicPPTrials(Selection(q)).PoorPressure.Right(j) == 0
%             % manual points
%             if strcmp(PPSettings.MaskChoice, 'Manual') || strcmp(PPSettings.MaskChoice, 'Manual (pink)') || strcmp(PPSettings.MaskChoice, 'Validation')
%                 if abs(DynamicPPTrials(Selection(q)).RProg.AngManual{j}) < 45 % dont perform QA on points if prog angle is > 45
%                     [HTout] = FindIntersect(DynamicPPTrials(Selection(q)).Rht.Manual{j}, DynamicPPTrials(Selection(q)).BoundingBox(j).Right, DynamicPPTrials(Selection(q)).RProg.AngManual{j});
%                     DynamicPPTrials(Selection(q)).Rht.Manual{j} = HTout;
%                     clearvars HTout
%                 end
%             end
%             % general image processing points
%             if abs(DynamicPPTrials(Selection(q)).RProg.AngGeneral{j}) < 45 % dont perform QA on points if prog angle is > 45
%                 [HTout] = FindIntersect(DynamicPPTrials(Selection(q)).Rht.General{j}, DynamicPPTrials(Selection(q)).BoundingBox(j).Right, DynamicPPTrials(Selection(q)).RProg.AngGeneral{j});
%                 DynamicPPTrials(Selection(q)).Rht.General{j} = HTout;
%                 clearvars HTout
%             end
%             % CoP points
%             if abs(DynamicPPTrials(Selection(q)).RProg.AngCoP{j}) < 45 % dont perform QA on points if prog angle is > 45
%                 [HTout] = FindIntersect(DynamicPPTrials(Selection(q)).Rht.CoP{j}, DynamicPPTrials(Selection(q)).BoundingBox(j).Right, DynamicPPTrials(Selection(q)).RProg.AngCoP{j});
%                 DynamicPPTrials(Selection(q)).Rht.CoP{j} = HTout;
%                 clearvars HTout
%             end
%             % 66% CoP points
%             if abs(DynamicPPTrials(Selection(q)).RProg.Ang66{j}) < 45 % dont perform QA on points if prog angle is > 45
%                 [HTout] = FindIntersect(DynamicPPTrials(Selection(q)).Rht.CoP66{j}, DynamicPPTrials(Selection(q)).BoundingBox(j).Right, DynamicPPTrials(Selection(q)).RProg.Ang66{j});
%                 DynamicPPTrials(Selection(q)).Rht.CoP66{j} = HTout;
%                 clearvars HTout
%             end
%             % IP points
%             if Block.IP == 0
%                 if abs(DynamicPPTrials(Selection(q)).RProg.AngIP{j}) < 45 % dont perform QA on points if prog angle is > 45
%                     [HTout] = FindIntersect(DynamicPPTrials(Selection(q)).Rht.IP{j}, DynamicPPTrials(Selection(q)).BoundingBox(j).Right, DynamicPPTrials(Selection(q)).RProg.AngIP{j});
%                     DynamicPPTrials(Selection(q)).Rht.IP{j} = HTout;
%                     clearvars HTout
%                 end
%             end
%             % Heel to centroid points
%             if Block.HC == 0
%                 if abs(DynamicPPTrials(Selection(q)).RProg.AngHeelCent{j}) < 45 % dont perform QA on points if prog angle is > 45
%                     [HTout] = FindIntersect(DynamicPPTrials(Selection(q)).Rht.HeelCent{j}, DynamicPPTrials(Selection(q)).BoundingBox(j).Right, DynamicPPTrials(Selection(q)).RProg.AngHeelCent{j});
%                     DynamicPPTrials(Selection(q)).Rht.HeelCent{j} = HTout;
%                     clearvars HTout
%                 end
%             end
%         end
%     end
% end
% 
% %% Search for borders of the foot
% clearvars i j k L F s S R q z M LRTest
% % parallel searching for the medial and lateral edges
% clc; disp('Searching for medial/lateral foot boundaries');
% [DynamicPPTrials] = SearchParallel(DynamicPPTrials, Selection, Block, PPSettings);
% 
% % perpendicular searching for anterior and posterior edges
% clc; disp('Searching for anterior/posterior foot boundaries');
% [DynamicPPTrials] = SearchPerpendicular(DynamicPPTrials, Selection, Block, PPSettings);
    
% %% Redefine heel-toe points at the shortened regions
% for q = 1:length(Selection)
%     %LEFT
%     for j = 1:DynamicPPTrials(Selection(q)).NumLeft
%         % manual points
%         if strcmp(PPSettings.MaskChoice, 'Manual') || strcmp(PPSettings.MaskChoice, 'Manual (pink)') || strcmp(PPSettings.MaskChoice, 'Validation')
%             % change y location to the mean of med and lat borders
%             DynamicPPTrials(Selection(q)).Lht.Manual{j}(1,2) = mean([DynamicPPTrials(Selection(q)).L.MedManual{j}(1,2), DynamicPPTrials(Selection(q)).L.LatManual{j}(1,2)]);
%             DynamicPPTrials(Selection(q)).Lht.Manual{j}(2,2) = mean([DynamicPPTrials(Selection(q)).L.MedManual{j}(2,2), DynamicPPTrials(Selection(q)).L.LatManual{j}(2,2)]);
%             % change x location to the mean of med and lat borders
%             DynamicPPTrials(Selection(q)).Lht.Manual{j}(1,1) = mean([DynamicPPTrials(Selection(q)).L.MedManual{j}(1,1), DynamicPPTrials(Selection(q)).L.LatManual{j}(1,1)]);
%             DynamicPPTrials(Selection(q)).Lht.Manual{j}(2,1) = mean([DynamicPPTrials(Selection(q)).L.MedManual{j}(2,1), DynamicPPTrials(Selection(q)).L.LatManual{j}(2,1)]);
%         end
%         % general image processing points
%         % change y location to the mean of med and lat borders
%         DynamicPPTrials(Selection(q)).Lht.General{j}(1,2) = mean([DynamicPPTrials(Selection(q)).L.MedGeneral{j}(1,2), DynamicPPTrials(Selection(q)).L.LatGeneral{j}(1,2)]);
%         DynamicPPTrials(Selection(q)).Lht.General{j}(2,2) = mean([DynamicPPTrials(Selection(q)).L.MedGeneral{j}(2,2), DynamicPPTrials(Selection(q)).L.LatGeneral{j}(2,2)]);
%         % change x location to the mean of med and lat borders
%         DynamicPPTrials(Selection(q)).Lht.General{j}(1,1) = mean([DynamicPPTrials(Selection(q)).L.MedGeneral{j}(1,1), DynamicPPTrials(Selection(q)).L.LatGeneral{j}(1,1)]);
%         DynamicPPTrials(Selection(q)).Lht.General{j}(2,1) = mean([DynamicPPTrials(Selection(q)).L.MedGeneral{j}(2,1), DynamicPPTrials(Selection(q)).L.LatGeneral{j}(2,1)]);
%         
%         % CoP points
%         % change y location to the mean of med and lat borders
%         DynamicPPTrials(Selection(q)).Lht.CoP{j}(1,2) = mean([DynamicPPTrials(Selection(q)).L.MedCoP{j}(1,2), DynamicPPTrials(Selection(q)).L.LatCoP{j}(1,2)]);
%         DynamicPPTrials(Selection(q)).Lht.CoP{j}(2,2) = mean([DynamicPPTrials(Selection(q)).L.MedCoP{j}(2,2), DynamicPPTrials(Selection(q)).L.LatCoP{j}(2,2)]);
%         % change x location to the mean of med and lat borders
%         DynamicPPTrials(Selection(q)).Lht.CoP{j}(1,1) = mean([DynamicPPTrials(Selection(q)).L.MedCoP{j}(1,1), DynamicPPTrials(Selection(q)).L.LatCoP{j}(1,1)]);
%         DynamicPPTrials(Selection(q)).Lht.CoP{j}(2,1) = mean([DynamicPPTrials(Selection(q)).L.MedCoP{j}(2,1), DynamicPPTrials(Selection(q)).L.LatCoP{j}(2,1)]);
%         
%         % 66% CoP points
%         % change y location to the mean of med and lat borders
%         DynamicPPTrials(Selection(q)).Lht.CoP66{j}(1,2) = mean([DynamicPPTrials(Selection(q)).L.Med66{j}(1,2), DynamicPPTrials(Selection(q)).L.Lat66{j}(1,2)]);
%         DynamicPPTrials(Selection(q)).Lht.CoP66{j}(2,2) = mean([DynamicPPTrials(Selection(q)).L.Med66{j}(2,2), DynamicPPTrials(Selection(q)).L.Lat66{j}(2,2)]);
%         % change x location to the mean of med and lat borders
%         DynamicPPTrials(Selection(q)).Lht.CoP66{j}(1,1) = mean([DynamicPPTrials(Selection(q)).L.Med66{j}(1,1), DynamicPPTrials(Selection(q)).L.Lat66{j}(1,1)]);
%         DynamicPPTrials(Selection(q)).Lht.CoP66{j}(2,1) = mean([DynamicPPTrials(Selection(q)).L.Med66{j}(2,1), DynamicPPTrials(Selection(q)).L.Lat66{j}(2,1)]);
%         
%         % IP points
%         if Block.IP == 0
%             % change y location to the mean of med and lat borders
%             DynamicPPTrials(Selection(q)).Lht.IP{j}(1,2) = mean([DynamicPPTrials(Selection(q)).L.MedIP{j}(1,2), DynamicPPTrials(Selection(q)).L.LatIP{j}(1,2)]);
%             DynamicPPTrials(Selection(q)).Lht.IP{j}(2,2) = mean([DynamicPPTrials(Selection(q)).L.MedIP{j}(2,2), DynamicPPTrials(Selection(q)).L.LatIP{j}(2,2)]);
%             % change x location to the mean of med and lat borders
%             DynamicPPTrials(Selection(q)).Lht.IP{j}(1,1) = mean([DynamicPPTrials(Selection(q)).L.MedIP{j}(1,1), DynamicPPTrials(Selection(q)).L.LatIP{j}(1,1)]);
%             DynamicPPTrials(Selection(q)).Lht.IP{j}(2,1) = mean([DynamicPPTrials(Selection(q)).L.MedIP{j}(2,1), DynamicPPTrials(Selection(q)).L.LatIP{j}(2,1)]);
%         end
%         
%         % Heel to centroid points
%         if Block.HC == 0
%             % change y location to the mean of med and lat borders
%             DynamicPPTrials(Selection(q)).Lht.HeelCent{j}(1,2) = mean([DynamicPPTrials(Selection(q)).L.MedHeelCent{j}(1,2), DynamicPPTrials(Selection(q)).L.LatHeelCent{j}(1,2)]);
%             DynamicPPTrials(Selection(q)).Lht.HeelCent{j}(2,2) = mean([DynamicPPTrials(Selection(q)).L.MedHeelCent{j}(2,2), DynamicPPTrials(Selection(q)).L.LatHeelCent{j}(2,2)]);
%             % change x location to the mean of med and lat borders
%             DynamicPPTrials(Selection(q)).Lht.HeelCent{j}(1,1) = mean([DynamicPPTrials(Selection(q)).L.MedHeelCent{j}(1,1), DynamicPPTrials(Selection(q)).L.LatHeelCent{j}(1,1)]);
%             DynamicPPTrials(Selection(q)).Lht.HeelCent{j}(2,1) = mean([DynamicPPTrials(Selection(q)).L.MedHeelCent{j}(2,1), DynamicPPTrials(Selection(q)).L.LatHeelCent{j}(2,1)]);
%         end
%     end
%     % Uncomment to verify re-definition of heel-toe locations
%     %         for j = 1: DynamicPPTrials(Selection(q)).NumLeft
%     %             subplot(1,NumDynTrials,Selection(q)); hold on;
%     %             plot(DynamicPPTrials(Selection(q)).Lht.{j}(:,1), DynamicPPTrials(Selection(q)).Lht.{j}(:,2), 'k*');
%     %         end
%     
%     % RIGHT
%     for j = 1: DynamicPPTrials(Selection(q)).NumRight
%         % manual points
%         if strcmp(PPSettings.MaskChoice, 'Manual') || strcmp(PPSettings.MaskChoice, 'Manual (pink)') || strcmp(PPSettings.MaskChoice, 'Validation')
%             % change y location to the mean of med and lat borders
%             DynamicPPTrials(Selection(q)).Rht.Manual{j}(1,2) = mean([DynamicPPTrials(Selection(q)).R.MedManual{j}(1,2), DynamicPPTrials(Selection(q)).R.LatManual{j}(1,2)]);
%             DynamicPPTrials(Selection(q)).Rht.Manual{j}(2,2) = mean([DynamicPPTrials(Selection(q)).R.MedManual{j}(2,2), DynamicPPTrials(Selection(q)).R.LatManual{j}(2,2)]);
%             % change x location to the mean of med and lat borders
%             DynamicPPTrials(Selection(q)).Rht.Manual{j}(1,1) = mean([DynamicPPTrials(Selection(q)).R.MedManual{j}(1,1), DynamicPPTrials(Selection(q)).R.LatManual{j}(1,1)]);
%             DynamicPPTrials(Selection(q)).Rht.Manual{j}(2,1) = mean([DynamicPPTrials(Selection(q)).R.MedManual{j}(2,1), DynamicPPTrials(Selection(q)).R.LatManual{j}(2,1)]);
%         end
%         
%         % General image processing points
%         % change y location to the mean of med and lat borders
%         DynamicPPTrials(Selection(q)).Rht.General{j}(1,2) = mean([DynamicPPTrials(Selection(q)).R.MedGeneral{j}(1,2), DynamicPPTrials(Selection(q)).R.LatGeneral{j}(1,2)]);
%         DynamicPPTrials(Selection(q)).Rht.General{j}(2,2) = mean([DynamicPPTrials(Selection(q)).R.MedGeneral{j}(2,2), DynamicPPTrials(Selection(q)).R.LatGeneral{j}(2,2)]);
%         % change x location to the mean of med and lat borders
%         DynamicPPTrials(Selection(q)).Rht.General{j}(1,1) = mean([DynamicPPTrials(Selection(q)).R.MedGeneral{j}(1,1), DynamicPPTrials(Selection(q)).R.LatGeneral{j}(1,1)]);
%         DynamicPPTrials(Selection(q)).Rht.General{j}(2,1) = mean([DynamicPPTrials(Selection(q)).R.MedGeneral{j}(2,1), DynamicPPTrials(Selection(q)).R.LatGeneral{j}(2,1)]);
%         
%         % CoP points
%         % change y location to the mean of med and lat borders
%         DynamicPPTrials(Selection(q)).Rht.CoP{j}(1,2) = mean([DynamicPPTrials(Selection(q)).R.MedCoP{j}(1,2), DynamicPPTrials(Selection(q)).R.LatCoP{j}(1,2)]);
%         DynamicPPTrials(Selection(q)).Rht.CoP{j}(2,2) = mean([DynamicPPTrials(Selection(q)).R.MedCoP{j}(2,2), DynamicPPTrials(Selection(q)).R.LatCoP{j}(2,2)]);
%         % change y location to the mean of med and lat borders
%         DynamicPPTrials(Selection(q)).Rht.CoP{j}(1,1) = mean([DynamicPPTrials(Selection(q)).R.MedCoP{j}(1,1), DynamicPPTrials(Selection(q)).R.LatCoP{j}(1,1)]);
%         DynamicPPTrials(Selection(q)).Rht.CoP{j}(2,1) = mean([DynamicPPTrials(Selection(q)).R.MedCoP{j}(2,1), DynamicPPTrials(Selection(q)).R.LatCoP{j}(2,1)]);
%         
%         % 66% CoP points
%         % change y location to the mean of med and lat borders
%         DynamicPPTrials(Selection(q)).Rht.CoP66{j}(1,2) = mean([DynamicPPTrials(Selection(q)).R.Med66{j}(1,2), DynamicPPTrials(Selection(q)).R.Lat66{j}(1,2)]);
%         DynamicPPTrials(Selection(q)).Rht.CoP66{j}(2,2) = mean([DynamicPPTrials(Selection(q)).R.Med66{j}(2,2), DynamicPPTrials(Selection(q)).R.Lat66{j}(2,2)]);
%         % change x location to the mean of med and lat borders
%         DynamicPPTrials(Selection(q)).Rht.CoP66{j}(1,1) = mean([DynamicPPTrials(Selection(q)).R.Med66{j}(1,1), DynamicPPTrials(Selection(q)).R.Lat66{j}(1,1)]);
%         DynamicPPTrials(Selection(q)).Rht.CoP66{j}(2,1) = mean([DynamicPPTrials(Selection(q)).R.Med66{j}(2,1), DynamicPPTrials(Selection(q)).R.Lat66{j}(2,1)]);
%         
%         % IP points
%         if Block.IP == 0
%             % change y location to the mean of med and lat borders
%             DynamicPPTrials(Selection(q)).Rht.IP{j}(1,2) = mean([DynamicPPTrials(Selection(q)).R.MedIP{j}(1,2), DynamicPPTrials(Selection(q)).R.LatIP{j}(1,2)]);
%             DynamicPPTrials(Selection(q)).Rht.IP{j}(2,2) = mean([DynamicPPTrials(Selection(q)).R.MedIP{j}(2,2), DynamicPPTrials(Selection(q)).R.LatIP{j}(2,2)]);
%             % change x location to the mean of med and lat borders
%             DynamicPPTrials(Selection(q)).Rht.IP{j}(1,1) = mean([DynamicPPTrials(Selection(q)).R.MedIP{j}(1,1), DynamicPPTrials(Selection(q)).R.LatIP{j}(1,1)]);
%             DynamicPPTrials(Selection(q)).Rht.IP{j}(2,1) = mean([DynamicPPTrials(Selection(q)).R.MedIP{j}(2,1), DynamicPPTrials(Selection(q)).R.LatIP{j}(2,1)]);
%         end
%         
%         % Heel to centroid points
%         if Block.HC == 0
%             % change y location to the mean of med and lat borders
%             DynamicPPTrials(Selection(q)).Rht.HeelCent{j}(1,2) = mean([DynamicPPTrials(Selection(q)).R.MedHeelCent{j}(1,2), DynamicPPTrials(Selection(q)).R.LatHeelCent{j}(1,2)]);
%             DynamicPPTrials(Selection(q)).Rht.HeelCent{j}(2,2) = mean([DynamicPPTrials(Selection(q)).R.MedHeelCent{j}(2,2), DynamicPPTrials(Selection(q)).R.LatHeelCent{j}(2,2)]);
%             % change y location to the mean of med and lat borders
%             DynamicPPTrials(Selection(q)).Rht.HeelCent{j}(1,1) = mean([DynamicPPTrials(Selection(q)).R.MedHeelCent{j}(1,1), DynamicPPTrials(Selection(q)).R.LatHeelCent{j}(1,1)]);
%             DynamicPPTrials(Selection(q)).Rht.HeelCent{j}(2,1) = mean([DynamicPPTrials(Selection(q)).R.MedHeelCent{j}(2,1), DynamicPPTrials(Selection(q)).R.LatHeelCent{j}(2,1)]);
%         end
%     end
%     % Uncomment to verify re-definition of heel-toe locations
%     %         for j = 1: DynamicPPTrials(Selection(q)).NumRight
%     %             subplot(1,NumDynTrials,Selection(q)); hold on;
%     %             plot(DynamicPPTrials(Selection(q)).Rht.General{j}(:,1), DynamicPPTrials(Selection(q)).Rht.General{j}(:,2), 'k*');
%     %         end
% end
% clearvars xH yH xT yT xFP yFP IntersectHFP IntersectTFP x q d iyH iyT iyFP ixH ixT m n p pT pFP pH Srch i j k Reg W z yTval yHval Num2Del PkThresh s S F
% 
% %% optimize C3D locations
if strcmp(PPSettings.C3DInput,'Yes')
    if exist('ToDel', 'var')
        C3Ddata(ToDel) = [];
    end
    clearvars Dist ToDel
    DelThresh = 1; % threshold for optimizing centroid location = 0.5 cm (1 half-cm)
    for z = 1:length(C3Ddata) % optimize C3D locations by finding difference between the foot center from C3D and the average of all masks
        j = 1;
        % LEFT
        if isempty(DynamicPPTrials(z).Lht) == 0
            if isempty(DynamicPPTrials(z).Lht.General) == 0
                if Block.IP == 0 && Block.HC == 0
                    AvgX = [mean(DynamicPPTrials(z).Lht.Manual{j}(:,1)), mean(DynamicPPTrials(z).Lht.General{j}(:,1)), mean(DynamicPPTrials(z).Lht.CoP{j}(:,1)),...
                        mean(DynamicPPTrials(z).Lht.CoP66{j}(:,1)), mean(DynamicPPTrials(z).Lht.IP{j}(:,1)), mean(DynamicPPTrials(z).Lht.HeelCent{j}(:,1))];
                    AvgY = [mean(DynamicPPTrials(z).Lht.Manual{j}(:,2)), mean(DynamicPPTrials(z).Lht.General{j}(:,2)), mean(DynamicPPTrials(z).Lht.CoP{j}(:,2)),...
                        mean(DynamicPPTrials(z).Lht.CoP66{j}(:,2)), mean(DynamicPPTrials(z).Lht.IP{j}(:,2)), mean(DynamicPPTrials(z).Lht.HeelCent{j}(:,2))];
                elseif Block.IP == 1 && Block.HC == 0
                     AvgX = [mean(DynamicPPTrials(z).Lht.Manual{j}(:,1)), mean(DynamicPPTrials(z).Lht.General{j}(:,1)), mean(DynamicPPTrials(z).Lht.CoP{j}(:,1)),...
                        mean(DynamicPPTrials(z).Lht.CoP66{j}(:,1)),  mean(DynamicPPTrials(z).Lht.HeelCent{j}(:,1))];
                    AvgY = [mean(DynamicPPTrials(z).Lht.Manual{j}(:,2)), mean(DynamicPPTrials(z).Lht.General{j}(:,2)), mean(DynamicPPTrials(z).Lht.CoP{j}(:,2)),...
                        mean(DynamicPPTrials(z).Lht.CoP66{j}(:,2)), mean(DynamicPPTrials(z).Lht.HeelCent{j}(:,2))];
                elseif Block.IP == 0 && Block.HC == 1
                     AvgX = [mean(DynamicPPTrials(z).Lht.Manual{j}(:,1)), mean(DynamicPPTrials(z).Lht.General{j}(:,1)), mean(DynamicPPTrials(z).Lht.CoP{j}(:,1)),...
                        mean(DynamicPPTrials(z).Lht.CoP66{j}(:,1)), mean(DynamicPPTrials(z).Lht.IP{j}(:,1))];
                    AvgY = [mean(DynamicPPTrials(z).Lht.Manual{j}(:,2)), mean(DynamicPPTrials(z).Lht.General{j}(:,2)), mean(DynamicPPTrials(z).Lht.CoP{j}(:,2)),...
                        mean(DynamicPPTrials(z).Lht.CoP66{j}(:,2)), mean(DynamicPPTrials(z).Lht.IP{j}(:,2))];
                elseif Block.IP == 1 && Block.HC == 1
                  AvgX = [mean(DynamicPPTrials(z).Lht.Manual{j}(:,1)), mean(DynamicPPTrials(z).Lht.General{j}(:,1)), mean(DynamicPPTrials(z).Lht.CoP{j}(:,1)),...
                        mean(DynamicPPTrials(z).Lht.CoP66{j}(:,1))];
                    AvgY = [mean(DynamicPPTrials(z).Lht.Manual{j}(:,2)), mean(DynamicPPTrials(z).Lht.General{j}(:,2)), mean(DynamicPPTrials(z).Lht.CoP{j}(:,2)),...
                        mean(DynamicPPTrials(z).Lht.CoP66{j}(:,2))];
                end
                Dist(z).LeftX = mean(AvgX) - C3Ddata(z).LCentB(j,2);
                Dist(z).LeftY = mean(AvgY) - C3Ddata(z).LCentB(j,1);
            end
        end
        % RIGHT
        if isempty(DynamicPPTrials(z).Rht) == 0
            if isempty(DynamicPPTrials(z).Rht.General) == 0
                if Block.IP == 0 && Block.HC == 0
                    AvgX = [mean(DynamicPPTrials(z).Rht.Manual{j}(:,1)), mean(DynamicPPTrials(z).Rht.General{j}(:,1)), mean(DynamicPPTrials(z).Rht.CoP{j}(:,1)),...
                        mean(DynamicPPTrials(z).Rht.CoP66{j}(:,1)), mean(DynamicPPTrials(z).Rht.IP{j}(:,1)), mean(DynamicPPTrials(z).Rht.HeelCent{j}(:,1))];
                    AvgY = [mean(DynamicPPTrials(z).Rht.Manual{j}(:,2)), mean(DynamicPPTrials(z).Rht.General{j}(:,2)), mean(DynamicPPTrials(z).Rht.CoP{j}(:,2)),...
                        mean(DynamicPPTrials(z).Rht.CoP66{j}(:,2)), mean(DynamicPPTrials(z).Rht.IP{j}(:,2)), mean(DynamicPPTrials(z).Rht.HeelCent{j}(:,2))];
                elseif Block.IP == 1 && Block.HC == 0
                    AvgX = [mean(DynamicPPTrials(z).Rht.Manual{j}(:,1)), mean(DynamicPPTrials(z).Rht.General{j}(:,1)), mean(DynamicPPTrials(z).Rht.CoP{j}(:,1)),...
                        mean(DynamicPPTrials(z).Rht.CoP66{j}(:,1)),  mean(DynamicPPTrials(z).Rht.HeelCent{j}(:,1))];
                    AvgY = [mean(DynamicPPTrials(z).Rht.Manual{j}(:,2)), mean(DynamicPPTrials(z).Rht.General{j}(:,2)), mean(DynamicPPTrials(z).Rht.CoP{j}(:,2)),...
                        mean(DynamicPPTrials(z).Rht.CoP66{j}(:,2)), mean(DynamicPPTrials(z).Rht.HeelCent{j}(:,2))];
                elseif Block.IP == 0 && Block.HC == 1
                    AvgX = [mean(DynamicPPTrials(z).Rht.Manual{j}(:,1)), mean(DynamicPPTrials(z).Rht.General{j}(:,1)), mean(DynamicPPTrials(z).Rht.CoP{j}(:,1)),...
                        mean(DynamicPPTrials(z).Rht.CoP66{j}(:,1)), mean(DynamicPPTrials(z).Rht.IP{j}(:,1))];
                    AvgY = [mean(DynamicPPTrials(z).Rht.Manual{j}(:,2)), mean(DynamicPPTrials(z).Rht.General{j}(:,2)), mean(DynamicPPTrials(z).Rht.CoP{j}(:,2)),...
                        mean(DynamicPPTrials(z).Rht.CoP66{j}(:,2)), mean(DynamicPPTrials(z).Rht.IP{j}(:,2))];
                elseif Block.IP == 1 && Block.HC == 1
                    AvgX = [mean(DynamicPPTrials(z).Rht.Manual{j}(:,1)), mean(DynamicPPTrials(z).Rht.General{j}(:,1)), mean(DynamicPPTrials(z).Rht.CoP{j}(:,1)),...
                        mean(DynamicPPTrials(z).Rht.CoP66{j}(:,1))];
                    AvgY = [mean(DynamicPPTrials(z).Rht.Manual{j}(:,2)), mean(DynamicPPTrials(z).Rht.General{j}(:,2)), mean(DynamicPPTrials(z).Rht.CoP{j}(:,2)),...
                        mean(DynamicPPTrials(z).Rht.CoP66{j}(:,2))];
                end
                Dist(z).RightX = mean(AvgX) - C3Ddata(z).RCentB(j,2);
                Dist(z).RightY = mean(AvgY) - C3Ddata(z).RCentB(j,1);
            end
        end

        % determine side-to-side offset between C3D and MoCap
        if isnan(Dist(z).LeftX) == 0
            X = nanmean(Dist(z).LeftX);
            if abs(X) > DelThresh % correct if mean difference is greater than threshold
                C3Ddata(z).Adjust.LeftX = X;
                if isempty(C3Ddata(z).LCent) == 0
                    for j = 1:size(C3Ddata(z).LHA, 1)
                        C3Ddata(z).LHA(j, 2) = C3Ddata(z).LHA(j, 2) + X;
                        C3Ddata(z).LFA(j, 2) = C3Ddata(z).LFA(j, 2) + X;
                        C3Ddata(z).LHAb(j, 2) = C3Ddata(z).LHAb(j, 2) + X;
                        C3Ddata(z).LFAb(j, 2) = C3Ddata(z).LFAb(j, 2) + X;
                        C3Ddata(z).LCent(j, 2) = C3Ddata(z).LCent(j, 2) + X;
                        C3Ddata(z).LCentB(j, 2) = C3Ddata(z).LCentB(j, 2) + X;
                        C3Ddata(z).LHEE.Adj(:,2) = C3Ddata(z).LHEE.Adj(:,2) + X;
                        C3Ddata(z).LHEE.base(j,2) = C3Ddata(z).LHEE.base(j,2) + X;
                        C3Ddata(z).LTOE.top(j,2) = C3Ddata(z).LTOE.top(j,2) + X;
                        C3Ddata(z).LTOE.step(j,2) = C3Ddata(z).LTOE.step(j,2) + X;
                        C3Ddata(z).LD1M.step(j,2) = C3Ddata(z).LD1M.step(j,2) + X;
                        C3Ddata(z).LD5M.step(j,2) = C3Ddata(z).LD5M.step(j,2) + X;
                        C3Ddata(z).LP1M.step(j,2) = C3Ddata(z).LP1M.step(j,2) + X;
                        C3Ddata(z).LP5M.step(j,2) = C3Ddata(z).LP5M.step(j,2) + X;
                        C3Ddata(z).LLCA.step(j,2) = C3Ddata(z).LLCA.step(j,2) + X;
                        C3Ddata(z).LMCA.step(j,2) = C3Ddata(z).LMCA.step(j,2) + X;
                    end
                end
            end
        end
        if isnan(Dist(z).RightX) == 0
            X = nanmean(Dist(z).RightX);
            if abs(X) > DelThresh % correct if mean difference is greater than threshold
                C3Ddata(z).Adjust.RightX = X;
                if isempty(C3Ddata(z).RCent) == 0
                    for j = 1:size(C3Ddata(z).RHA,1)
                        C3Ddata(z).RHA(j, 2) = C3Ddata(z).RHA(j, 2) + X;
                        C3Ddata(z).RFA(j, 2) = C3Ddata(z).RFA(j, 2) + X;
                        C3Ddata(z).RHAb(j, 2) = C3Ddata(z).RHAb(j, 2) + X;
                        C3Ddata(z).RFAb(j, 2) = C3Ddata(z).RFAb(j, 2) + X;
                        C3Ddata(z).RCent(j, 2) = C3Ddata(z).RCent(j, 2) + X;
                        C3Ddata(z).RCentB(j, 2) = C3Ddata(z).RCentB(j, 2) + X;
                        C3Ddata(z).RHEE.Adj(:,2) = C3Ddata(z).RHEE.Adj(:,2) + X;
                        C3Ddata(z).RHEE.base(j,2) = C3Ddata(z).RHEE.base(j,2) + X;
                        C3Ddata(z).RTOE.top(j,2) = C3Ddata(z).RTOE.top(j,2) + X;
                        C3Ddata(z).RTOE.step(j,2) = C3Ddata(z).RTOE.step(j,2) + X;
                        C3Ddata(z).RD1M.step(j,2) = C3Ddata(z).RD1M.step(j,2) + X;
                        C3Ddata(z).RD5M.step(j,2) = C3Ddata(z).RD5M.step(j,2) + X;
                        C3Ddata(z).RP1M.step(j,2) = C3Ddata(z).RP1M.step(j,2) + X;
                        C3Ddata(z).RP5M.step(j,2) = C3Ddata(z).RP5M.step(j,2) + X;
                        C3Ddata(z).RLCA.step(j,2) = C3Ddata(z).RLCA.step(j,2) + X;
                        C3Ddata(z).RMCA.step(j,2) = C3Ddata(z).RMCA.step(j,2) + X;
                    end
                end
            end
        end
        
        % determine up and down offset
        if isnan(Dist(z).LeftY) == 0
            Y = nanmean(Dist(z).LeftY);
            if abs(Y) > DelThresh % correct if mean difference is greater than threshold
                C3Ddata(z).Adjust.LeftY = Y;
                if isempty( C3Ddata(z).LCent) == 0
                    for j = 1:size(C3Ddata(z).LHA)
                        C3Ddata(z).LHA(j, 1) = C3Ddata(z).LHA(j, 1) + Y;
                        C3Ddata(z).LFA(j, 1) = C3Ddata(z).LFA(j, 1) + Y;
                        C3Ddata(z).LHAb(j, 1) = C3Ddata(z).LHAb(j, 1) + Y;
                        C3Ddata(z).LFAb(j, 1) = C3Ddata(z).LFAb(j, 1) + Y;
                        C3Ddata(z).LCent(j, 1) = C3Ddata(z).LCent(j, 1) + Y;
                        C3Ddata(z).LCentB(j, 1) = C3Ddata(z).LCentB(j, 1) + Y;
                        C3Ddata(z).LHEE.Adj(:,1) = C3Ddata(z).LHEE.Adj(:,1) + Y;
                        C3Ddata(z).LHEE.base(j,1) = C3Ddata(z).LHEE.base(j,1) + Y;
                        C3Ddata(z).LTOE.top(j,1) = C3Ddata(z).LTOE.top(j,1) + Y;
                        C3Ddata(z).LTOE.step(j,1) = C3Ddata(z).LTOE.step(j,1) + Y;
                        C3Ddata(z).LD1M.step(j,1) = C3Ddata(z).LD1M.step(j,1) + Y;
                        C3Ddata(z).LD5M.step(j,1) = C3Ddata(z).LD5M.step(j,1) + Y;
                        C3Ddata(z).LP1M.step(j,1) = C3Ddata(z).LP1M.step(j,1) + Y;
                        C3Ddata(z).LP5M.step(j,1) = C3Ddata(z).LP5M.step(j,1) + Y;
                        C3Ddata(z).LLCA.step(j,1) = C3Ddata(z).LLCA.step(j,1) + Y;
                        C3Ddata(z).LMCA.step(j,1) = C3Ddata(z).LMCA.step(j,1) + Y;
                    end
                end
            end
        end
        if isnan(Dist(z).RightY) == 0
            Y = nanmean(Dist(z).RightY);
            if abs(Y) > DelThresh % correct if mean difference is greater than threshold
                C3Ddata(z).Adjust.RightY = Y;
                if isempty(C3Ddata(z).RCent) == 0
                    for j = 1:size(C3Ddata(z).RHA)
                        C3Ddata(z).RHA(j, 1) = C3Ddata(z).RHA(j, 1) + Y;
                        C3Ddata(z).RFA(j, 1) = C3Ddata(z).RFA(j, 1) + Y;
                        C3Ddata(z).RHAb(j, 1) = C3Ddata(z).RHAb(j, 1) + Y;
                        C3Ddata(z).RFAb(j, 1) = C3Ddata(z).RFAb(j, 1) + Y;
                        C3Ddata(z).RCent(j, 1) = C3Ddata(z).RCent(j, 1) + Y;
                        C3Ddata(z).RCentB(j, 1) = C3Ddata(z).RCentB(j, 1) + Y;
                        C3Ddata(z).RHEE.Adj(:,1) = C3Ddata(z).RHEE.Adj(:,1) + Y;
                        C3Ddata(z).RHEE.base(j,1) = C3Ddata(z).RHEE.base(j,1) + Y;
                        C3Ddata(z).RTOE.top(j,1) = C3Ddata(z).RTOE.top(j,1) + Y;
                        C3Ddata(z).RTOE.step(j,1) = C3Ddata(z).RTOE.step(j,1) + Y;
                        C3Ddata(z).RD1M.step(j,1) = C3Ddata(z).RD1M.step(j,1) + Y;
                        C3Ddata(z).RD5M.step(j,1) = C3Ddata(z).RD5M.step(j,1) + Y;
                        C3Ddata(z).RP1M.step(j,1) = C3Ddata(z).RP1M.step(j,1) + Y;
                        C3Ddata(z).RP5M.step(j,1) = C3Ddata(z).RP5M.step(j,1) + Y;
                        C3Ddata(z).RLCA.step(j,1) = C3Ddata(z).RLCA.step(j,1) + Y;
                        C3Ddata(z).RMCA.step(j,1) = C3Ddata(z).RMCA.step(j,1) + Y;
                    end
                end
            end
        end
        clearvars X Y
    end
end
% 
% %% 1/3 sectioning
% clc; disp('Segmenting pressures into Heel/Arch/Forefoot and Medial/Lateral regions.');
% for q = 1:length(Selection)
%     % LEFT
%     for i = 1:DynamicPPTrials(Selection(q)).NumLeft
%         % find the 1/3 of foot prog line points
%         if strcmp(PPSettings.MaskChoice, 'Manual') || strcmp(PPSettings.MaskChoice, 'Manual (pink)') || strcmp(PPSettings.MaskChoice, 'Validation')
%             DynamicPPTrials(Selection(q)).LProg.HAManual{i} = [DynamicPPTrials(Selection(q)).Lht.Manual{i}(1,1) + ((DynamicPPTrials(Selection(q)).Lht.Manual{i}(2,1) - DynamicPPTrials(Selection(q)).Lht.Manual{i}(1,1))*0.3),  DynamicPPTrials(Selection(q)).Lht.Manual{i}(1,2) + ((DynamicPPTrials(Selection(q)).Lht.Manual{i}(2,2) - DynamicPPTrials(Selection(q)).Lht.Manual{i}(1,2))*0.3)];
%             DynamicPPTrials(Selection(q)).LProg.FAManual{i} = [DynamicPPTrials(Selection(q)).Lht.Manual{i}(1,1) + 2*((DynamicPPTrials(Selection(q)).Lht.Manual{i}(2,1) - DynamicPPTrials(Selection(q)).Lht.Manual{i}(1,1))*0.3),  DynamicPPTrials(Selection(q)).Lht.Manual{i}(1,2) + 2*((DynamicPPTrials(Selection(q)).Lht.Manual{i}(2,2) - DynamicPPTrials(Selection(q)).Lht.Manual{i}(1,2))*0.3)];
%         end
%         DynamicPPTrials(Selection(q)).LProg.HAGeneral{i} = [DynamicPPTrials(Selection(q)).Lht.General{i}(1,1) + ((DynamicPPTrials(Selection(q)).Lht.General{i}(2,1) - DynamicPPTrials(Selection(q)).Lht.General{i}(1,1))*0.3),  DynamicPPTrials(Selection(q)).Lht.General{i}(1,2) + ((DynamicPPTrials(Selection(q)).Lht.General{i}(2,2) - DynamicPPTrials(Selection(q)).Lht.General{i}(1,2))*0.3)];
%         DynamicPPTrials(Selection(q)).LProg.FAGeneral{i} = [DynamicPPTrials(Selection(q)).Lht.General{i}(1,1) + 2*((DynamicPPTrials(Selection(q)).Lht.General{i}(2,1) - DynamicPPTrials(Selection(q)).Lht.General{i}(1,1))*0.3),  DynamicPPTrials(Selection(q)).Lht.General{i}(1,2) + 2*((DynamicPPTrials(Selection(q)).Lht.General{i}(2,2) - DynamicPPTrials(Selection(q)).Lht.General{i}(1,2))*0.3)];
%         DynamicPPTrials(Selection(q)).LProg.HACoP{i} = [DynamicPPTrials(Selection(q)).Lht.CoP{i}(1,1) + ((DynamicPPTrials(Selection(q)).Lht.CoP{i}(2,1) - DynamicPPTrials(Selection(q)).Lht.CoP{i}(1,1))*0.3),  DynamicPPTrials(Selection(q)).Lht.CoP{i}(1,2) + ((DynamicPPTrials(Selection(q)).Lht.CoP{i}(2,2) - DynamicPPTrials(Selection(q)).Lht.CoP{i}(1,2))*0.3)];
%         DynamicPPTrials(Selection(q)).LProg.FACoP{i} = [DynamicPPTrials(Selection(q)).Lht.CoP{i}(1,1) + 2*((DynamicPPTrials(Selection(q)).Lht.CoP{i}(2,1) - DynamicPPTrials(Selection(q)).Lht.CoP{i}(1,1))*0.3),  DynamicPPTrials(Selection(q)).Lht.CoP{i}(1,2) + 2*((DynamicPPTrials(Selection(q)).Lht.CoP{i}(2,2) - DynamicPPTrials(Selection(q)).Lht.CoP{i}(1,2))*0.3)];
%         DynamicPPTrials(Selection(q)).LProg.HA66{i} = [DynamicPPTrials(Selection(q)).Lht.CoP66{i}(1,1) + ((DynamicPPTrials(Selection(q)).Lht.CoP66{i}(2,1) - DynamicPPTrials(Selection(q)).Lht.CoP66{i}(1,1))*0.3),  DynamicPPTrials(Selection(q)).Lht.CoP66{i}(1,2) + ((DynamicPPTrials(Selection(q)).Lht.CoP66{i}(2,2) - DynamicPPTrials(Selection(q)).Lht.CoP66{i}(1,2))*0.3)];
%         DynamicPPTrials(Selection(q)).LProg.FA66{i} = [DynamicPPTrials(Selection(q)).Lht.CoP66{i}(1,1) + 2*((DynamicPPTrials(Selection(q)).Lht.CoP66{i}(2,1) - DynamicPPTrials(Selection(q)).Lht.CoP66{i}(1,1))*0.3),  DynamicPPTrials(Selection(q)).Lht.CoP66{i}(1,2) + 2*((DynamicPPTrials(Selection(q)).Lht.CoP66{i}(2,2) - DynamicPPTrials(Selection(q)).Lht.CoP66{i}(1,2))*0.3)];
%         if Block.IP == 0
%             DynamicPPTrials(Selection(q)).LProg.HAIP{i} = [DynamicPPTrials(Selection(q)).Lht.IP{i}(1,1) + ((DynamicPPTrials(Selection(q)).Lht.IP{i}(2,1) - DynamicPPTrials(Selection(q)).Lht.IP{i}(1,1))*0.3),  DynamicPPTrials(Selection(q)).Lht.IP{i}(1,2) + ((DynamicPPTrials(Selection(q)).Lht.IP{i}(2,2) - DynamicPPTrials(Selection(q)).Lht.IP{i}(1,2))*0.3)];
%             DynamicPPTrials(Selection(q)).LProg.FAIP{i} = [DynamicPPTrials(Selection(q)).Lht.IP{i}(1,1) + 2*((DynamicPPTrials(Selection(q)).Lht.IP{i}(2,1) - DynamicPPTrials(Selection(q)).Lht.IP{i}(1,1))*0.3),  DynamicPPTrials(Selection(q)).Lht.IP{i}(1,2) + 2*((DynamicPPTrials(Selection(q)).Lht.IP{i}(2,2) - DynamicPPTrials(Selection(q)).Lht.IP{i}(1,2))*0.3)];
%         end
%         if Block.HC == 0
%             DynamicPPTrials(Selection(q)).LProg.HAHeelCent{i} = [DynamicPPTrials(Selection(q)).Lht.HeelCent{i}(1,1) + ((DynamicPPTrials(Selection(q)).Lht.HeelCent{i}(2,1) - DynamicPPTrials(Selection(q)).Lht.HeelCent{i}(1,1))*0.3),  DynamicPPTrials(Selection(q)).Lht.HeelCent{i}(1,2) + ((DynamicPPTrials(Selection(q)).Lht.HeelCent{i}(2,2) - DynamicPPTrials(Selection(q)).Lht.HeelCent{i}(1,2))*0.3)];
%             DynamicPPTrials(Selection(q)).LProg.FAHeelCent{i} = [DynamicPPTrials(Selection(q)).Lht.HeelCent{i}(1,1) + 2*((DynamicPPTrials(Selection(q)).Lht.HeelCent{i}(2,1) - DynamicPPTrials(Selection(q)).Lht.HeelCent{i}(1,1))*0.3),  DynamicPPTrials(Selection(q)).Lht.HeelCent{i}(1,2) + 2*((DynamicPPTrials(Selection(q)).Lht.HeelCent{i}(2,2) - DynamicPPTrials(Selection(q)).Lht.HeelCent{i}(1,2))*0.3)];
%         end
%         % define the inverse slope for the 1/3 sectioning lines
%         if strcmp(PPSettings.MaskChoice, 'Manual') || strcmp(PPSettings.MaskChoice, 'Manual (pink)') || strcmp(PPSettings.MaskChoice, 'Validation')
%             DynamicPPTrials(Selection(q)).LProg.HAlineManual{i} = [DynamicPPTrials(Selection(q)).LProg.HAManual{i} + DynamicPPTrials(Selection(q)).Lht.ManualInV{i}; DynamicPPTrials(Selection(q)).LProg.HAManual{i} - DynamicPPTrials(Selection(q)).Lht.ManualInV{i}];
%             DynamicPPTrials(Selection(q)).LProg.FAlineManual{i} = [DynamicPPTrials(Selection(q)).LProg.FAManual{i} + DynamicPPTrials(Selection(q)).Lht.ManualInV{i}; DynamicPPTrials(Selection(q)).LProg.FAManual{i} - DynamicPPTrials(Selection(q)).Lht.ManualInV{i}];
%         end
%         DynamicPPTrials(Selection(q)).LProg.HAlineGeneral{i} = [DynamicPPTrials(Selection(q)).LProg.HAGeneral{i} + DynamicPPTrials(Selection(q)).Lht.GeneralInV{i}; DynamicPPTrials(Selection(q)).LProg.HAGeneral{i} - DynamicPPTrials(Selection(q)).Lht.GeneralInV{i}];
%         DynamicPPTrials(Selection(q)).LProg.FAlineGeneral{i} = [DynamicPPTrials(Selection(q)).LProg.FAGeneral{i} + DynamicPPTrials(Selection(q)).Lht.GeneralInV{i}; DynamicPPTrials(Selection(q)).LProg.FAGeneral{i} - DynamicPPTrials(Selection(q)).Lht.GeneralInV{i}];
%         DynamicPPTrials(Selection(q)).LProg.HAlineCoP{i} = [DynamicPPTrials(Selection(q)).LProg.HACoP{i} + DynamicPPTrials(Selection(q)).Lht.InVCoP{i}; DynamicPPTrials(Selection(q)).LProg.HACoP{i} - DynamicPPTrials(Selection(q)).Lht.InVCoP{i}];
%         DynamicPPTrials(Selection(q)).LProg.FAlineCoP{i} = [DynamicPPTrials(Selection(q)).LProg.FACoP{i} + DynamicPPTrials(Selection(q)).Lht.InVCoP{i}; DynamicPPTrials(Selection(q)).LProg.FACoP{i} - DynamicPPTrials(Selection(q)).Lht.InVCoP{i}];
%         DynamicPPTrials(Selection(q)).LProg.HAline66{i} = [DynamicPPTrials(Selection(q)).LProg.HA66{i} + DynamicPPTrials(Selection(q)).Lht.CoPInV66{i}; DynamicPPTrials(Selection(q)).LProg.HA66{i} - DynamicPPTrials(Selection(q)).Lht.CoPInV66{i}];
%         DynamicPPTrials(Selection(q)).LProg.FAline66{i} = [DynamicPPTrials(Selection(q)).LProg.FA66{i} + DynamicPPTrials(Selection(q)).Lht.CoPInV66{i}; DynamicPPTrials(Selection(q)).LProg.FA66{i} - DynamicPPTrials(Selection(q)).Lht.CoPInV66{i}];
%         if Block.IP == 0
%             DynamicPPTrials(Selection(q)).LProg.HAlineIP{i} = [DynamicPPTrials(Selection(q)).LProg.HAIP{i} + DynamicPPTrials(Selection(q)).Lht.InVIP{i}; DynamicPPTrials(Selection(q)).LProg.HAIP{i} - DynamicPPTrials(Selection(q)).Lht.InVIP{i}];
%             DynamicPPTrials(Selection(q)).LProg.FAlineIP{i} = [DynamicPPTrials(Selection(q)).LProg.FAIP{i} + DynamicPPTrials(Selection(q)).Lht.InVIP{i}; DynamicPPTrials(Selection(q)).LProg.FAIP{i} - DynamicPPTrials(Selection(q)).Lht.InVIP{i}];
%         end
%         if Block.HC == 0
%             DynamicPPTrials(Selection(q)).LProg.HAlineHeelCent{i} = [DynamicPPTrials(Selection(q)).LProg.HAHeelCent{i} + DynamicPPTrials(Selection(q)).Lht.InVHeelCent{i}; DynamicPPTrials(Selection(q)).LProg.HAHeelCent{i} - DynamicPPTrials(Selection(q)).Lht.InVHeelCent{i}];
%             DynamicPPTrials(Selection(q)).LProg.FAlineHeelCent{i} = [DynamicPPTrials(Selection(q)).LProg.FAHeelCent{i} + DynamicPPTrials(Selection(q)).Lht.InVHeelCent{i}; DynamicPPTrials(Selection(q)).LProg.FAHeelCent{i} - DynamicPPTrials(Selection(q)).Lht.InVHeelCent{i}];
%         end
%         % find intersections between third sectioning lines and lat/med borders
%         if strcmp(PPSettings.MaskChoice, 'Manual') || strcmp(PPSettings.MaskChoice, 'Manual (pink)') || strcmp(PPSettings.MaskChoice, 'Validation')
%             DynamicPPTrials(Selection(q)).LProg.HALatManual{i} = linlinintersect([DynamicPPTrials(Selection(q)).LProg.HAlineManual{i};DynamicPPTrials(Selection(q)).L.LatManual{i}]);
%             DynamicPPTrials(Selection(q)).LProg.HAMedManual{i} = linlinintersect([DynamicPPTrials(Selection(q)).LProg.HAlineManual{i};DynamicPPTrials(Selection(q)).L.MedManual{i}]);
%             DynamicPPTrials(Selection(q)).LProg.FALatManual{i} = linlinintersect([DynamicPPTrials(Selection(q)).LProg.FAlineManual{i};DynamicPPTrials(Selection(q)).L.LatManual{i}]);
%             DynamicPPTrials(Selection(q)).LProg.FAMedManual{i} = linlinintersect([DynamicPPTrials(Selection(q)).LProg.FAlineManual{i};DynamicPPTrials(Selection(q)).L.MedManual{i}]);
%         end
%         DynamicPPTrials(Selection(q)).LProg.HALatGeneral{i} = linlinintersect([DynamicPPTrials(Selection(q)).LProg.HAlineGeneral{i};DynamicPPTrials(Selection(q)).L.LatGeneral{i}]);
%         DynamicPPTrials(Selection(q)).LProg.HAMedGeneral{i} = linlinintersect([DynamicPPTrials(Selection(q)).LProg.HAlineGeneral{i};DynamicPPTrials(Selection(q)).L.MedGeneral{i}]);
%         DynamicPPTrials(Selection(q)).LProg.FALatGeneral{i} = linlinintersect([DynamicPPTrials(Selection(q)).LProg.FAlineGeneral{i};DynamicPPTrials(Selection(q)).L.LatGeneral{i}]);
%         DynamicPPTrials(Selection(q)).LProg.FAMedGeneral{i} = linlinintersect([DynamicPPTrials(Selection(q)).LProg.FAlineGeneral{i};DynamicPPTrials(Selection(q)).L.MedGeneral{i}]);
%         DynamicPPTrials(Selection(q)).LProg.HALatCoP{i} = linlinintersect([DynamicPPTrials(Selection(q)).LProg.HAlineCoP{i};DynamicPPTrials(Selection(q)).L.LatCoP{i}]);
%         DynamicPPTrials(Selection(q)).LProg.HAMedCoP{i} = linlinintersect([DynamicPPTrials(Selection(q)).LProg.HAlineCoP{i};DynamicPPTrials(Selection(q)).L.MedCoP{i}]);
%         DynamicPPTrials(Selection(q)).LProg.FALatCoP{i} = linlinintersect([DynamicPPTrials(Selection(q)).LProg.FAlineCoP{i};DynamicPPTrials(Selection(q)).L.LatCoP{i}]);
%         DynamicPPTrials(Selection(q)).LProg.FAMedCoP{i} = linlinintersect([DynamicPPTrials(Selection(q)).LProg.FAlineCoP{i};DynamicPPTrials(Selection(q)).L.MedCoP{i}]);
%         DynamicPPTrials(Selection(q)).LProg.HALat66{i} = linlinintersect([DynamicPPTrials(Selection(q)).LProg.HAline66{i};DynamicPPTrials(Selection(q)).L.Lat66{i}]);
%         DynamicPPTrials(Selection(q)).LProg.HAMed66{i} = linlinintersect([DynamicPPTrials(Selection(q)).LProg.HAline66{i};DynamicPPTrials(Selection(q)).L.Med66{i}]);
%         DynamicPPTrials(Selection(q)).LProg.FALat66{i} = linlinintersect([DynamicPPTrials(Selection(q)).LProg.FAline66{i};DynamicPPTrials(Selection(q)).L.Lat66{i}]);
%         DynamicPPTrials(Selection(q)).LProg.FAMed66{i} = linlinintersect([DynamicPPTrials(Selection(q)).LProg.FAline66{i};DynamicPPTrials(Selection(q)).L.Med66{i}]);
%         if Block.IP == 0
%             DynamicPPTrials(Selection(q)).LProg.HALatIP{i} = linlinintersect([DynamicPPTrials(Selection(q)).LProg.HAlineIP{i};DynamicPPTrials(Selection(q)).L.LatIP{i}]);
%             DynamicPPTrials(Selection(q)).LProg.HAMedIP{i} = linlinintersect([DynamicPPTrials(Selection(q)).LProg.HAlineIP{i};DynamicPPTrials(Selection(q)).L.MedIP{i}]);
%             DynamicPPTrials(Selection(q)).LProg.FALatIP{i} = linlinintersect([DynamicPPTrials(Selection(q)).LProg.FAlineIP{i};DynamicPPTrials(Selection(q)).L.LatIP{i}]);
%             DynamicPPTrials(Selection(q)).LProg.FAMedIP{i} = linlinintersect([DynamicPPTrials(Selection(q)).LProg.FAlineIP{i};DynamicPPTrials(Selection(q)).L.MedIP{i}]);
%         end
%         if Block.HC == 0
%             DynamicPPTrials(Selection(q)).LProg.HALatHeelCent{i} = linlinintersect([DynamicPPTrials(Selection(q)).LProg.HAlineHeelCent{i};DynamicPPTrials(Selection(q)).L.LatHeelCent{i}]);
%             DynamicPPTrials(Selection(q)).LProg.HAMedHeelCent{i} = linlinintersect([DynamicPPTrials(Selection(q)).LProg.HAlineHeelCent{i};DynamicPPTrials(Selection(q)).L.MedHeelCent{i}]);
%             DynamicPPTrials(Selection(q)).LProg.FALatHeelCent{i} = linlinintersect([DynamicPPTrials(Selection(q)).LProg.FAlineHeelCent{i};DynamicPPTrials(Selection(q)).L.LatHeelCent{i}]);
%             DynamicPPTrials(Selection(q)).LProg.FAMedHeelCent{i} = linlinintersect([DynamicPPTrials(Selection(q)).LProg.FAlineHeelCent{i};DynamicPPTrials(Selection(q)).L.MedHeelCent{i}]);
%         end
% %         % show the third-sectioning lines
% %         subplot(1,NumDynTrials, Selection(q));
% %         line([DynamicPPTrials(Selection(q)).LProg.HALat{i}(1) DynamicPPTrials(Selection(q)).LProg.HAMed{i}(1) ],[DynamicPPTrials(Selection(q)).LProg.HALat{i}(2) DynamicPPTrials(Selection(q)).LProg.HAMed{i}(2)],'Color','k');
% %         line([DynamicPPTrials(Selection(q)).LProg.FALat{i}(1) DynamicPPTrials(Selection(q)).LProg.FAMed{i}(1)],[DynamicPPTrials(Selection(q)).LProg.FALat{i}(2) DynamicPPTrials(Selection(q)).LProg.FAMed{i}(2)],'Color','k');
% %         % plot c3d data foot coordinates if input
% %         if strcmp(PPSettings.C3DInput,'Yes')
% %             if isempty(DynamicPPTrials(Selection(q)).LProg.HA) == 0
% %                 subplot(1,NumDynTrials, (Selection(q)));
% %                 hold on;
% %                 plot(C3Ddata(Selection(q)).LHA(i,2), C3Ddata(Selection(q)).LHA(i,1),'om');
% %                 plot(C3Ddata(Selection(q)).LFA(i,2), C3Ddata(Selection(q)).LFA(i,1),'om');
% %             end
% %         end
%     end
%     % RIGHT
%     for i = 1:DynamicPPTrials(Selection(q)).NumRight
%         % find the 1/3 of foot prog line points
%         if strcmp(PPSettings.MaskChoice, 'Manual') || strcmp(PPSettings.MaskChoice, 'Manual (pink)') || strcmp(PPSettings.MaskChoice, 'Validation')
%             DynamicPPTrials(Selection(q)).RProg.HAManual{i} = [DynamicPPTrials(Selection(q)).Rht.Manual{i}(1,1) + ((DynamicPPTrials(Selection(q)).Rht.Manual{i}(2,1) - DynamicPPTrials(Selection(q)).Rht.Manual{i}(1,1))*0.3),  DynamicPPTrials(Selection(q)).Rht.Manual{i}(1,2) + ((DynamicPPTrials(Selection(q)).Rht.Manual{i}(2,2) - DynamicPPTrials(Selection(q)).Rht.Manual{i}(1,2))*0.3)];
%             DynamicPPTrials(Selection(q)).RProg.FAManual{i} = [DynamicPPTrials(Selection(q)).Rht.Manual{i}(1,1) + 2*((DynamicPPTrials(Selection(q)).Rht.Manual{i}(2,1) - DynamicPPTrials(Selection(q)).Rht.Manual{i}(1,1))*0.3),  DynamicPPTrials(Selection(q)).Rht.Manual{i}(1,2) + 2*((DynamicPPTrials(Selection(q)).Rht.Manual{i}(2,2) - DynamicPPTrials(Selection(q)).Rht.Manual{i}(1,2))*0.3)];
%         end
%         DynamicPPTrials(Selection(q)).RProg.HAGeneral{i} = [DynamicPPTrials(Selection(q)).Rht.General{i}(1,1) + ((DynamicPPTrials(Selection(q)).Rht.General{i}(2,1) - DynamicPPTrials(Selection(q)).Rht.General{i}(1,1))*0.3),  DynamicPPTrials(Selection(q)).Rht.General{i}(1,2) + ((DynamicPPTrials(Selection(q)).Rht.General{i}(2,2) - DynamicPPTrials(Selection(q)).Rht.General{i}(1,2))*0.3)];
%         DynamicPPTrials(Selection(q)).RProg.FAGeneral{i} = [DynamicPPTrials(Selection(q)).Rht.General{i}(1,1) + 2*((DynamicPPTrials(Selection(q)).Rht.General{i}(2,1) - DynamicPPTrials(Selection(q)).Rht.General{i}(1,1))*0.3),  DynamicPPTrials(Selection(q)).Rht.General{i}(1,2) + 2*((DynamicPPTrials(Selection(q)).Rht.General{i}(2,2) - DynamicPPTrials(Selection(q)).Rht.General{i}(1,2))*0.3)];
%         DynamicPPTrials(Selection(q)).RProg.HACoP{i} = [DynamicPPTrials(Selection(q)).Rht.CoP{i}(1,1) + ((DynamicPPTrials(Selection(q)).Rht.CoP{i}(2,1) - DynamicPPTrials(Selection(q)).Rht.CoP{i}(1,1))*0.3),  DynamicPPTrials(Selection(q)).Rht.CoP{i}(1,2) + ((DynamicPPTrials(Selection(q)).Rht.CoP{i}(2,2) - DynamicPPTrials(Selection(q)).Rht.CoP{i}(1,2))*0.3)];
%         DynamicPPTrials(Selection(q)).RProg.FACoP{i} = [DynamicPPTrials(Selection(q)).Rht.CoP{i}(1,1) + 2*((DynamicPPTrials(Selection(q)).Rht.CoP{i}(2,1) - DynamicPPTrials(Selection(q)).Rht.CoP{i}(1,1))*0.3),  DynamicPPTrials(Selection(q)).Rht.CoP{i}(1,2) + 2*((DynamicPPTrials(Selection(q)).Rht.CoP{i}(2,2) - DynamicPPTrials(Selection(q)).Rht.CoP{i}(1,2))*0.3)];
%         DynamicPPTrials(Selection(q)).RProg.HA66{i} = [DynamicPPTrials(Selection(q)).Rht.CoP66{i}(1,1) + ((DynamicPPTrials(Selection(q)).Rht.CoP66{i}(2,1) - DynamicPPTrials(Selection(q)).Rht.CoP66{i}(1,1))*0.3),  DynamicPPTrials(Selection(q)).Rht.CoP66{i}(1,2) + ((DynamicPPTrials(Selection(q)).Rht.CoP66{i}(2,2) - DynamicPPTrials(Selection(q)).Rht.CoP66{i}(1,2))*0.3)];
%         DynamicPPTrials(Selection(q)).RProg.FA66{i} = [DynamicPPTrials(Selection(q)).Rht.CoP66{i}(1,1) + 2*((DynamicPPTrials(Selection(q)).Rht.CoP66{i}(2,1) - DynamicPPTrials(Selection(q)).Rht.CoP66{i}(1,1))*0.3),  DynamicPPTrials(Selection(q)).Rht.CoP66{i}(1,2) + 2*((DynamicPPTrials(Selection(q)).Rht.CoP66{i}(2,2) - DynamicPPTrials(Selection(q)).Rht.CoP66{i}(1,2))*0.3)];
%         if Block.IP == 0
%             DynamicPPTrials(Selection(q)).RProg.HAIP{i} = [DynamicPPTrials(Selection(q)).Rht.IP{i}(1,1) + ((DynamicPPTrials(Selection(q)).Rht.IP{i}(2,1) - DynamicPPTrials(Selection(q)).Rht.IP{i}(1,1))*0.3),  DynamicPPTrials(Selection(q)).Rht.IP{i}(1,2) + ((DynamicPPTrials(Selection(q)).Rht.IP{i}(2,2) - DynamicPPTrials(Selection(q)).Rht.IP{i}(1,2))*0.3)];
%             DynamicPPTrials(Selection(q)).RProg.FAIP{i} = [DynamicPPTrials(Selection(q)).Rht.IP{i}(1,1) + 2*((DynamicPPTrials(Selection(q)).Rht.IP{i}(2,1) - DynamicPPTrials(Selection(q)).Rht.IP{i}(1,1))*0.3),  DynamicPPTrials(Selection(q)).Rht.IP{i}(1,2) + 2*((DynamicPPTrials(Selection(q)).Rht.IP{i}(2,2) - DynamicPPTrials(Selection(q)).Rht.IP{i}(1,2))*0.3)];
%         end
%         if Block.HC == 0
%             DynamicPPTrials(Selection(q)).RProg.HAHeelCent{i} = [DynamicPPTrials(Selection(q)).Rht.HeelCent{i}(1,1) + ((DynamicPPTrials(Selection(q)).Rht.HeelCent{i}(2,1) - DynamicPPTrials(Selection(q)).Rht.HeelCent{i}(1,1))*0.3),  DynamicPPTrials(Selection(q)).Rht.HeelCent{i}(1,2) + ((DynamicPPTrials(Selection(q)).Rht.HeelCent{i}(2,2) - DynamicPPTrials(Selection(q)).Rht.HeelCent{i}(1,2))*0.3)];
%             DynamicPPTrials(Selection(q)).RProg.FAHeelCent{i} = [DynamicPPTrials(Selection(q)).Rht.HeelCent{i}(1,1) + 2*((DynamicPPTrials(Selection(q)).Rht.HeelCent{i}(2,1) - DynamicPPTrials(Selection(q)).Rht.HeelCent{i}(1,1))*0.3),  DynamicPPTrials(Selection(q)).Rht.HeelCent{i}(1,2) + 2*((DynamicPPTrials(Selection(q)).Rht.HeelCent{i}(2,2) - DynamicPPTrials(Selection(q)).Rht.HeelCent{i}(1,2))*0.3)];
%         end
%         % define the inverse slope for the 1/3 sectioning lines
%         if strcmp(PPSettings.MaskChoice, 'Manual') || strcmp(PPSettings.MaskChoice, 'Manual (pink)') || strcmp(PPSettings.MaskChoice, 'Validation')
%             DynamicPPTrials(Selection(q)).RProg.HAlineManual{i} = [DynamicPPTrials(Selection(q)).RProg.HAManual{i} + DynamicPPTrials(Selection(q)).Rht.ManualInV{i}; DynamicPPTrials(Selection(q)).RProg.HAManual{i} - DynamicPPTrials(Selection(q)).Rht.ManualInV{i}];
%             DynamicPPTrials(Selection(q)).RProg.FAlineManual{i} = [DynamicPPTrials(Selection(q)).RProg.FAManual{i} + DynamicPPTrials(Selection(q)).Rht.ManualInV{i}; DynamicPPTrials(Selection(q)).RProg.FAManual{i} - DynamicPPTrials(Selection(q)).Rht.ManualInV{i}];
%         end
%         DynamicPPTrials(Selection(q)).RProg.HAlineGeneral{i} = [DynamicPPTrials(Selection(q)).RProg.HAGeneral{i} + DynamicPPTrials(Selection(q)).Rht.GeneralInV{i}; DynamicPPTrials(Selection(q)).RProg.HAGeneral{i} - DynamicPPTrials(Selection(q)).Rht.GeneralInV{i}];
%         DynamicPPTrials(Selection(q)).RProg.FAlineGeneral{i} = [DynamicPPTrials(Selection(q)).RProg.FAGeneral{i} + DynamicPPTrials(Selection(q)).Rht.GeneralInV{i}; DynamicPPTrials(Selection(q)).RProg.FAGeneral{i} - DynamicPPTrials(Selection(q)).Rht.GeneralInV{i}];
%         DynamicPPTrials(Selection(q)).RProg.HAlineCoP{i} = [DynamicPPTrials(Selection(q)).RProg.HACoP{i} + DynamicPPTrials(Selection(q)).Rht.InVCoP{i}; DynamicPPTrials(Selection(q)).RProg.HACoP{i} - DynamicPPTrials(Selection(q)).Rht.InVCoP{i}];
%         DynamicPPTrials(Selection(q)).RProg.FAlineCoP{i} = [DynamicPPTrials(Selection(q)).RProg.FACoP{i} + DynamicPPTrials(Selection(q)).Rht.InVCoP{i}; DynamicPPTrials(Selection(q)).RProg.FACoP{i} - DynamicPPTrials(Selection(q)).Rht.InVCoP{i}];
%         DynamicPPTrials(Selection(q)).RProg.HAline66{i} = [DynamicPPTrials(Selection(q)).RProg.HA66{i} + DynamicPPTrials(Selection(q)).Rht.CoPInV66{i}; DynamicPPTrials(Selection(q)).RProg.HA66{i} - DynamicPPTrials(Selection(q)).Rht.CoPInV66{i}];
%         DynamicPPTrials(Selection(q)).RProg.FAline66{i} = [DynamicPPTrials(Selection(q)).RProg.FA66{i} + DynamicPPTrials(Selection(q)).Rht.CoPInV66{i}; DynamicPPTrials(Selection(q)).RProg.FA66{i} - DynamicPPTrials(Selection(q)).Rht.CoPInV66{i}];
%         if Block.IP == 0
%             DynamicPPTrials(Selection(q)).RProg.HAlineIP{i} = [DynamicPPTrials(Selection(q)).RProg.HAIP{i} + DynamicPPTrials(Selection(q)).Rht.InVIP{i}; DynamicPPTrials(Selection(q)).RProg.HAIP{i} - DynamicPPTrials(Selection(q)).Rht.InVIP{i}];
%             DynamicPPTrials(Selection(q)).RProg.FAlineIP{i} = [DynamicPPTrials(Selection(q)).RProg.FAIP{i} + DynamicPPTrials(Selection(q)).Rht.InVIP{i}; DynamicPPTrials(Selection(q)).RProg.FAIP{i} - DynamicPPTrials(Selection(q)).Rht.InVIP{i}];
%         end
%         if Block.HC == 0
%             DynamicPPTrials(Selection(q)).RProg.HAlineHeelCent{i} = [DynamicPPTrials(Selection(q)).RProg.HAHeelCent{i} + DynamicPPTrials(Selection(q)).Rht.InVHeelCent{i}; DynamicPPTrials(Selection(q)).RProg.HAHeelCent{i} - DynamicPPTrials(Selection(q)).Rht.InVHeelCent{i}];
%             DynamicPPTrials(Selection(q)).RProg.FAlineHeelCent{i} = [DynamicPPTrials(Selection(q)).RProg.FAHeelCent{i} + DynamicPPTrials(Selection(q)).Rht.InVHeelCent{i}; DynamicPPTrials(Selection(q)).RProg.FAHeelCent{i} - DynamicPPTrials(Selection(q)).Rht.InVHeelCent{i}];
%         end
%         % find intersections between third sectioning lines and lat/med borders
%         if strcmp(PPSettings.MaskChoice, 'Manual') || strcmp(PPSettings.MaskChoice, 'Manual (pink)') || strcmp(PPSettings.MaskChoice, 'Validation')
%             DynamicPPTrials(Selection(q)).RProg.HALatManual{i} = linlinintersect([DynamicPPTrials(Selection(q)).RProg.HAlineManual{i};DynamicPPTrials(Selection(q)).R.LatManual{i}]);
%             DynamicPPTrials(Selection(q)).RProg.HAMedManual{i} = linlinintersect([DynamicPPTrials(Selection(q)).RProg.HAlineManual{i};DynamicPPTrials(Selection(q)).R.MedManual{i}]);
%             DynamicPPTrials(Selection(q)).RProg.FALatManual{i} = linlinintersect([DynamicPPTrials(Selection(q)).RProg.FAlineManual{i};DynamicPPTrials(Selection(q)).R.LatManual{i}]);
%             DynamicPPTrials(Selection(q)).RProg.FAMedManual{i} = linlinintersect([DynamicPPTrials(Selection(q)).RProg.FAlineManual{i};DynamicPPTrials(Selection(q)).R.MedManual{i}]);
%         end
%         DynamicPPTrials(Selection(q)).RProg.HALatGeneral{i} = linlinintersect([DynamicPPTrials(Selection(q)).RProg.HAlineGeneral{i};DynamicPPTrials(Selection(q)).R.LatGeneral{i}]);
%         DynamicPPTrials(Selection(q)).RProg.HAMedGeneral{i} = linlinintersect([DynamicPPTrials(Selection(q)).RProg.HAlineGeneral{i};DynamicPPTrials(Selection(q)).R.MedGeneral{i}]);
%         DynamicPPTrials(Selection(q)).RProg.FALatGeneral{i} = linlinintersect([DynamicPPTrials(Selection(q)).RProg.FAlineGeneral{i};DynamicPPTrials(Selection(q)).R.LatGeneral{i}]);
%         DynamicPPTrials(Selection(q)).RProg.FAMedGeneral{i} = linlinintersect([DynamicPPTrials(Selection(q)).RProg.FAlineGeneral{i};DynamicPPTrials(Selection(q)).R.MedGeneral{i}]);
%         DynamicPPTrials(Selection(q)).RProg.HALatCoP{i} = linlinintersect([DynamicPPTrials(Selection(q)).RProg.HAlineCoP{i};DynamicPPTrials(Selection(q)).R.LatCoP{i}]);
%         DynamicPPTrials(Selection(q)).RProg.HAMedCoP{i} = linlinintersect([DynamicPPTrials(Selection(q)).RProg.HAlineCoP{i};DynamicPPTrials(Selection(q)).R.MedCoP{i}]);
%         DynamicPPTrials(Selection(q)).RProg.FALatCoP{i} = linlinintersect([DynamicPPTrials(Selection(q)).RProg.FAlineCoP{i};DynamicPPTrials(Selection(q)).R.LatCoP{i}]);
%         DynamicPPTrials(Selection(q)).RProg.FAMedCoP{i} = linlinintersect([DynamicPPTrials(Selection(q)).RProg.FAlineCoP{i};DynamicPPTrials(Selection(q)).R.MedCoP{i}]);
%         DynamicPPTrials(Selection(q)).RProg.HALat66{i} = linlinintersect([DynamicPPTrials(Selection(q)).RProg.HAline66{i};DynamicPPTrials(Selection(q)).R.Lat66{i}]);
%         DynamicPPTrials(Selection(q)).RProg.HAMed66{i} = linlinintersect([DynamicPPTrials(Selection(q)).RProg.HAline66{i};DynamicPPTrials(Selection(q)).R.Med66{i}]);
%         DynamicPPTrials(Selection(q)).RProg.FALat66{i} = linlinintersect([DynamicPPTrials(Selection(q)).RProg.FAline66{i};DynamicPPTrials(Selection(q)).R.Lat66{i}]);
%         DynamicPPTrials(Selection(q)).RProg.FAMed66{i} = linlinintersect([DynamicPPTrials(Selection(q)).RProg.FAline66{i};DynamicPPTrials(Selection(q)).R.Med66{i}]);
%         if Block.IP == 0
%             DynamicPPTrials(Selection(q)).RProg.HALatIP{i} = linlinintersect([DynamicPPTrials(Selection(q)).RProg.HAlineIP{i};DynamicPPTrials(Selection(q)).R.LatIP{i}]);
%             DynamicPPTrials(Selection(q)).RProg.HAMedIP{i} = linlinintersect([DynamicPPTrials(Selection(q)).RProg.HAlineIP{i};DynamicPPTrials(Selection(q)).R.MedIP{i}]);
%             DynamicPPTrials(Selection(q)).RProg.FALatIP{i} = linlinintersect([DynamicPPTrials(Selection(q)).RProg.FAlineIP{i};DynamicPPTrials(Selection(q)).R.LatIP{i}]);
%             DynamicPPTrials(Selection(q)).RProg.FAMedIP{i} = linlinintersect([DynamicPPTrials(Selection(q)).RProg.FAlineIP{i};DynamicPPTrials(Selection(q)).R.MedIP{i}]);
%         end
%         if Block.HC == 0
%             DynamicPPTrials(Selection(q)).RProg.HALatHeelCent{i} = linlinintersect([DynamicPPTrials(Selection(q)).RProg.HAlineHeelCent{i};DynamicPPTrials(Selection(q)).R.LatHeelCent{i}]);
%             DynamicPPTrials(Selection(q)).RProg.HAMedHeelCent{i} = linlinintersect([DynamicPPTrials(Selection(q)).RProg.HAlineHeelCent{i};DynamicPPTrials(Selection(q)).R.MedHeelCent{i}]);
%             DynamicPPTrials(Selection(q)).RProg.FALatHeelCent{i} = linlinintersect([DynamicPPTrials(Selection(q)).RProg.FAlineHeelCent{i};DynamicPPTrials(Selection(q)).R.LatHeelCent{i}]);
%             DynamicPPTrials(Selection(q)).RProg.FAMedHeelCent{i} = linlinintersect([DynamicPPTrials(Selection(q)).RProg.FAlineHeelCent{i};DynamicPPTrials(Selection(q)).R.MedHeelCent{i}]);
%         end
% %         % show the third-sectioning lines
% %         subplot(1,NumDynTrials, Selection(q));
% %         line([DynamicPPTrials(Selection(q)).RProg.HALatGeneral{i}(1) DynamicPPTrials(Selection(q)).RProg.HAMedGeneral{i}(1) ],[DynamicPPTrials(Selection(q)).RProg.HALatGeneral{i}(2) DynamicPPTrials(Selection(q)).RProg.HAMedGeneral{i}(2)],'Color','k');
% %         line([DynamicPPTrials(Selection(q)).RProg.FALatGeneral{i}(1) DynamicPPTrials(Selection(q)).RProg.FAMedGeneral{i}(1)],[DynamicPPTrials(Selection(q)).RProg.FALatGeneral{i}(2) DynamicPPTrials(Selection(q)).RProg.FAMedGenearl{i}(2)],'Color','k');
% %         % plot c3d data foot coordinates if input
% %         if strcmp(PPSettings.C3DInput,'Yes')
% %             if isempty(DynamicPPTrials(Selection(q)).RProg.HA) == 0
% %                 subplot(1,NumDynTrials, (Selection(q)));
% %                 hold on;
% %                 plot(C3Ddata(Selection(q)).RHA(i,2), C3Ddata(Selection(q)).RHA(i,1),'om');
% %                 plot(C3Ddata(Selection(q)).RFA(i,2), C3Ddata(Selection(q)).RFA(i,1),'om');
% %             end
% %         end
%     end
% end
% 
% %% Plot each trial and save
% close;
% W = 1.5; % set line width
% figure('Position', [DisplayDim(1) DisplayDim(2) (DisplayDim(3)/NumDynTrials)+100 DisplayDim(4)]);
% for i = 1:NumDynTrials
%     contour(DynamicPPTrials(i).SumTM,100); hold on; 
%     axis equal;
%     title(['Dynamic Trial ' num2str(i)], 'FontSize',8);
%     % LEFT 
%     if isempty(DynamicPPTrials(i).Lht) == 0
%         for j = 1:DynamicPPTrials(i).NumLeft
%             if strcmp(PPSettings.MaskChoice, 'Manual') || strcmp(PPSettings.MaskChoice, 'Manual (pink)') || strcmp(PPSettings.MaskChoice, 'Validation')
%                 line(DynamicPPTrials(i).Lht.Manual{j}(:,1),DynamicPPTrials(i).Lht.Manual{j}(:,2), 'Color',rgb('HotPink'), 'LineWidth', W);
%             end
%             if strcmp(PPSettings.MaskChoice, 'General') || strcmp(PPSettings.MaskChoice, 'General (yellow)') || strcmp(PPSettings.MaskChoice, 'Validation')
%                 line(DynamicPPTrials(i).Lht.General{j}(:,1),DynamicPPTrials(i).Lht.General{j}(:,2), 'Color',rgb('Gold'), 'LineWidth', W);
%             end
%             if strcmp(PPSettings.MaskChoice, 'CoP') || strcmp(PPSettings.MaskChoice, 'CoP (purple)') || strcmp(PPSettings.MaskChoice, 'Validation')
%                 line(DynamicPPTrials(i).Lht.CoP{j}(:,1),DynamicPPTrials(i).Lht.CoP{j}(:,2), 'Color',rgb('Purple'), 'LineWidth', W);
%             end
%             if strcmp(PPSettings.MaskChoice, '66% CoP') || strcmp(PPSettings.MaskChoice, '66% CoP (cyan)') || strcmp(PPSettings.MaskChoice, 'Validation')
%                 line(DynamicPPTrials(i).Lht.CoP66{j}(:,1),DynamicPPTrials(i).Lht.CoP66{j}(:,2), 'Color',rgb('Cyan'), 'LineWidth', W);
%             end
%             if Block.IP == 0
%                 if strcmp(PPSettings.MaskChoice, 'Inter-Peak') || strcmp(PPSettings.MaskChoice, 'Inter-Peak (green)') || strcmp(PPSettings.MaskChoice, 'Validation')
%                     line(DynamicPPTrials(i).Lht.IP{j}(:,1),DynamicPPTrials(i).Lht.IP{j}(:,2), 'Color',rgb('ForestGreen'), 'LineWidth', W);
%                 end
%             end
%             if Block.HC == 0
%                 if strcmp(PPSettings.MaskChoice, 'Heel-Centroid') || strcmp(PPSettings.MaskChoice, 'Heel-Centroid (red)') || strcmp(PPSettings.MaskChoice, 'Validation')
%                     line(DynamicPPTrials(i).Lht.HeelCent{j}(:,1),DynamicPPTrials(i).Lht.HeelCent{j}(:,2), 'Color',rgb('Red'), 'LineWidth', W);
%                 end
%             end
%         end
%     end
%     
%     % RIGHT
%     if isempty(DynamicPPTrials(i).Rht) == 0
%         for j = 1:DynamicPPTrials(i).NumRight
%             if strcmp(PPSettings.MaskChoice, 'Manual') || strcmp(PPSettings.MaskChoice, 'Manual (pink)') || strcmp(PPSettings.MaskChoice, 'Validation')
%                 line(DynamicPPTrials(i).Rht.Manual{j}(:,1),DynamicPPTrials(i).Rht.Manual{j}(:,2), 'Color',rgb('HotPink'), 'LineWidth', W);
%             end
%             if strcmp(PPSettings.MaskChoice, 'General') || strcmp(PPSettings.MaskChoice, 'General (yellow)') || strcmp(PPSettings.MaskChoice, 'Validation')
%                 line(DynamicPPTrials(i).Rht.General{j}(:,1),DynamicPPTrials(i).Rht.General{j}(:,2), 'Color',rgb('Gold'), 'LineWidth', W);
%             end
%             if strcmp(PPSettings.MaskChoice, 'CoP') || strcmp(PPSettings.MaskChoice, 'CoP (purple)') ||  strcmp(PPSettings.MaskChoice, 'Validation')
%                 line(DynamicPPTrials(i).Rht.CoP{j}(:,1),DynamicPPTrials(i).Rht.CoP{j}(:,2), 'Color',rgb('Purple'), 'LineWidth', W);
%             end
%             if strcmp(PPSettings.MaskChoice, '66% CoP') || strcmp(PPSettings.MaskChoice, '66% CoP (cyan)') || strcmp(PPSettings.MaskChoice, 'Validation')
%                 line(DynamicPPTrials(i).Rht.CoP66{j}(:,1),DynamicPPTrials(i).Rht.CoP66{j}(:,2), 'Color',rgb('Cyan'), 'LineWidth', W);
%             end
%             if Block.IP == 0
%                 if strcmp(PPSettings.MaskChoice, 'Inter-Peak') || strcmp(PPSettings.MaskChoice, 'Inter-Peak (green)') || strcmp(PPSettings.MaskChoice, 'Validation')
%                     line(DynamicPPTrials(i).Rht.IP{j}(:,1),DynamicPPTrials(i).Rht.IP{j}(:,2), 'Color',rgb('ForestGreen'), 'LineWidth', W);
%                 end
%             end
%             if Block.HC == 0
%                 if strcmp(PPSettings.MaskChoice, 'Heel-Centroid') || strcmp(PPSettings.MaskChoice, 'Heel-Centroid (red)') || strcmp(PPSettings.MaskChoice, 'Validation')
%                     line(DynamicPPTrials(i).Rht.HeelCent{j}(:,1),DynamicPPTrials(i).Rht.HeelCent{j}(:,2), 'Color',rgb('Red'), 'LineWidth', W);
%                 end
%             end
%         end
%     end
%     % add in CoPs
%     for j = 1:length(Regions{1,i})
%         for k = 1:length(Regions{1,i}(j).StepCoP)
%             x = Regions{1,i}(j).BoundingBox(1) + Regions{1,i}(j).StepCoP(k,1);
%             y = Regions{1,i}(j).BoundingBox(2) + Regions{1,i}(j).StepCoP(k,2);
%             plot(x,y,'.k','MarkerSize', 8);
%         end
%     end
%     
%    % plot 3D markers
%    if strcmp(PPSettings.C3DInput, 'Yes')
%        if isempty(C3Ddata(i).LHA) == 0
%            for j = 1:DynamicPPTrials(i).NumLeft
%                % plot marker points
%                plot(C3Ddata(i).LHEE.base(j,2), C3Ddata(i).LHEE.base(j,1), '*k');
%                plot(C3Ddata(i).LTOE.top(j,2), C3Ddata(i).LTOE.top(j,1), 'ok');
%                plot(C3Ddata(i).LTOE.step(j,2), C3Ddata(i).LTOE.step(j,1), '*k');
%                %                     plot(C3Ddata(i).LHAb(j,2), C3Ddata(i).LHAb(j,1), 'ok');
%                %                     plot(C3Ddata(i).LFAb(j,2), C3Ddata(i).LFAb(j,1), 'ok');
%                plot(C3Ddata(i).LLCA.step(j,2), C3Ddata(i).LLCA.step(j,1), '*k');
%                plot(C3Ddata(i).LMCA.step(j,2), C3Ddata(i).LMCA.step(j,1), '*k');
%                plot(C3Ddata(i).LD1M.step(j,2), C3Ddata(i).LD1M.step(j,1), '*k');
%                plot(C3Ddata(i).LD5M.step(j,2), C3Ddata(i).LD5M.step(j,1), '*k');
%                % plot trisection lines
%                line([C3Ddata(i).LHEE.base(j,2); C3Ddata(i).LTOE.top(j,2)], [C3Ddata(i).LHEE.base(j,1); C3Ddata(i).LTOE.top(j,1)], 'Color', 'k','LineWidth',W); % foot prog angle
%                pt = [C3Ddata(i).LHA(j,2) + C3Ddata(i).LeftContact(j).UnitVector(1), C3Ddata(i).LHA(j,1) - C3Ddata(i).LeftContact(j).UnitVector(2)];
%                HALat = linlinintersect([pt(1), pt(2); [C3Ddata(i).LHA(j,2), C3Ddata(i).LHA(j,1)]; DynamicPPTrials(i).L.LatGeneral{j}(1,:); DynamicPPTrials(i).L.LatGeneral{j}(2,:)]);
%                HAMed = linlinintersect([pt(1), pt(2); [C3Ddata(i).LHA(j,2), C3Ddata(i).LHA(j,1)]; DynamicPPTrials(i).L.MedGeneral{j}(1,:); DynamicPPTrials(i).L.MedGeneral{j}(2,:)]);
%                line([HALat(1); HAMed(1)], [HALat(2); HAMed(2)], 'Color', 'k','LineWidth',W);
%                pt = [C3Ddata(i).LFA(j,2) + C3Ddata(i).LeftContact(j).UnitVector(1), C3Ddata(i).LFA(j,1) - C3Ddata(i).LeftContact(j).UnitVector(2)];
%                FALat = linlinintersect([pt(1), pt(2); [C3Ddata(i).LFA(j,2), C3Ddata(i).LFA(j,1)]; DynamicPPTrials(i).L.LatGeneral{j}(1,:); DynamicPPTrials(i).L.LatGeneral{j}(2,:)]);
%                FAMed = linlinintersect([pt(1), pt(2); [C3Ddata(i).LFA(j,2), C3Ddata(i).LFA(j,1)]; DynamicPPTrials(i).L.MedGeneral{j}(1,:); DynamicPPTrials(i).L.MedGeneral{j}(2,:)]);
%                line([FALat(1); FAMed(1)], [FALat(2); FAMed(2)], 'Color', 'k','LineWidth',W);
%                clearvars HALat FALat HA Med FAMed
%            end
%        end
%        if isempty(C3Ddata(i).RHA) == 0
%            for j = 1:DynamicPPTrials(i).NumRight
%                % plot marker points
%                plot(C3Ddata(i).RHEE.base(j,2), C3Ddata(i).RHEE.base(j,1), '*k');
%                plot(C3Ddata(i).RTOE.top(j,2), C3Ddata(i).RTOE.top(j,1), 'ok');
%                plot(C3Ddata(i).RTOE.step(j,2), C3Ddata(i).RTOE.step(j,1), '*k');
%                %                     plot(C3Ddata(i).RHAb(j,2), C3Ddata(i).RHAb(j,1), 'ok');
%                %                     plot(C3Ddata(i).RFAb(j,2), C3Ddata(i).RFAb(j,1), 'ok');
%                plot(C3Ddata(i).RLCA.step(j,2), C3Ddata(i).RLCA.step(j,1), '*k');
%                plot(C3Ddata(i).RMCA.step(j,2), C3Ddata(i).RMCA.step(j,1), '*k');
%                plot(C3Ddata(i).RD1M.step(j,2), C3Ddata(i).RD1M.step(j,1), '*k');
%                plot(C3Ddata(i).RD5M.step(j,2), C3Ddata(i).RD5M.step(j,1), '*k');
%                % plot trisection lines
%                line([C3Ddata(i).RHEE.base(j,2); C3Ddata(i).RTOE.top(j,2)], [C3Ddata(i).RHEE.base(j,1); C3Ddata(i).RTOE.top(j,1)], 'Color', 'k','LineWidth',W); % foot prog angle
%                pt = [C3Ddata(i).RHA(j,2) + C3Ddata(i).RightContact(j).UnitVector(1), C3Ddata(i).RHA(j,1) - C3Ddata(i).RightContact(j).UnitVector(2)];
%                HALat = linlinintersect([pt(1), pt(2); [C3Ddata(i).RHA(j,2), C3Ddata(i).RHA(j,1)]; DynamicPPTrials(i).R.LatGeneral{j}(1,:); DynamicPPTrials(i).R.LatGeneral{j}(2,:)]);
%                HAMed = linlinintersect([pt(1), pt(2); [C3Ddata(i).RHA(j,2), C3Ddata(i).RHA(j,1)]; DynamicPPTrials(i).R.MedGeneral{j}(1,:); DynamicPPTrials(i).R.MedGeneral{j}(2,:)]);
%                line([HALat(1); HAMed(1)], [HALat(2); HAMed(2)], 'Color', 'k','LineWidth',W);
%                pt = [C3Ddata(i).RFA(j,2) + C3Ddata(i).RightContact(j).UnitVector(1), C3Ddata(i).RFA(j,1) - C3Ddata(i).RightContact(j).UnitVector(2)];
%                FALat = linlinintersect([pt(1), pt(2); [C3Ddata(i).RFA(j,2), C3Ddata(i).RFA(j,1)]; DynamicPPTrials(i).R.LatGeneral{j}(1,:); DynamicPPTrials(i).R.LatGeneral{j}(2,:)]);
%                FAMed = linlinintersect([pt(1), pt(2); [C3Ddata(i).RFA(j,2), C3Ddata(i).RFA(j,1)]; DynamicPPTrials(i).R.MedGeneral{j}(1,:); DynamicPPTrials(i).R.MedGeneral{j}(2,:)]);
%                line([FALat(1); FAMed(1)], [FALat(2); FAMed(2)], 'Color', 'k','LineWidth',W);
%                clearvars HALat FALat HA Med FAMed
%            end
%        end
%    end
%     
%     set(gca, 'XTick', []);
%     set(gca, 'YTick', []);
%     
%     savefig(strcat(folder,'\',filename(1:end-4),'_Trial_',num2str(i),'.fig'));
%     clf; 
% end
% clearvars Press2Delete i j k DblEdit DblCheck h pos
% close all

% %% Determine masking type if still undetermined
% clc; disp('Masking pressures');
% while Block.IP == 1 && strcmp(PPSettings.MaskChoice, 'Inter-Peak')
%     MaskList = {'General (yellow)','CoP (purple)','66% CoP (cyan)','Inter-Peak (green)','Heel-Centroid (red)','Manual (pink)'};
%     A = listdlg( 'PromptString','Choose a different mask','ListString', MaskList, 'SelectionMode','single');
%     PPSettings.MaskChoice = MaskList{A};
% end
% 
% while Block.HC == 1 && strcmp(PPSettings.MaskChoice, 'Heel-Centroid')
%     MaskList = {'General (yellow)','CoP (purple)','66% CoP (cyan)','Inter-Peak (green)','Heel-Centroid (red)','Manual (pink)'};
%     A = listdlg( 'PromptString','Choose a different mask','ListString', MaskList, 'SelectionMode','single');
%     PPSettings.MaskChoice = MaskList{A};
% end
% 
% if strcmp(PPSettings.MaskChoice, 'Validation')
%     MaskList = {'General (yellow)','CoP (purple)','66% CoP (cyan)','Inter-Peak (green)','Heel-Centroid (red)','Manual (pink)'};
%     A = listdlg( 'PromptString','Choose a different mask','ListString', MaskList, 'SelectionMode','single');
%     PPSettings.MaskChoice = MaskList{A};
% end

%% Apply mask and save points
% Assign intersection points for zone classification
% Split up the LEFT Foot into the 6 zones
for q = 1:length(Selection)
    for i = 1:DynamicPPTrials(Selection(q)).NumLeft
        if strcmp(PPSettings.MaskChoice, 'Manual (pink)') || strcmp(PPSettings.MaskChoice, 'Manual') || strcmp(PPSettings.MaskChoice, 'Manual (none)')
            % Lateral and Medial Heel
            DynamicPPTrials(Selection(q)).Zones.L(i).LatHeel = round([DynamicPPTrials(Selection(q)).L.LatManual{i}(1,:); DynamicPPTrials(Selection(q)).Lht.Manual{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HAManual{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HALatManual{i}(1,:)]);
            DynamicPPTrials(Selection(q)).Zones.L(i).MedHeel = round([DynamicPPTrials(Selection(q)).L.MedManual{i}(1,:); DynamicPPTrials(Selection(q)).Lht.Manual{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HAManual{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HAMedManual{i}(1,:)]);
            % Lateral and Medial Arch
            DynamicPPTrials(Selection(q)).Zones.L(i).LatArch = round([DynamicPPTrials(Selection(q)).LProg.HALatManual{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HAManual{i}(1,:); DynamicPPTrials(Selection(q)).LProg.FAManual{i}(1,:); DynamicPPTrials(Selection(q)).LProg.FALatManual{i}(1,:)]);
            DynamicPPTrials(Selection(q)).Zones.L(i).MedArch = round([DynamicPPTrials(Selection(q)).LProg.HAMedManual{i}(1,:); DynamicPPTrials(Selection(q)).LProg.FAMedManual{i}(1,:); DynamicPPTrials(Selection(q)).LProg.FAManual{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HAManual{i}(1,:)]);
            % Lateral and Medial Forefoot
            DynamicPPTrials(Selection(q)).Zones.L(i).LatFore = round([DynamicPPTrials(Selection(q)).LProg.FAManual{i}(1,:); DynamicPPTrials(Selection(q)).Lht.Manual{i}(2,:); DynamicPPTrials(Selection(q)).L.LatManual{i}(2,:); DynamicPPTrials(Selection(q)).LProg.FALatManual{i}(1,:)]);
            DynamicPPTrials(Selection(q)).Zones.L(i).MedFore = round([DynamicPPTrials(Selection(q)).LProg.FAMedManual{i}(1,:); DynamicPPTrials(Selection(q)).L.MedManual{i}(2,:); DynamicPPTrials(Selection(q)).Lht.Manual{i}(2,:); DynamicPPTrials(Selection(q)).LProg.FAManual{i}(1,:)]);  
            % Whole Foot
             DynamicPPTrials(Selection(q)).Zones.L(i).Whole = round([DynamicPPTrials(Selection(q)).L.LatManual{i}(1,:); DynamicPPTrials(Selection(q)).L.LatManual{i}(2,:); DynamicPPTrials(Selection(q)).L.MedManual{i}(2,:); DynamicPPTrials(Selection(q)).L.MedManual{i}(1,:);]);
            % Prog Angle
            DynamicPPTrials(Selection(q)).LProg.Ang{i} = DynamicPPTrials(Selection(q)).LProg.AngManual{i};
        elseif strcmp(PPSettings.MaskChoice, 'General (yellow)') || strcmp(PPSettings.MaskChoice, 'General') 
            % Lateral and Medial Heel
            DynamicPPTrials(Selection(q)).Zones.L(i).LatHeel = round([DynamicPPTrials(Selection(q)).L.LatGeneral{i}(1,:); DynamicPPTrials(Selection(q)).Lht.General{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HAGeneral{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HALatGeneral{i}(1,:)]);
            DynamicPPTrials(Selection(q)).Zones.L(i).MedHeel = round([DynamicPPTrials(Selection(q)).L.MedGeneral{i}(1,:); DynamicPPTrials(Selection(q)).Lht.General{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HAGeneral{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HAMedGeneral{i}(1,:)]);
            % Lateral and Medial Arch
            DynamicPPTrials(Selection(q)).Zones.L(i).LatArch = round([DynamicPPTrials(Selection(q)).LProg.HALatGeneral{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HAGeneral{i}(1,:); DynamicPPTrials(Selection(q)).LProg.FAGeneral{i}(1,:); DynamicPPTrials(Selection(q)).LProg.FALatGeneral{i}(1,:)]);
            DynamicPPTrials(Selection(q)).Zones.L(i).MedArch = round([DynamicPPTrials(Selection(q)).LProg.HAMedGeneral{i}(1,:); DynamicPPTrials(Selection(q)).LProg.FAMedGeneral{i}(1,:); DynamicPPTrials(Selection(q)).LProg.FAGeneral{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HAGeneral{i}(1,:)]);
            % Lateral and Medial Forefoot
            DynamicPPTrials(Selection(q)).Zones.L(i).LatFore = round([DynamicPPTrials(Selection(q)).LProg.FAGeneral{i}(1,:); DynamicPPTrials(Selection(q)).Lht.General{i}(2,:); DynamicPPTrials(Selection(q)).L.LatGeneral{i}(2,:); DynamicPPTrials(Selection(q)).LProg.FALatGeneral{i}(1,:)]);
            DynamicPPTrials(Selection(q)).Zones.L(i).MedFore = round([DynamicPPTrials(Selection(q)).LProg.FAMedGeneral{i}(1,:); DynamicPPTrials(Selection(q)).L.MedGeneral{i}(2,:); DynamicPPTrials(Selection(q)).Lht.General{i}(2,:); DynamicPPTrials(Selection(q)).LProg.FAGeneral{i}(1,:)]);
            % Whole Foot
            DynamicPPTrials(Selection(q)).Zones.L(i).Whole = round([DynamicPPTrials(Selection(q)).L.LatGeneral{i}(1,:); DynamicPPTrials(Selection(q)).L.LatGeneral{i}(2,:); DynamicPPTrials(Selection(q)).L.MedGeneral{i}(2,:); DynamicPPTrials(Selection(q)).L.MedGeneral{i}(1,:);]);
            % Prog Angle
            DynamicPPTrials(Selection(q)).LProg.Ang{i} = DynamicPPTrials(Selection(q)).LProg.AngGeneral{i};
        elseif strcmp(PPSettings.MaskChoice, 'CoP (purple)') || strcmp(PPSettings.MaskChoice, 'CoP')
            % Lateral and Medial Heel
            DynamicPPTrials(Selection(q)).Zones.L(i).LatHeel = round([DynamicPPTrials(Selection(q)).L.LatCoP{i}(1,:); DynamicPPTrials(Selection(q)).Lht.CoP{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HACoP{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HALatCoP{i}(1,:)]);
            DynamicPPTrials(Selection(q)).Zones.L(i).MedHeel = round([DynamicPPTrials(Selection(q)).L.MedCoP{i}(1,:); DynamicPPTrials(Selection(q)).Lht.CoP{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HACoP{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HAMedCoP{i}(1,:)]);
            % Lateral and Medial Arch
            DynamicPPTrials(Selection(q)).Zones.L(i).LatArch = round([DynamicPPTrials(Selection(q)).LProg.HALatCoP{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HACoP{i}(1,:); DynamicPPTrials(Selection(q)).LProg.FACoP{i}(1,:); DynamicPPTrials(Selection(q)).LProg.FALatCoP{i}(1,:)]);
            DynamicPPTrials(Selection(q)).Zones.L(i).MedArch = round([DynamicPPTrials(Selection(q)).LProg.HAMedCoP{i}(1,:); DynamicPPTrials(Selection(q)).LProg.FAMedCoP{i}(1,:); DynamicPPTrials(Selection(q)).LProg.FACoP{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HACoP{i}(1,:)]);
            % Lateral and Medial Forefoot
            DynamicPPTrials(Selection(q)).Zones.L(i).LatFore = round([DynamicPPTrials(Selection(q)).LProg.FACoP{i}(1,:); DynamicPPTrials(Selection(q)).Lht.CoP{i}(2,:); DynamicPPTrials(Selection(q)).L.LatCoP{i}(2,:); DynamicPPTrials(Selection(q)).LProg.FALatCoP{i}(1,:)]);
            DynamicPPTrials(Selection(q)).Zones.L(i).MedFore = round([DynamicPPTrials(Selection(q)).LProg.FAMedCoP{i}(1,:); DynamicPPTrials(Selection(q)).L.MedCoP{i}(2,:); DynamicPPTrials(Selection(q)).Lht.CoP{i}(2,:); DynamicPPTrials(Selection(q)).LProg.FACoP{i}(1,:)]);
             % Whole Foot
             DynamicPPTrials(Selection(q)).Zones.L(i).Whole = round([DynamicPPTrials(Selection(q)).L.LatCoP{i}(1,:); DynamicPPTrials(Selection(q)).L.LatCoP{i}(2,:); DynamicPPTrials(Selection(q)).L.MedCoP{i}(2,:); DynamicPPTrials(Selection(q)).L.MedCoP{i}(1,:);]);
            % Prog Angle
            DynamicPPTrials(Selection(q)).LProg.Ang{i} = DynamicPPTrials(Selection(q)).LProg.AngCoP{i};
        elseif strcmp(PPSettings.MaskChoice, '66% CoP (cyan)') || strcmp(PPSettings.MaskChoice, '66% CoP')
            % Lateral and Medial Heel
            DynamicPPTrials(Selection(q)).Zones.L(i).LatHeel = round([DynamicPPTrials(Selection(q)).L.Lat66{i}(1,:); DynamicPPTrials(Selection(q)).Lht.CoP66{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HA66{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HALat66{i}(1,:)]);
            DynamicPPTrials(Selection(q)).Zones.L(i).MedHeel = round([DynamicPPTrials(Selection(q)).L.Med66{i}(1,:); DynamicPPTrials(Selection(q)).Lht.CoP66{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HA66{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HAMed66{i}(1,:)]);
            % Lateral and Medial Arch
            DynamicPPTrials(Selection(q)).Zones.L(i).LatArch = round([DynamicPPTrials(Selection(q)).LProg.HALat66{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HA66{i}(1,:); DynamicPPTrials(Selection(q)).LProg.FA66{i}(1,:); DynamicPPTrials(Selection(q)).LProg.FALat66{i}(1,:)]);
            DynamicPPTrials(Selection(q)).Zones.L(i).MedArch = round([DynamicPPTrials(Selection(q)).LProg.HAMed66{i}(1,:); DynamicPPTrials(Selection(q)).LProg.FAMed66{i}(1,:); DynamicPPTrials(Selection(q)).LProg.FA66{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HA66{i}(1,:)]);
            % Lateral and Medial Forefoot
            DynamicPPTrials(Selection(q)).Zones.L(i).LatFore = round([DynamicPPTrials(Selection(q)).LProg.FA66{i}(1,:); DynamicPPTrials(Selection(q)).Lht.CoP66{i}(2,:); DynamicPPTrials(Selection(q)).L.Lat66{i}(2,:); DynamicPPTrials(Selection(q)).LProg.FALat66{i}(1,:)]);
            DynamicPPTrials(Selection(q)).Zones.L(i).MedFore = round([DynamicPPTrials(Selection(q)).LProg.FAMed66{i}(1,:); DynamicPPTrials(Selection(q)).L.Med66{i}(2,:); DynamicPPTrials(Selection(q)).Lht.CoP66{i}(2,:); DynamicPPTrials(Selection(q)).LProg.FA66{i}(1,:)]);
             % Whole Foot
             DynamicPPTrials(Selection(q)).Zones.L(i).Whole = round([DynamicPPTrials(Selection(q)).L.LatCoP66{i}(1,:); DynamicPPTrials(Selection(q)).L.LatCoP66{i}(2,:); DynamicPPTrials(Selection(q)).L.MedCoP66{i}(2,:); DynamicPPTrials(Selection(q)).L.MedCoP66{i}(1,:);]);
            % Prog Angle
            DynamicPPTrials(Selection(q)).LProg.Ang{i} = DynamicPPTrials(Selection(q)).LProg.Ang66{i};
        elseif strcmp(PPSettings.MaskChoice, 'Inter-Peak (green)') || strcmp(PPSettings.MaskChoice, 'Inter-Peak')
            if Block.IP == 0
                % Lateral and Medial Heel
                DynamicPPTrials(Selection(q)).Zones.L(i).LatHeel = round([DynamicPPTrials(Selection(q)).L.LatIP{i}(1,:); DynamicPPTrials(Selection(q)).Lht.IP{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HAIP{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HALatIP{i}(1,:)]);
                DynamicPPTrials(Selection(q)).Zones.L(i).MedHeel = round([DynamicPPTrials(Selection(q)).L.MedIP{i}(1,:); DynamicPPTrials(Selection(q)).Lht.IP{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HAIP{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HAMedIP{i}(1,:)]);
                % Lateral and Medial Arch
                DynamicPPTrials(Selection(q)).Zones.L(i).LatArch = round([DynamicPPTrials(Selection(q)).LProg.HALatIP{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HAIP{i}(1,:); DynamicPPTrials(Selection(q)).LProg.FAIP{i}(1,:); DynamicPPTrials(Selection(q)).LProg.FALatIP{i}(1,:)]);
                DynamicPPTrials(Selection(q)).Zones.L(i).MedArch = round([DynamicPPTrials(Selection(q)).LProg.HAMedIP{i}(1,:); DynamicPPTrials(Selection(q)).LProg.FAMedIP{i}(1,:); DynamicPPTrials(Selection(q)).LProg.FAIP{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HAIP{i}(1,:)]);
                % Lateral and Medial Forefoot
                DynamicPPTrials(Selection(q)).Zones.L(i).LatFore = round([DynamicPPTrials(Selection(q)).LProg.FAIP{i}(1,:); DynamicPPTrials(Selection(q)).Lht.IP{i}(2,:); DynamicPPTrials(Selection(q)).L.LatIP{i}(2,:); DynamicPPTrials(Selection(q)).LProg.FALatIP{i}(1,:)]);
                DynamicPPTrials(Selection(q)).Zones.L(i).MedFore = round([DynamicPPTrials(Selection(q)).LProg.FAMedIP{i}(1,:); DynamicPPTrials(Selection(q)).L.MedIP{i}(2,:); DynamicPPTrials(Selection(q)).Lht.IP{i}(2,:); DynamicPPTrials(Selection(q)).LProg.FAIP{i}(1,:)]);
                % Whole Foot
                DynamicPPTrials(Selection(q)).Zones.L(i).Whole = round([DynamicPPTrials(Selection(q)).L.LatIP{i}(1,:); DynamicPPTrials(Selection(q)).L.LatIP{i}(2,:); DynamicPPTrials(Selection(q)).L.MedIP{i}(2,:); DynamicPPTrials(Selection(q)).L.MedIP{i}(1,:);]);
                % Prog Angle
                DynamicPPTrials(Selection(q)).LProg.Ang{i} = DynamicPPTrials(Selection(q)).LProg.AngIP{i};
            end
        elseif strcmp(PPSettings.MaskChoice, 'Heel-Centroid (red)') || strcmp(PPSettings.MaskChoice, 'Heel-Centroid')
            if Block.HC == 0
                % Lateral and Medial Heel
                DynamicPPTrials(Selection(q)).Zones.L(i).LatHeel = round([DynamicPPTrials(Selection(q)).L.LatHeelCent{i}(1,:); DynamicPPTrials(Selection(q)).Lht.HeelCent{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HAHeelCent{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HALatHeelCent{i}(1,:)]);
                DynamicPPTrials(Selection(q)).Zones.L(i).MedHeel = round([DynamicPPTrials(Selection(q)).L.MedHeelCent{i}(1,:); DynamicPPTrials(Selection(q)).Lht.HeelCent{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HAHeelCent{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HAMedHeelCent{i}(1,:)]);
                % Lateral and Medial Arch
                DynamicPPTrials(Selection(q)).Zones.L(i).LatArch = round([DynamicPPTrials(Selection(q)).LProg.HALatHeelCent{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HAHeelCent{i}(1,:); DynamicPPTrials(Selection(q)).LProg.FAHeelCent{i}(1,:); DynamicPPTrials(Selection(q)).LProg.FALatHeelCent{i}(1,:)]);
                DynamicPPTrials(Selection(q)).Zones.L(i).MedArch = round([DynamicPPTrials(Selection(q)).LProg.HAMedHeelCent{i}(1,:); DynamicPPTrials(Selection(q)).LProg.FAMedHeelCent{i}(1,:); DynamicPPTrials(Selection(q)).LProg.FAHeelCent{i}(1,:); DynamicPPTrials(Selection(q)).LProg.HAHeelCent{i}(1,:)]);
                % Lateral and Medial Forefoot
                DynamicPPTrials(Selection(q)).Zones.L(i).LatFore = round([DynamicPPTrials(Selection(q)).LProg.FAHeelCent{i}(1,:); DynamicPPTrials(Selection(q)).Lht.HeelCent{i}(2,:); DynamicPPTrials(Selection(q)).L.LatHeelCent{i}(2,:); DynamicPPTrials(Selection(q)).LProg.FALatHeelCent{i}(1,:)]);
                DynamicPPTrials(Selection(q)).Zones.L(i).MedFore = round([DynamicPPTrials(Selection(q)).LProg.FAMedHeelCent{i}(1,:); DynamicPPTrials(Selection(q)).L.MedHeelCent{i}(2,:); DynamicPPTrials(Selection(q)).Lht.HeelCent{i}(2,:); DynamicPPTrials(Selection(q)).LProg.FAHeelCent{i}(1,:)]);
                % Whole Foot
                DynamicPPTrials(Selection(q)).Zones.L(i).Whole = round([DynamicPPTrials(Selection(q)).L.LatHeelCent{i}(1,:); DynamicPPTrials(Selection(q)).L.LatHeelCent{i}(2,:); DynamicPPTrials(Selection(q)).L.MedHeelCent{i}(2,:); DynamicPPTrials(Selection(q)).L.MedHeelCent{i}(1,:);]);
                % Prog Angle
                DynamicPPTrials(Selection(q)).LProg.Ang{i} = DynamicPPTrials(Selection(q)).LProg.AngHeelCent{i};
            end
        end
    end
    
    % Split up the RIGHT Foot into the 6 zones
    for i = 1:DynamicPPTrials(Selection(q)).NumRight
        if strcmp(PPSettings.MaskChoice, 'Manual (pink)') || strcmp(PPSettings.MaskChoice, 'Manual')  || strcmp(PPSettings.MaskChoice, 'Manual (none)')
             % Lateral and Medial Heel
            DynamicPPTrials(Selection(q)).Zones.R(i).LatHeel = round([DynamicPPTrials(Selection(q)).R.LatManual{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HALatManual{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HAManual{i}(1,:); DynamicPPTrials(Selection(q)).Rht.Manual{i}(1,:)]);
            DynamicPPTrials(Selection(q)).Zones.R(i).MedHeel = round([DynamicPPTrials(Selection(q)).R.MedManual{i}(1,:); DynamicPPTrials(Selection(q)).Rht.Manual{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HAManual{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HAMedManual{i}(1,:)]);
            % Lateral and Medial Arch
            DynamicPPTrials(Selection(q)).Zones.R(i).LatArch = round([DynamicPPTrials(Selection(q)).RProg.HALatManual{i}(1,:); DynamicPPTrials(Selection(q)).RProg.FALatManual{i}(1,:); DynamicPPTrials(Selection(q)).RProg.FAManual{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HAManual{i}(1,:)]);
            DynamicPPTrials(Selection(q)).Zones.R(i).MedArch = round([DynamicPPTrials(Selection(q)).RProg.HAManual{i}(1,:); DynamicPPTrials(Selection(q)).RProg.FAManual{i}(1,:); DynamicPPTrials(Selection(q)).RProg.FAMedManual{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HAMedManual{i}(1,:)]);
            % Lateral and Medial Forefoot
            DynamicPPTrials(Selection(q)).Zones.R(i).LatFore = round([DynamicPPTrials(Selection(q)).RProg.FALatManual{i}(1,:); DynamicPPTrials(Selection(q)).R.LatManual{i}(2,:); DynamicPPTrials(Selection(q)).Rht.Manual{i}(2,:); DynamicPPTrials(Selection(q)).RProg.FAManual{i}(1,:)]);
            DynamicPPTrials(Selection(q)).Zones.R(i).MedFore = round([DynamicPPTrials(Selection(q)).RProg.FAMedManual{i}(1,:);  DynamicPPTrials(Selection(q)).RProg.FAManual{i}(1,:); DynamicPPTrials(Selection(q)).Rht.Manual{i}(2,:); DynamicPPTrials(Selection(q)).R.MedManual{i}(2,:)]);
            % Whole Foot
            DynamicPPTrials(Selection(q)).Zones.R(i).Whole = round([DynamicPPTrials(Selection(q)).R.LatManual{i}(1,:); DynamicPPTrials(Selection(q)).R.LatManual{i}(2,:); DynamicPPTrials(Selection(q)).R.MedManual{i}(2,:); DynamicPPTrials(Selection(q)).R.MedManual{i}(1,:);]);
            % Prog Angle
            DynamicPPTrials(Selection(q)).RProg.Ang{i} = DynamicPPTrials(Selection(q)).RProg.AngManual{i};
        elseif strcmp(PPSettings.MaskChoice, 'General (yellow)') ||  strcmp(PPSettings.MaskChoice, 'General')
            % Lateral and Medial Heel
            DynamicPPTrials(Selection(q)).Zones.R(i).LatHeel = round([DynamicPPTrials(Selection(q)).R.LatGeneral{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HALatGeneral{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HAGeneral{i}(1,:); DynamicPPTrials(Selection(q)).Rht.General{i}(1,:)]);
            DynamicPPTrials(Selection(q)).Zones.R(i).MedHeel = round([DynamicPPTrials(Selection(q)).R.MedGeneral{i}(1,:); DynamicPPTrials(Selection(q)).Rht.General{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HAGeneral{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HAMedGeneral{i}(1,:)]);
            % Lateral and Medial Arch
            DynamicPPTrials(Selection(q)).Zones.R(i).LatArch = round([DynamicPPTrials(Selection(q)).RProg.HALatGeneral{i}(1,:); DynamicPPTrials(Selection(q)).RProg.FALatGeneral{i}(1,:); DynamicPPTrials(Selection(q)).RProg.FAGeneral{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HAGeneral{i}(1,:)]);
            DynamicPPTrials(Selection(q)).Zones.R(i).MedArch = round([DynamicPPTrials(Selection(q)).RProg.HAGeneral{i}(1,:); DynamicPPTrials(Selection(q)).RProg.FAGeneral{i}(1,:); DynamicPPTrials(Selection(q)).RProg.FAMedGeneral{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HAMedGeneral{i}(1,:)]);
            % Lateral and Medial Forefoot
            DynamicPPTrials(Selection(q)).Zones.R(i).LatFore = round([DynamicPPTrials(Selection(q)).RProg.FALatGeneral{i}(1,:); DynamicPPTrials(Selection(q)).R.LatGeneral{i}(2,:); DynamicPPTrials(Selection(q)).Rht.General{i}(2,:); DynamicPPTrials(Selection(q)).RProg.FAGeneral{i}(1,:)]);
            DynamicPPTrials(Selection(q)).Zones.R(i).MedFore = round([DynamicPPTrials(Selection(q)).RProg.FAMedGeneral{i}(1,:);  DynamicPPTrials(Selection(q)).RProg.FAGeneral{i}(1,:); DynamicPPTrials(Selection(q)).Rht.General{i}(2,:); DynamicPPTrials(Selection(q)).R.MedGeneral{i}(2,:)]);
             % Whole Foot
            DynamicPPTrials(Selection(q)).Zones.R(i).Whole = round([DynamicPPTrials(Selection(q)).R.LatGeneral{i}(1,:); DynamicPPTrials(Selection(q)).R.LatGeneral{i}(2,:); DynamicPPTrials(Selection(q)).R.MedGeneral{i}(2,:); DynamicPPTrials(Selection(q)).R.MedGeneral{i}(1,:);]);
            % Prog Angle
            DynamicPPTrials(Selection(q)).RProg.Ang{i} = DynamicPPTrials(Selection(q)).RProg.AngGeneral{i};
        elseif strcmp(PPSettings.MaskChoice, 'CoP (purple)')  || strcmp(PPSettings.MaskChoice, 'CoP')
            % Lateral and Medial Heel
            DynamicPPTrials(Selection(q)).Zones.R(i).LatHeel = round([DynamicPPTrials(Selection(q)).R.LatCoP{i}(1,:); DynamicPPTrials(Selection(q)).Rht.CoP{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HACoP{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HALatCoP{i}(1,:)]);
            DynamicPPTrials(Selection(q)).Zones.R(i).MedHeel = round([DynamicPPTrials(Selection(q)).R.MedCoP{i}(1,:); DynamicPPTrials(Selection(q)).Rht.CoP{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HACoP{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HAMedCoP{i}(1,:)]);
            % Lateral and Medial Arch
            DynamicPPTrials(Selection(q)).Zones.R(i).LatArch = round([DynamicPPTrials(Selection(q)).RProg.HALatCoP{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HACoP{i}(1,:); DynamicPPTrials(Selection(q)).RProg.FACoP{i}(1,:); DynamicPPTrials(Selection(q)).RProg.FALatCoP{i}(1,:)]);
            DynamicPPTrials(Selection(q)).Zones.R(i).MedArch = round([DynamicPPTrials(Selection(q)).RProg.HAMedCoP{i}(1,:); DynamicPPTrials(Selection(q)).RProg.FAMedCoP{i}(1,:); DynamicPPTrials(Selection(q)).RProg.FACoP{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HACoP{i}(1,:)]);
            % Lateral and Medial Forefoot
            DynamicPPTrials(Selection(q)).Zones.R(i).LatFore = round([DynamicPPTrials(Selection(q)).RProg.FACoP{i}(1,:); DynamicPPTrials(Selection(q)).Rht.CoP{i}(2,:); DynamicPPTrials(Selection(q)).R.LatCoP{i}(2,:); DynamicPPTrials(Selection(q)).RProg.FALatCoP{i}(1,:)]);
            DynamicPPTrials(Selection(q)).Zones.R(i).MedFore = round([DynamicPPTrials(Selection(q)).RProg.FAMedCoP{i}(1,:); DynamicPPTrials(Selection(q)).R.MedCoP{i}(2,:); DynamicPPTrials(Selection(q)).Rht.CoP{i}(2,:); DynamicPPTrials(Selection(q)).RProg.FACoP{i}(1,:)]);
             % Whole Foot
            DynamicPPTrials(Selection(q)).Zones.R(i).Whole = round([DynamicPPTrials(Selection(q)).R.LatCoP{i}(1,:); DynamicPPTrials(Selection(q)).R.LatCoP{i}(2,:); DynamicPPTrials(Selection(q)).R.MedCoP{i}(2,:); DynamicPPTrials(Selection(q)).R.MedCoP{i}(1,:);]);
            % Prog Angle
            DynamicPPTrials(Selection(q)).RProg.Ang{i} = DynamicPPTrials(Selection(q)).RProg.AngCoP{i};
        elseif strcmp(PPSettings.MaskChoice, '66% CoP (cyan)')  || strcmp(PPSettings.MaskChoice, '66% CoP')
            % Lateral and Medial Heel
            DynamicPPTrials(Selection(q)).Zones.R(i).LatHeel = round([DynamicPPTrials(Selection(q)).R.Lat66{i}(1,:); DynamicPPTrials(Selection(q)).Rht.CoP66{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HA66{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HALat66{i}(1,:)]);
            DynamicPPTrials(Selection(q)).Zones.R(i).MedHeel = round([DynamicPPTrials(Selection(q)).R.Med66{i}(1,:); DynamicPPTrials(Selection(q)).Rht.CoP66{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HA66{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HAMed66{i}(1,:)]);
            % Lateral and Medial Arch
            DynamicPPTrials(Selection(q)).Zones.R(i).LatArch = round([DynamicPPTrials(Selection(q)).RProg.HALat66{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HA66{i}(1,:); DynamicPPTrials(Selection(q)).RProg.FA66{i}(1,:); DynamicPPTrials(Selection(q)).RProg.FALat66{i}(1,:)]);
            DynamicPPTrials(Selection(q)).Zones.R(i).MedArch = round([DynamicPPTrials(Selection(q)).RProg.HAMed66{i}(1,:); DynamicPPTrials(Selection(q)).RProg.FAMed66{i}(1,:); DynamicPPTrials(Selection(q)).RProg.FA66{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HA66{i}(1,:)]);
            % Lateral and Medial Forefoot
            DynamicPPTrials(Selection(q)).Zones.R(i).LatFore = round([DynamicPPTrials(Selection(q)).RProg.FA66{i}(1,:); DynamicPPTrials(Selection(q)).Rht.CoP66{i}(2,:); DynamicPPTrials(Selection(q)).R.Lat66{i}(2,:); DynamicPPTrials(Selection(q)).RProg.FALat66{i}(1,:)]);
            DynamicPPTrials(Selection(q)).Zones.R(i).MedFore = round([DynamicPPTrials(Selection(q)).RProg.FAMed66{i}(1,:); DynamicPPTrials(Selection(q)).R.Med66{i}(2,:); DynamicPPTrials(Selection(q)).Rht.CoP66{i}(2,:); DynamicPPTrials(Selection(q)).RProg.FA66{i}(1,:)]);
             % Whole Foot
            DynamicPPTrials(Selection(q)).Zones.R(i).Whole = round([DynamicPPTrials(Selection(q)).R.LatCoP66{i}(1,:); DynamicPPTrials(Selection(q)).R.LatCoP66{i}(2,:); DynamicPPTrials(Selection(q)).R.MedCoP66{i}(2,:); DynamicPPTrials(Selection(q)).R.MedCoP66{i}(1,:);]);
            % Prog Angle
            DynamicPPTrials(Selection(q)).RProg.Ang{i} = DynamicPPTrials(Selection(q)).RProg.Ang66{i};
        elseif strcmp(PPSettings.MaskChoice, 'Inter-Peak (green)') || strcmp(PPSettings.MaskChoice, 'Inter-Peak')
            if Block.IP == 0
                % Lateral and Medial Heel
                DynamicPPTrials(Selection(q)).Zones.R(i).LatHeel = round([DynamicPPTrials(Selection(q)).R.LatIP{i}(1,:); DynamicPPTrials(Selection(q)).Rht.IP{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HAIP{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HALatIP{i}(1,:)]);
                DynamicPPTrials(Selection(q)).Zones.R(i).MedHeel = round([DynamicPPTrials(Selection(q)).R.MedIP{i}(1,:); DynamicPPTrials(Selection(q)).Rht.IP{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HAIP{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HAMedIP{i}(1,:)]);
                % Lateral and Medial Arch
                DynamicPPTrials(Selection(q)).Zones.R(i).LatArch = round([DynamicPPTrials(Selection(q)).RProg.HALatIP{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HAIP{i}(1,:); DynamicPPTrials(Selection(q)).RProg.FAIP{i}(1,:); DynamicPPTrials(Selection(q)).RProg.FALatIP{i}(1,:)]);
                DynamicPPTrials(Selection(q)).Zones.R(i).MedArch = round([DynamicPPTrials(Selection(q)).RProg.HAMedIP{i}(1,:); DynamicPPTrials(Selection(q)).RProg.FAMedIP{i}(1,:); DynamicPPTrials(Selection(q)).RProg.FAIP{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HAIP{i}(1,:)]);
                % Lateral and Medial Forefoot
                DynamicPPTrials(Selection(q)).Zones.R(i).LatFore = round([DynamicPPTrials(Selection(q)).RProg.FAIP{i}(1,:); DynamicPPTrials(Selection(q)).Rht.IP{i}(2,:); DynamicPPTrials(Selection(q)).R.LatIP{i}(2,:); DynamicPPTrials(Selection(q)).RProg.FALatIP{i}(1,:)]);
                DynamicPPTrials(Selection(q)).Zones.R(i).MedFore = round([DynamicPPTrials(Selection(q)).RProg.FAMedIP{i}(1,:); DynamicPPTrials(Selection(q)).R.MedIP{i}(2,:); DynamicPPTrials(Selection(q)).Rht.IP{i}(2,:); DynamicPPTrials(Selection(q)).RProg.FAIP{i}(1,:)]);
                 % Whole Foot
            DynamicPPTrials(Selection(q)).Zones.R(i).Whole = round([DynamicPPTrials(Selection(q)).R.LatIP{i}(1,:); DynamicPPTrials(Selection(q)).R.LatIP{i}(2,:); DynamicPPTrials(Selection(q)).R.MedIP{i}(2,:); DynamicPPTrials(Selection(q)).R.MedIP{i}(1,:);]);
            % Prog Angle
                DynamicPPTrials(Selection(q)).RProg.Ang{i} = DynamicPPTrials(Selection(q)).RProg.AngIP{i};
            end
        elseif strcmp(PPSettings.MaskChoice, 'Heel-Centroid (red)') || strcmp(PPSettings.MaskChoice, 'Heel-Centroid')
            if Block.HC == 0
                % Lateral and Medial Heel
                DynamicPPTrials(Selection(q)).Zones.R(i).LatHeel = round([DynamicPPTrials(Selection(q)).R.LatHeelCent{i}(1,:); DynamicPPTrials(Selection(q)).Rht.HeelCent{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HAHeelCent{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HALatHeelCent{i}(1,:)]);
                DynamicPPTrials(Selection(q)).Zones.R(i).MedHeel = round([DynamicPPTrials(Selection(q)).R.MedHeelCent{i}(1,:); DynamicPPTrials(Selection(q)).Rht.HeelCent{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HAHeelCent{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HAMedHeelCent{i}(1,:)]);
                % Lateral and Medial Arch
                DynamicPPTrials(Selection(q)).Zones.R(i).LatArch = round([DynamicPPTrials(Selection(q)).RProg.HALatHeelCent{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HAHeelCent{i}(1,:); DynamicPPTrials(Selection(q)).RProg.FAHeelCent{i}(1,:); DynamicPPTrials(Selection(q)).RProg.FALatHeelCent{i}(1,:)]);
                DynamicPPTrials(Selection(q)).Zones.R(i).MedArch = round([DynamicPPTrials(Selection(q)).RProg.HAMedHeelCent{i}(1,:); DynamicPPTrials(Selection(q)).RProg.FAMedHeelCent{i}(1,:); DynamicPPTrials(Selection(q)).RProg.FAHeelCent{i}(1,:); DynamicPPTrials(Selection(q)).RProg.HAHeelCent{i}(1,:)]);
                % Lateral and Medial Forefoot
                DynamicPPTrials(Selection(q)).Zones.R(i).LatFore = round([DynamicPPTrials(Selection(q)).RProg.FAHeelCent{i}(1,:); DynamicPPTrials(Selection(q)).Rht.HeelCent{i}(2,:); DynamicPPTrials(Selection(q)).R.LatHeelCent{i}(2,:); DynamicPPTrials(Selection(q)).RProg.FALatHeelCent{i}(1,:)]);
                DynamicPPTrials(Selection(q)).Zones.R(i).MedFore = round([DynamicPPTrials(Selection(q)).RProg.FAMedHeelCent{i}(1,:); DynamicPPTrials(Selection(q)).R.MedHeelCent{i}(2,:); DynamicPPTrials(Selection(q)).Rht.HeelCent{i}(2,:); DynamicPPTrials(Selection(q)).RProg.FAHeelCent{i}(1,:)]);
                % Whole Foot
                DynamicPPTrials(Selection(q)).Zones.R(i).Whole = round([DynamicPPTrials(Selection(q)).R.LatHeelCent{i}(1,:); DynamicPPTrials(Selection(q)).R.LatHeelCent{i}(2,:); DynamicPPTrials(Selection(q)).R.MedHeelCent{i}(2,:); DynamicPPTrials(Selection(q)).R.MedHeelCent{i}(1,:);]);
                % Prog Angle
                DynamicPPTrials(Selection(q)).RProg.Ang{i} = DynamicPPTrials(Selection(q)).RProg.AngHeelCent{i};
            end
        end
    end
end
% % *uncomment to check to make sure the polygon points are correct
% % plot(R.LatHeel(:,1),R.LatHeel(:,2), 'om');
% % plot(R.MedFore(:,1),R.MedFore(:,2), '*m'); %
clearvars i q A

%% create plot of 4 trials up close for display in report
clc; disp('Creating final image of all trials');
filename = strrep(filename,'.mat','');
% saveas(strcat(folder,'\',filename,'_AllTrials.png'));
if strcmp(PPSettings.ExportReport, 'Yes') == 1
    close;
    if NumDynTrials > 4
        %         Chosen4 = listdlg('PromptString','Select 4 trials to display','ListString',ListStr);\
        Chosen4 = Selection(1:4);
    else
        Chosen4 = Selection;
    end
    if strcmp(PPSettings.C3DInput,'No')
        C3Ddata = 0;
    end
    PlotPressures(DynamicPPTrials, C3Ddata, PPSettings, Chosen4);
end

%% Output figure of all trials to Report
if strcmp(PPSettings.ExportReport, 'Yes') == 1
    
%     [fid msg] = fopen('path\to\excel\file\filename.xls','a');
% if fid==-1
%     disp('The file is already open.')
% else
%     fclose(fid);
%     disp('If the file did not already exist, it has been created.')
% end
%     system('taskkill /F /IM EXCEL.EXE');
    
    PPPlots.TempSpat = 'Yes';
    PPPlots.CloseUp = 'Yes';
    PPPlots.AllForces = 'Yes';
    % copy and paste template for overwriting
    InputFileName = fullfile(pwd, 'Impressions Output Template.xlsx');
    OutputFileName = fullfile(pwd, 'ImpressionsOutput.xlsx');
    copyfile(InputFileName, OutputFileName);
    clearvars InputFileName OutputFileName
    %     subplotsqueeze(Cont,1.1);
    % export figure to excel
    xlsPasteTo('ImpressionsOutput.xlsx','Sheet1',680,600, 'A17');
end
clearvars i j locs q Pt1 Pt2 ToeDrag Xadj Yadj Xoffset Yoffset Pix_SS PtOnOff question z Reason Static xdist ydist yFPval yThresh xThresh LoadNew ans F S s CC CCC Cent LFootAdj RFootAdj


%% calculate step width and step length
clc; disp('Temporal Spatial calculations');
% we assume that the participant walks nearly straight on the mat
% determine number of steps
for i = 1:length(DynamicPPTrials)
    if isempty(DynamicPPTrials(i).Lht) || isempty(DynamicPPTrials(i).Rht)
        Empty(i) = 1;
    else
        Empty(i) = 0;
    end
end
if sum(~Empty) == 0
    TempSpat = 0;
else
    TempSpat = 1;
    TStrials = find(~Empty);
    Order = [1 1 2 2 3 3 4 4 5 5 6 6]; % order for indexing
    for j = 1:length(TStrials)
        NumSteps = DynamicPPTrials(TStrials(j)).NumLeft + DynamicPPTrials(TStrials(j)).NumRight - 1;
        if strcmp([Regions{1,TStrials(j)}(1).Side], 'Left') == 1 % if LEFT foot is first
            for i = 1:NumSteps+1
                if mod(i,2) == 1 % determines even or odd
                    Stride = DynamicPPTrials(TStrials(j)).Lht.General{Order(i)};
                else
                    Stride = DynamicPPTrials(TStrials(j)).Rht.General{Order(i)};
                end
                StepHeels(i,1:2) = Stride(1,:); % places heel corrdinates in array
            end
        else % if RIGHT foot is first
            for i = 1:NumSteps+1
                if mod(i,2) == 1 % determines even or odd
                    Stride = DynamicPPTrials(TStrials(j)).Rht.General{Order(i)};
                else
                    Stride = DynamicPPTrials(TStrials(j)).Lht.General{Order(i)};
                end
                StepHeels(i,1:2) = Stride(1,:); % places heel corrdinates in array
            end
        end
        
        for i = 1:NumSteps
            StepWidths(i) = abs(StepHeels(i+1,1) - StepHeels(i,1)) .* Adj.Width;
            StrideLengths(i) = abs(StepHeels(i+1,2) - StepHeels(i,2)) .* Adj.Length;
        end
        clearvars Stride StepHeels
        % Calculate Step Widths and Stride Lengths
        if exist('StepWidths','var') == 1
            DynamicPPTrials(TStrials(j)).StepWidth.Avg = nanmean(StepWidths(:));
            DynamicPPTrials(TStrials(j)).StepWidth.AvgRel = DynamicPPTrials(TStrials(j)).StepWidth.Avg / Subject.Height * 100;
            DynamicPPTrials(TStrials(j)).StepWidth.Std = nanstd(StepWidths(:));
            DynamicPPTrials(TStrials(j)).StepWidth.StdRel = DynamicPPTrials(TStrials(j)).StepWidth.Std / Subject.Height *100;
            DynamicPPTrials(TStrials(j)).StepWidth.Measures = StepWidths(:);
            
            DynamicPPTrials(TStrials(j)).StrideLength.Avg = nanmean(StrideLengths(:));
            DynamicPPTrials(TStrials(j)).StrideLength.AvgRel = DynamicPPTrials(TStrials(j)).StrideLength.Avg / Subject.Height * 100;
            DynamicPPTrials(TStrials(j)).StrideLength.Std = nanstd(StrideLengths(:));
            DynamicPPTrials(TStrials(j)).StrideLength.StdRel = DynamicPPTrials(TStrials(j)).StrideLength.Std / Subject.Height *100;
            DynamicPPTrials(TStrials(j)).StrideLength.Measures = StrideLengths;
        else
            DynamicPPTrials(TStrials(j)).StepWidth.Avg = NaN;
            DynamicPPTrials(TStrials(j)).StepWidth.AvgRel = NaN;
            DynamicPPTrials(TStrials(j)).StepWidth.Std = NaN;
            DynamicPPTrials(TStrials(j)).StepWidth.StdRel = NaN;
            DynamicPPTrials(TStrials(j)).StepWidth.Measures = NaN;
            
            DynamicPPTrials(TStrials(j)).StrideLength.Avg = NaN;
            DynamicPPTrials(TStrials(j)).StrideLength.AvgRel = NaN;
            DynamicPPTrials(TStrials(j)).StrideLength.Std = NaN;
            DynamicPPTrials(TStrials(j)).StrideLength.StdRel = NaN;
            DynamicPPTrials(TStrials(j)).StrideLength.Measures = NaN;
        end
        clearvars StepWidths StrideLengths
    end
    
    clearvars Order j i z q
    
    %% Load Control Temporal Spatial Data
    TSN = GetTempSpatNorms(Subject.Age);
    
    %% Temporal spatial calculations
    % Calculate gait speed and cadence
    for j = 1:length(TStrials)
        MpS(j) = ((Regions{1,TStrials(j)}(end).BoundingBox(2) - Regions{1,TStrials(j)}(1).BoundingBox(2))*Adj.Length) /  (100*(Regions{1,TStrials(j)}(end).TimeStrike - Regions{1,TStrials(j)}(1).TimeStrike));
        Cadence(j) = 60*(length(Regions{1,TStrials(j)}) -1) / (Regions{1,TStrials(j)}(end).TimeStrike - Regions{1,TStrials(j)}(1).TimeStrike);
        StepW(j) = DynamicPPTrials(TStrials(j)).StepWidth.Avg;
        StrideL(j) = DynamicPPTrials(TStrials(j)).StrideLength.Avg;
    end
    
    TS.GaitSpeed.MpS.Avg = nanmean(MpS); % averages
    TS.GaitSpeed.MpM.Avg = TS.GaitSpeed.MpS.Avg * 60;
    TS.Cadence.Avg = nanmean(Cadence);
    TS.StepWidth.Avg = nanmean(StepW);
    TS.StrideLength.Avg = nanmean(StrideL);
    TS.GaitSpeed.MpS.Measures = MpS; % measures
    TS.GaitSpeed.MpM.Measures = TS.GaitSpeed.MpS.Measures * 60;
    TS.Cadence.Measures = Cadence;
    TS.StepWidth.Measures = StepW;
    TS.StrideLength.Measures = StrideL;
    TS.GaitSpeed.MpS.Std = nanstd(MpS); % standard deviations
    TS.GaitSpeed.MpM.Std = nanstd(TS.GaitSpeed.MpM.Measures);
    TS.Cadence.Std = nanstd(Cadence);
    TS.StepWidth.Std = nanstd(StepW);
    TS.StrideLength.Std = nanstd(StrideL);
    
    % Initial and final double support and single limb support
    for j = 1:length(TStrials)
        for k = 1:length(Regions{1,TStrials(j)})
            if isempty(Regions{1,TStrials(j)}(k).Off)
                Regions{1,TStrials(j)}(k).Off = NaN;
            end
            if isempty(Regions{1,TStrials(j)}(k).TimeOff)
                Regions{1,TStrials(j)}(k).TimeOff = NaN;
            end
        end
        if length(Regions{1,TStrials(j)}) < 3
            % skip because initial and double supports metrics cannot be calculated without at least 3 steps
            % define as NaN so that it is not defined as 0, and averaged
            Stance(j) = NaN;
            Swing(j) = NaN;
            IDS(j) = NaN;
            SS(j) = NaN;
            SDS(j) = NaN;
        elseif length(Regions{1,TStrials(j)}) == 3
            Stance(j) = (Regions{1,TStrials(j)}(1).Off - Regions{1,TStrials(j)}(1).Strike) / (Regions{1,TStrials(j)}(3).Strike - Regions{1,TStrials(j)}(1).Strike);
            Swing(j) = 1 - Stance(j);
            IDS(j) = abs(Regions{1,TStrials(j)}(1).Off - Regions{1,TStrials(j)}(2).Strike) / (Regions{1,TStrials(j)}(3).Strike - Regions{1,TStrials(j)}(1).Strike);
            SS(j) = abs(Regions{1,TStrials(j)}(3).Strike - Regions{1,TStrials(j)}(1).Off) / (Regions{1,TStrials(j)}(3).Strike - Regions{1,TStrials(j)}(1).Strike);
            SDS(j) = abs(Regions{1,TStrials(j)}(2).Off - Regions{1,TStrials(j)}(3).Strike) / (Regions{1,TStrials(j)}(3).Strike - Regions{1,TStrials(j)}(1).Strike);
        elseif length(Regions{1,TStrials(j)}) == 4
            %         for i = 1:length(Regions{1,TStrials(j)})
            Stance1 = (Regions{1,TStrials(j)}(1).Off - Regions{1,TStrials(j)}(1).Strike) / (Regions{1,TStrials(j)}(3).Strike - Regions{1,TStrials(j)}(1).Strike);
            Stance2 = (Regions{1,TStrials(j)}(2).Off - Regions{1,TStrials(j)}(2).Strike) / (Regions{1,TStrials(j)}(4).Strike - Regions{1,TStrials(j)}(2).Strike);
            Stance(j) = (Stance1 + Stance2)/2;
            Swing(j) = 1 - Stance(j);
            IDS1 = (Regions{1,TStrials(j)}(1).Off - Regions{1,TStrials(j)}(2).Strike) / (Regions{1,TStrials(j)}(3).Strike - Regions{1,TStrials(j)}(1).Strike);
            IDS2 = (Regions{1,TStrials(j)}(2).Off - Regions{1,TStrials(j)}(3).Strike) / (Regions{1,TStrials(j)}(4).Strike - Regions{1,TStrials(j)}(2).Strike);
            IDS(j) = (IDS1 + IDS2) / 2;
            SS1 = (Regions{1,TStrials(j)}(3).Strike - Regions{1,TStrials(j)}(1).Off) / (Regions{1,TStrials(j)}(3).Strike - Regions{1,TStrials(j)}(1).Strike);
            SS2 = (Regions{1,TStrials(j)}(4).Strike - Regions{1,TStrials(j)}(2).Off) / (Regions{1,TStrials(j)}(4).Strike - Regions{1,TStrials(j)}(2).Strike);
            SS(j) = (SS1 + SS2) /2;
            SDS1 = (Regions{1,TStrials(j)}(2).Off - Regions{1,TStrials(j)}(3).Strike) / (Regions{1,TStrials(j)}(3).Strike - Regions{1,TStrials(j)}(1).Strike);
            SDS2 = (Regions{1,TStrials(j)}(3).Off - Regions{1,TStrials(j)}(4).Strike) / (Regions{1,TStrials(j)}(4).Strike - Regions{1,TStrials(j)}(2).Strike);
            SDS(j) = (SDS1 + SDS2) / 2;
            clearvars Stance1 Stance2 IDS1 IDS2 SS1 SS2 SDS1 SDS2
        elseif length(Regions{1,TStrials(j)}) == 5
            Stance1 = (Regions{1,TStrials(j)}(1).Off - Regions{1,TStrials(j)}(1).Strike) / (Regions{1,TStrials(j)}(3).Strike - Regions{1,TStrials(j)}(1).Strike);
            Stance2 = (Regions{1,TStrials(j)}(2).Off - Regions{1,TStrials(j)}(2).Strike) / (Regions{1,TStrials(j)}(4).Strike - Regions{1,TStrials(j)}(2).Strike);
            Stance3 = (Regions{1,TStrials(j)}(3).Off - Regions{1,TStrials(j)}(3).Strike) / (Regions{1,TStrials(j)}(5).Strike - Regions{1,TStrials(j)}(3).Strike);
            Stance(j) = (Stance1 + Stance2 + Stance3) / 3;
            Swing(j) = 1 - Stance(j);
            IDS1 = (Regions{1,TStrials(j)}(1).Off - Regions{1,TStrials(j)}(2).Strike) / (Regions{1,TStrials(j)}(3).Strike - Regions{1,TStrials(j)}(1).Strike);
            IDS2 = (Regions{1,TStrials(j)}(2).Off - Regions{1,TStrials(j)}(3).Strike) / (Regions{1,TStrials(j)}(4).Strike - Regions{1,TStrials(j)}(2).Strike);
            IDS3 = (Regions{1,TStrials(j)}(3).Off - Regions{1,TStrials(j)}(4).Strike) / (Regions{1,TStrials(j)}(5).Strike - Regions{1,TStrials(j)}(3).Strike);
            IDS(j) = (IDS1 + IDS2 + IDS3) / 3;
            SS1 = (Regions{1,TStrials(j)}(3).Strike - Regions{1,TStrials(j)}(1).Off) / (Regions{1,TStrials(j)}(3).Strike - Regions{1,TStrials(j)}(1).Strike);
            SS2 = (Regions{1,TStrials(j)}(4).Strike - Regions{1,TStrials(j)}(2).Off) / (Regions{1,TStrials(j)}(4).Strike - Regions{1,TStrials(j)}(2).Strike);
            SS3 = (Regions{1,TStrials(j)}(5).Strike - Regions{1,TStrials(j)}(3).Off) / (Regions{1,TStrials(j)}(5).Strike - Regions{1,TStrials(j)}(3).Strike);
            SS(j) = (SS1 + SS2 + SS3) / 3;
            SDS1 = (Regions{1,TStrials(j)}(2).Off - Regions{1,TStrials(j)}(3).Strike) / (Regions{1,TStrials(j)}(3).Strike - Regions{1,TStrials(j)}(1).Strike);
            SDS2 = (Regions{1,TStrials(j)}(3).Off - Regions{1,TStrials(j)}(4).Strike) / (Regions{1,TStrials(j)}(4).Strike - Regions{1,TStrials(j)}(2).Strike);
            SDS3 = (Regions{1,TStrials(j)}(4).Off - Regions{1,TStrials(j)}(5).Strike) / (Regions{1,TStrials(j)}(5).Strike - Regions{1,TStrials(j)}(3).Strike);
            SDS(j) = (SDS1 + SDS2 + SDS3) / 3;
            clearvars Stance1 Stance2 IDS1 IDS2 SS1 SS2 SDS1 SDS2 IDS3 SS3 Stance3 SDS3
        end
    end
    
    % calculate the average support times for plotting
    TS.Stance.Avg = nanmean(Stance); % averages
    TS.Swing.Avg = nanmean(Swing);
    TS.IDS.Avg = nanmean(IDS);
    TS.SS.Avg = nanmean(SS);
    TS.SDS.Avg = nanmean(SDS);
    TS.Stance.Measures = Stance; % measures
    TS.Swing.Measures = Swing;
    TS.IDS.Measures = IDS;
    TS.SS.Measures = SS;
    TS.SDS.Measures = SDS;
    TS.Stance.Std = nanstd(Stance); % standard deviations
    TS.Swing.Std = nanstd(Swing);
    TS.IDS.Std = nanstd(IDS);
    TS.SS.Std = nanstd(SS);
    TS.SDS.Std = nanstd(SDS);
    
    clearvars Stance Swing IDS SS SDS StepW StrideL Cadence
end

%% Compile Accuracy Metrics
clearvars ExportAccuracy ExportProg ExportPtDiff
if strcmp(PPSettings.C3DInput, 'Yes')
    for z = 1:length(Selection)
        ExportAccuracy(z).FileName = C3Ddata(z).FullFileName;
        ExportAccuracy(z).Direction = C3Ddata(z).Direction;
        ExportProg(z).FileName = C3Ddata(z).FullFileName;
        ExportProg(z).Direction = C3Ddata(z).Direction;
        ExportPtDiff(z).FileName = C3Ddata(z).FullFileName;
        ExportPtDiff(z).Direction = C3Ddata(z).Direction;
        
        if isempty(C3Ddata(z).LHA) == 0
            i = 1;
            ExportAccuracy(z).LHAxb = C3Ddata(z).LHAb(i,2);
            ExportAccuracy(z).LHAyb = C3Ddata(z).LHAb(i,1);
            ExportAccuracy(z).LFAxb = C3Ddata(z).LFAb(i,2);
            ExportAccuracy(z).LFAyb = C3Ddata(z).LFAb(i,1);
%             
%             for i = 1:size(C3Ddata(z).LHA, 1)
                x = C3Ddata(z).LFA(i,2) - C3Ddata(z).LHA(i,2);
                y = C3Ddata(z).LFA(i,1) - C3Ddata(z).LHA(i,1);
                LProgC3D(i) = (atan2(y, x) * 180 / pi) - 90;
%             end
            ExportProg(z).LProgC3D = nanmean(LProgC3D);
        end
        
        if isempty(C3Ddata(z).RHA) == 0
            i = 1; 
            ExportAccuracy(z).RHAxb = C3Ddata(z).RHAb(i,2);
            ExportAccuracy(z).RHAyb = C3Ddata(z).RHAb(i,1);
            ExportAccuracy(z).RFAxb = C3Ddata(z).RFAb(i,2);
            ExportAccuracy(z).RFAyb = C3Ddata(z).RFAb(i,1);
            
%             for i = 1:size(C3Ddata(z).RHA, 1)
                x = C3Ddata(z).RFA(i,2) - C3Ddata(z).RHA(i,2);
                y = C3Ddata(z).RFA(i,1) - C3Ddata(z).RHA(i,1);
                RProgC3D(i) = -((atan2(y, x) * 180 / pi) - 90);
%             end
            ExportProg(z).RProgC3D = nanmean(RProgC3D);
        end
        
        if isempty(DynamicPPTrials(z).LProg) == 0
            i = 1; 
            ExportProg(z).LProgAngManual = cell2mat(DynamicPPTrials(z).LProg.AngActualManual(i));
            ExportProg(z).LProgAngGeneral = cell2mat(DynamicPPTrials(z).LProg.AngActualGeneral(i));
            ExportProg(z).LProgAngCoP = cell2mat(DynamicPPTrials(z).LProg.AngActualCoP(i));
            ExportProg(z).LProgAng66 = cell2mat(DynamicPPTrials(z).LProg.AngActual66(i));
            if Block.IP == 0
                ExportProg(z).LProgAngIP = cell2mat(DynamicPPTrials(z).LProg.AngActualIP(i));
            end
            if Block.HC == 0
                ExportProg(z).LProgAngHeelCent = cell2mat(DynamicPPTrials(z).LProg.AngActualHeelCent(i));
            end
            
            ExportAccuracy(z).LHAManualx = DynamicPPTrials(z).LProg.HAManual{i}(1); % Manual
            ExportAccuracy(z).LHAManualy = DynamicPPTrials(z).LProg.HAManual{i}(2);
            ExportAccuracy(z).LFAManualx = DynamicPPTrials(z).LProg.FAManual{i}(1);
            ExportAccuracy(z).LFAManualy = DynamicPPTrials(z).LProg.FAManual{i}(2);
            ExportAccuracy(z).LHAGeneralx = DynamicPPTrials(z).LProg.HAGeneral{i}(1); % general image processing
            ExportAccuracy(z).LHAGeneraly = DynamicPPTrials(z).LProg.HAGeneral{i}(2);
            ExportAccuracy(z).LFAGeneralx = DynamicPPTrials(z).LProg.FAGeneral{i}(1);
            ExportAccuracy(z).LFAGeneraly = DynamicPPTrials(z).LProg.FAGeneral{i}(2);
            ExportAccuracy(z).LHACoPx = DynamicPPTrials(z).LProg.HACoP{i}(1); % CoP
            ExportAccuracy(z).LHACoPy = DynamicPPTrials(z).LProg.HACoP{i}(2);
            ExportAccuracy(z).LFACoPx = DynamicPPTrials(z).LProg.FACoP{i}(1);
            ExportAccuracy(z).LFACoPy = DynamicPPTrials(z).LProg.FACoP{i}(2);
            ExportAccuracy(z).LHA66x = DynamicPPTrials(z).LProg.HA66{i}(1); % 66% of CoP
            ExportAccuracy(z).LHA66y = DynamicPPTrials(z).LProg.HA66{i}(2);
            ExportAccuracy(z).LFA66x = DynamicPPTrials(z).LProg.FA66{i}(1);
            ExportAccuracy(z).LFA66y = DynamicPPTrials(z).LProg.FA66{i}(2);
            if Block.IP == 0
                ExportAccuracy(z).LHAIPx = DynamicPPTrials(z).LProg.HAIP{i}(1); % Inter-peak
                ExportAccuracy(z).LHAIPy = DynamicPPTrials(z).LProg.HAIP{i}(2);
                ExportAccuracy(z).LFAIPx = DynamicPPTrials(z).LProg.FAIP{i}(1);
                ExportAccuracy(z).LFAIPy = DynamicPPTrials(z).LProg.FAIP{i}(2);
            end
            if Block.HC == 0
                ExportAccuracy(z).LHAHeelCentx = DynamicPPTrials(z).LProg.HAHeelCent{i}(1); % heel to centroid
                ExportAccuracy(z).LHAHeelCenty = DynamicPPTrials(z).LProg.HAHeelCent{i}(2);
                ExportAccuracy(z).LFAHeelCentx = DynamicPPTrials(z).LProg.FAHeelCent{i}(1);
                ExportAccuracy(z).LFAHeelCenty = DynamicPPTrials(z).LProg.FAHeelCent{i}(2);
            end
        end
        if isempty(DynamicPPTrials(z).RProg) == 0
            i = 1; 
            ExportProg(z).RProgAngManual = -cell2mat(DynamicPPTrials(z).RProg.AngActualManual(i));
            ExportProg(z).RProgAngGeneral = -cell2mat(DynamicPPTrials(z).RProg.AngActualGeneral(i));
            ExportProg(z).RProgAngCoP = -cell2mat(DynamicPPTrials(z).RProg.AngActualCoP(i));
            ExportProg(z).RProgAng66 = -cell2mat(DynamicPPTrials(z).RProg.AngActual66(i));
            if Block.IP == 0
                ExportProg(z).RProgAngIP = -cell2mat(DynamicPPTrials(z).RProg.AngActualIP(i));
            end
            if Block.HC == 0
                ExportProg(z).RProgAngHeelCent = -cell2mat(DynamicPPTrials(z).RProg.AngActualHeelCent(i));
            end
            
            ExportAccuracy(z).RHAManualx = DynamicPPTrials(z).RProg.HAManual{i}(1); % manual
            ExportAccuracy(z).RHAManualy = DynamicPPTrials(z).RProg.HAManual{i}(2);
            ExportAccuracy(z).RFAManualx = DynamicPPTrials(z).RProg.FAManual{i}(1);
            ExportAccuracy(z).RFAManualy = DynamicPPTrials(z).RProg.FAManual{i}(2);
            ExportAccuracy(z).RHAGeneralx = DynamicPPTrials(z).RProg.HAGeneral{i}(1); % general image processing
            ExportAccuracy(z).RHAGeneraly = DynamicPPTrials(z).RProg.HAGeneral{i}(2);
            ExportAccuracy(z).RFAGeneralx = DynamicPPTrials(z).RProg.FAGeneral{i}(1);
            ExportAccuracy(z).RFAGeneraly = DynamicPPTrials(z).RProg.FAGeneral{i}(2);
            ExportAccuracy(z).RHACoPx = DynamicPPTrials(z).RProg.HACoP{i}(1); % CoP
            ExportAccuracy(z).RHACoPy = DynamicPPTrials(z).RProg.HACoP{i}(2);
            ExportAccuracy(z).RFACoPx = DynamicPPTrials(z).RProg.FACoP{i}(1);
            ExportAccuracy(z).RFACoPy = DynamicPPTrials(z).RProg.FACoP{i}(2);
            ExportAccuracy(z).RHA66x = DynamicPPTrials(z).RProg.HA66{i}(1); % 66
            ExportAccuracy(z).RHA66y = DynamicPPTrials(z).RProg.HA66{i}(2);
            ExportAccuracy(z).RFA66x = DynamicPPTrials(z).RProg.FA66{i}(1);
            ExportAccuracy(z).RFA66y = DynamicPPTrials(z).RProg.FA66{i}(2);
            if Block.IP == 0
                ExportAccuracy(z).RHAIPx = DynamicPPTrials(z).RProg.HAIP{i}(1); % IP
                ExportAccuracy(z).RHAIPy = DynamicPPTrials(z).RProg.HAIP{i}(2);
                ExportAccuracy(z).RFAIPx = DynamicPPTrials(z).RProg.FAIP{i}(1);
                ExportAccuracy(z).RFAIPy = DynamicPPTrials(z).RProg.FAIP{i}(2);
            end
            if Block.HC == 0
                ExportAccuracy(z).RHAHeelCentx = DynamicPPTrials(z).RProg.HAHeelCent{i}(1); % HeelCent
                ExportAccuracy(z).RHAHeelCenty = DynamicPPTrials(z).RProg.HAHeelCent{i}(2);
                ExportAccuracy(z).RFAHeelCentx = DynamicPPTrials(z).RProg.FAHeelCent{i}(1);
                ExportAccuracy(z).RFAHeelCenty = DynamicPPTrials(z).RProg.FAHeelCent{i}(2);
            end
        end
        
        % save prog and accuracy differences
        if isempty(DynamicPPTrials(z).LProg) == 0
            ExportProg(z).LProgDiffManual = abs(ExportProg(z).LProgC3D - ExportProg(z).LProgAngManual);
            ExportProg(z).LProgDiffGeneral = abs(ExportProg(z).LProgC3D - ExportProg(z).LProgAngGeneral);
            ExportProg(z).LProgDiffCoP = abs(ExportProg(z).LProgC3D - ExportProg(z).LProgAngCoP);
            ExportProg(z).LProgDiff66 = abs(ExportProg(z).LProgC3D - ExportProg(z).LProgAng66);
            if Block.IP == 0
                ExportProg(z).LProgDiffIP =abs( ExportProg(z).LProgC3D - ExportProg(z).LProgAngIP);
            end
            if Block.HC == 0
                ExportProg(z).LProgDiffHeelCent = abs(ExportProg(z).LProgC3D - ExportProg(z).LProgAngHeelCent);
            end
            
            ExportPtDiff(z).LHAdiffManual = sqrt((ExportAccuracy(z).LHAxb - ExportAccuracy(z).LHAManualx)^2 + (ExportAccuracy(z).LHAyb - ExportAccuracy(z).LHAManualy)^2); % manual
            ExportPtDiff(z).LFAdiffManual = sqrt((ExportAccuracy(z).LFAxb - ExportAccuracy(z).LFAManualx)^2 + (ExportAccuracy(z).LFAyb - ExportAccuracy(z).LFAManualy)^2);
            ExportPtDiff(z).LHAdiffGeneral = sqrt((ExportAccuracy(z).LHAxb - ExportAccuracy(z).LHAGeneralx)^2 + (ExportAccuracy(z).LHAyb - ExportAccuracy(z).LHAGeneraly)^2); % general
            ExportPtDiff(z).LFAdiffGeneral = sqrt((ExportAccuracy(z).LFAxb - ExportAccuracy(z).LFAGeneralx)^2 + (ExportAccuracy(z).LFAyb - ExportAccuracy(z).LFAGeneraly)^2);
            ExportPtDiff(z).LHAdiffCoP = sqrt((ExportAccuracy(z).LHAxb - ExportAccuracy(z).LHACoPx)^2 + (ExportAccuracy(z).LHAyb - ExportAccuracy(z).LHACoPy)^2); % CoP
            ExportPtDiff(z).LFAdiffCoP = sqrt((ExportAccuracy(z).LFAxb - ExportAccuracy(z).LFACoPx)^2 + (ExportAccuracy(z).LFAyb - ExportAccuracy(z).LFACoPy)^2);
            ExportPtDiff(z).LHAdiff66 = sqrt((ExportAccuracy(z).LHAxb - ExportAccuracy(z).LHA66x)^2 + (ExportAccuracy(z).LHAyb - ExportAccuracy(z).LHA66y)^2); % 66
            ExportPtDiff(z).LFAdiff66 = sqrt((ExportAccuracy(z).LFAxb - ExportAccuracy(z).LFA66x)^2 + (ExportAccuracy(z).LFAyb - ExportAccuracy(z).LFA66y)^2);
            if Block.IP == 0
                ExportPtDiff(z).LHAdiffIP = sqrt((ExportAccuracy(z).LHAxb - ExportAccuracy(z).LHAIPx)^2 + (ExportAccuracy(z).LHAyb - ExportAccuracy(z).LHAIPy)^2);% IP
                ExportPtDiff(z).LFAdiffIP = sqrt((ExportAccuracy(z).LFAxb - ExportAccuracy(z).LFAIPx)^2 + (ExportAccuracy(z).LFAyb - ExportAccuracy(z).LFAIPy)^2);
            end
            if Block.HC == 0
                ExportPtDiff(z).LHAdiffHeelCent = sqrt((ExportAccuracy(z).LHAxb - ExportAccuracy(z).LHAHeelCentx)^2 + (ExportAccuracy(z).LHAyb - ExportAccuracy(z).LHAHeelCenty)^2); % HeelCent
                ExportPtDiff(z).LFAdiffHeelCent = sqrt((ExportAccuracy(z).LFAxb - ExportAccuracy(z).LFAHeelCentx)^2 + (ExportAccuracy(z).LFAyb - ExportAccuracy(z).LFAHeelCenty)^2);
            end
        end
        
        % RIGHT
        if isempty(DynamicPPTrials(z).RProg) == 0
            ExportProg(z).RProgDiffManual = abs(ExportProg(z).RProgC3D - ExportProg(z).RProgAngManual);
            ExportProg(z).RProgDiffGeneral = abs(ExportProg(z).RProgC3D - ExportProg(z).RProgAngGeneral);
            ExportProg(z).RProgDiffCoP = abs(ExportProg(z).RProgC3D - ExportProg(z).RProgAngCoP);
            ExportProg(z).RProgDiff66 = abs(ExportProg(z).RProgC3D - ExportProg(z).RProgAng66);
            if Block.IP == 0
                ExportProg(z).RProgDiffIP = abs(ExportProg(z).RProgC3D - ExportProg(z).RProgAngIP);
            end
            if Block.HC == 0
                ExportProg(z).RProgDiffHeelCent = abs(ExportProg(z).RProgC3D - ExportProg(z).RProgAngHeelCent);
            end
            
            ExportPtDiff(z).RHAdiffManual = sqrt((ExportAccuracy(z).RHAxb - ExportAccuracy(z).RHAManualx)^2 + (ExportAccuracy(z).RHAyb - ExportAccuracy(z).RHAManualy)^2); % manual
            ExportPtDiff(z).RFAdiffManual = sqrt((ExportAccuracy(z).RFAxb - ExportAccuracy(z).RFAManualx)^2 + (ExportAccuracy(z).RFAyb - ExportAccuracy(z).RFAManualy)^2);
            ExportPtDiff(z).RHAdiffGeneral = sqrt((ExportAccuracy(z).RHAxb - ExportAccuracy(z).RHAGeneralx)^2 + (ExportAccuracy(z).RHAyb - ExportAccuracy(z).RHAGeneraly)^2); % general
            ExportPtDiff(z).RFAdiffGeneral = sqrt((ExportAccuracy(z).RFAxb - ExportAccuracy(z).RFAGeneralx)^2 + (ExportAccuracy(z).RFAyb - ExportAccuracy(z).RFAGeneraly)^2);
            ExportPtDiff(z).RHAdiffCoP = sqrt((ExportAccuracy(z).RHAxb - ExportAccuracy(z).RHACoPx)^2 + (ExportAccuracy(z).RHAyb - ExportAccuracy(z).RHACoPy)^2); % CoP
            ExportPtDiff(z).RFAdiffCoP = sqrt((ExportAccuracy(z).RFAxb - ExportAccuracy(z).RFACoPx)^2 + (ExportAccuracy(z).RFAyb - ExportAccuracy(z).RFACoPy)^2);
            ExportPtDiff(z).RHAdiff66 = sqrt((ExportAccuracy(z).RHAxb - ExportAccuracy(z).RHA66x)^2 + (ExportAccuracy(z).RHAyb - ExportAccuracy(z).RHA66y)^2); % 66
            ExportPtDiff(z).RFAdiff66 = sqrt((ExportAccuracy(z).RFAxb - ExportAccuracy(z).RFA66x)^2 + (ExportAccuracy(z).RFAyb - ExportAccuracy(z).RFA66y)^2);
            if Block.IP == 0
                ExportPtDiff(z).RHAdiffIP = sqrt((ExportAccuracy(z).RHAxb - ExportAccuracy(z).RHAIPx)^2 + (ExportAccuracy(z).RHAyb - ExportAccuracy(z).RHAIPy)^2);% IP
                ExportPtDiff(z).RFAdiffIP = sqrt((ExportAccuracy(z).RFAxb - ExportAccuracy(z).RFAIPx)^2 + (ExportAccuracy(z).RFAyb - ExportAccuracy(z).RFAIPy)^2);
            end
            if Block.HC == 0
                ExportPtDiff(z).RHAdiffHeelCent = sqrt((ExportAccuracy(z).RHAxb - ExportAccuracy(z).RHAHeelCentx)^2 + (ExportAccuracy(z).RHAyb - ExportAccuracy(z).RHAHeelCenty)^2); % HeelCent
                ExportPtDiff(z).RFAdiffHeelCent = sqrt((ExportAccuracy(z).RFAxb - ExportAccuracy(z).RFAHeelCentx)^2 + (ExportAccuracy(z).RFAyb - ExportAccuracy(z).RFAHeelCenty)^2);
            end
        end
        
        %         ExportTS(z).CadenceC3D = C3Ddata(z).TempSpat(1).mean;
        %         ExportTS(z).WalkSpeedC3D = C3Ddata(z).TempSpat(2).mean;
        %         ExportTS(z).StepLengthC3D = C3Ddata(z).TempSpat(4).mean;
        %         ExportTS(z).StepWidthC3D = C3Ddata(z).TempSpat(5).mean;
        %         ExportTS(z).StanceC3D = C3Ddata(z).TempSpat(7).mean;
        %
        %         CoP
        %         ExportTS(z).CadencePP = TS.Cadence.Measures(z);
        %         ExportTS(z).WalkSpeedPP = TS.GaitSpeed.MpS.Measures(z);
        %         ExportTS(z).StepLengthPP = TS.StrideLength.Measures(z) / 100;
        %         ExportTS(z).StepWidthPP = TS.StepWidth.Measures(z) / 100;
        %         ExportTS(z).StancePP = TS.Stance.Measures(z);
        
        clearvars LProgC3D RProgC3D
    end
    
    % Exporting TS data and Foot prog angles for comparisons to MoCap
    ExportAccuracy = struct2table(ExportAccuracy);
    ExportProg = struct2table(ExportProg);
    ExportPtDiff = struct2table(ExportPtDiff);
    % ExportTS = [struct2table(ExportTS), TrimSpat(TStrials,:)];
    % Export = [ExportAccuracy, ExportTS];
    
    % Save Figure
    %     Names = strsplit(C3Ddata(1).FullFileName);
    %     FigFileName = strcat(Names{1}, 'Trials');
    %     savefig(Cont,FigFileName,'compact');
end

% Save Data for later use
% filename2 = strcat(filename(1:end-4), '_HalfData.mat');
% save(filename2);

%% Plot similar to WAAAG outputs with bullet graphs
if TempSpat == 1
    if strcmp(PPPlots.TempSpat, 'Yes')
        TSplots = figure;
        DisplayDim(1) = DisplayDim(1)+10;
        DisplayDim(2) = DisplayDim(2)+10;
        set(TSplots,  'Position', DisplayDim);
        % gait speed
        subplot(251);
        BulletGraph('V', 75, 90, TS.GaitSpeed.MpM.Avg/TSN.TimeDist(3)*100, 100, [0 140], rgb('RoyalBlue'),1, 'Yes', 'No');
        errorbar(1, TS.GaitSpeed.MpM.Avg/TSN.TimeDist(3)*100, TS.GaitSpeed.MpM.Std/TSN.TimeDist(3)*100,'k.')
        title('Gait Speed', 'FontSize',8);
        ylabel('% of Normal', 'FontSize',8);
        % cadence
        subplot(252);
        BulletGraph('V', 75, 90,  TS.Cadence.Avg/TSN.TimeDist(1)*100, 100, [0 140], rgb('RoyalBlue'), 1, 'Yes', 'No');
        errorbar(1, TS.Cadence.Avg/TSN.TimeDist(1)*100, TS.Cadence.Std/TSN.TimeDist(1)*100,'k.')
        title('Cadence', 'FontSize',8);
        % number of trials and steps!
        subplot(253);
        RGB = imread('Mushy_GUI_Foot_Logo.png');
        Bottom = 255*ones(400, 300,3);
        Bottom = uint8(Bottom);
        IM = vertcat(RGB,Bottom);
        imshow(IM, 'InitialMagnification', 'fit');
        txt = strcat('# Trials = ',num2str(length(TStrials)));
        text(10, 500, txt, 'FontSize', 7);
        
        NumLeftStrides = zeros(1,length(TStrials));
        NumRightStrides = zeros(1,length(TStrials));
        for j = 1:length(TStrials)
            NumLeftStrides(j) = DynamicPPTrials(TStrials(j)).NumLeft;
            NumRightStrides(j) = DynamicPPTrials(TStrials(j)).NumRight;
        end
        txt = strcat('# L Steps = ',num2str(sum(NumLeftStrides)));...
            text(10, 550, txt, 'FontSize', 7);
        txt = strcat('# R Steps = ',num2str(sum(NumRightStrides)));
        text(10, 600, txt, 'FontSize', 7);
        text(25, 420, 'IMPRESSIONS', 'FontSize', 7);
        % step length
        subplot(254);
        BulletGraph('V', 75, 90,  TS.StrideLength.Avg/TSN.TimeDist(4)*100, 100, [0 140], rgb('Crimson'), 1, 'Yes', 'No');
        errorbar(1, TS.StrideLength.Avg/TSN.TimeDist(4)*100, TS.StrideLength.Std/TSN.TimeDist(4)*100,'k.')
        title('Stride Length', 'FontSize',8);
        % step width
        subplot(255);
        BulletGraph('ReverseV', 125, 110,  TS.StepWidth.Avg/TSN.TimeDist(6)*100,100, [0 300], rgb('Crimson'), 1, 'Yes', 'No');
        errorbar(1, TS.StepWidth.Avg/TSN.TimeDist(6)*100, TS.StepWidth.Std/TSN.TimeDist(6)*100,'k.')
        title('Step Width', 'FontSize',8);
        % Stance Phase
        subplot(256);
        if isnan(TS.Stance.Avg)
            text(0.44, 90, 'Not enough continuous', 'FontSize',6)
            text(0.44, 85, 'steps for timing calcs.', 'FontSize',6)
        end
        BulletGraph('V', 50, 58, TS.Stance.Avg*100, 62, [0 100], rgb('ForestGreen'), 1, 'Yes', 'No');
        errorbar(1, TS.Stance.Avg*100, TS.Stance.Std*100,'k.')
        title('Stance', 'FontSize',8);
        ylabel('% of Gait Cycle', 'FontSize',8);
        % Swing Phase
        subplot(257);
        BulletGraph('V', 20, 30,  TS.Swing.Avg*100, 38, [0 100], rgb('ForestGreen'), 1, 'Yes', 'No');
        errorbar(1, TS.Swing.Avg*100, TS.Swing.Std*100,'k.')
        title('Swing', 'FontSize',8);
        % initial double support
        subplot(258);
        BulletGraph('ReverseV', 25, 15,  TS.IDS.Avg*100, 12, [0 50], rgb('Violet'), 1, 'Yes', 'No');
        errorbar(1, TS.IDS.Avg*100, TS.IDS.Std*100,'k.')
        title('1DS', 'FontSize',8);
        %  single limb support
        subplot(259);
        BulletGraph('V', 25, 30,  TS.SS.Avg*100, 38, [0 50], rgb('Violet'), 1, 'Yes', 'No');
        errorbar(1, TS.SS.Avg*100, TS.SS.Std*100,'k.')
        title('SLS', 'FontSize',8);
        %  secondary double support
        subplot(2,5,10);
        BulletGraph('ReverseV', 25, 15,  TS.SDS.Avg*100, 12, [0 50], rgb('Violet'), 1, 'Yes', 'No');
        errorbar(1, TS.SDS.Avg*100, TS.SDS.Std*100,'k.')
        title('2DS', 'FontSize',8);
        
        for i = [1 2 4 5 6 7 8 9 10]
            subplot(2,5,i);
            ax = gca;
            ax.FontSize = 7;
        end
        clearvars StepWidth StrideLength
        
        subplotsqueeze(TSplots, 1.02);
        % export figure to excel
        if strcmp(PPSettings.ExportReport, 'Yes') == 1
            % write to output file
            xlsPasteTo('ImpressionsOutput.xlsx','Sheet1',680,580, 'O17');
        end
    end
end

%% Extract Left and Right times from Regioned data
for j = 1:length(Selection)
    % LEFT foot timing
    for i = 1:length(Regions{1,Selection(j)})
        if strcmp(Regions{1,Selection(j)}(i).Side, 'Left') == 1
            DynamicPPTrials(Selection(j)).Strike.Left(i) = Regions{1,Selection(j)}(i).Strike;
            if isempty(Regions{1,Selection(j)}(i).Off) == 1
                DynamicPPTrials(Selection(j)).Off.Left(i) = 248;
            else
                DynamicPPTrials(Selection(j)).Off.Left(i) = Regions{1,Selection(j)}(i).Off;
            end
        end
    end
    % RIGHT foot timing
    for i = 1:length(Regions{1,Selection(j)})
        if strcmp(Regions{1,Selection(j)}(i).Side, 'Right') == 1
            DynamicPPTrials(Selection(j)).Strike.Right(i) = Regions{1,Selection(j)}(i).Strike;
            if isempty(Regions{1,Selection(j)}(i).Off) == 1
                DynamicPPTrials(Selection(j)).Off.Right(i) = 248;
            else
                DynamicPPTrials(Selection(j)).Off.Right(i) = Regions{1,Selection(j)}(i).Off;
            end
        end
    end
    % delete duplicates
    if isfield(DynamicPPTrials(Selection(j)).Strike, 'Left')
        for i = 1:length(DynamicPPTrials(Selection(j)).Strike.Left)-1
            if DynamicPPTrials(Selection(j)).Strike.Left(i) == DynamicPPTrials(Selection(j)).Strike.Left(i+1)
                DynamicPPTrials(Selection(j)).Strike.Left(i+1) = [];
                break
            end
        end
        for i = 1:length(DynamicPPTrials(Selection(j)).Off.Left)-1
            if DynamicPPTrials(Selection(j)).Off.Left(i) == DynamicPPTrials(Selection(j)).Off.Left(i+1)
                DynamicPPTrials(Selection(j)).Off.Left(i+1) = [];
                break
            end
        end
        % delete zeros
        for i = 1:length(DynamicPPTrials(Selection(j)).Strike.Left) % if an empty matrix is input, delete it
            if DynamicPPTrials(Selection(j)).Strike.Left(i) == 0
                DynamicPPTrials(Selection(j)).Strike.Left(i) = [];
                DynamicPPTrials(Selection(j)).Off.Left(i) = [];
                break
            end
        end
    end
    if isfield(DynamicPPTrials(Selection(j)).Strike, 'Right')
        for i = 1:length(DynamicPPTrials(Selection(j)).Strike.Right)-1
            if DynamicPPTrials(Selection(j)).Strike.Right(i) == DynamicPPTrials(Selection(j)).Strike.Right(i+1)
                DynamicPPTrials(Selection(j)).Strike.Right(i+1) = [];
                break
            end
        end
        for i = 1:length(DynamicPPTrials(Selection(j)).Off.Right)-1
            if DynamicPPTrials(Selection(j)).Off.Right(i) == DynamicPPTrials(Selection(j)).Off.Right(i+1)
                DynamicPPTrials(Selection(j)).Off.Right(i+1) = [];
                break
            end
        end
        % delete zeros
        for i = 1:length(DynamicPPTrials(Selection(j)).Strike.Right) % if an empty matrix is input, delete it
            if DynamicPPTrials(Selection(j)).Strike.Right(i) == 0
                DynamicPPTrials(Selection(j)).Strike.Right(i) = [];
                DynamicPPTrials(Selection(j)).Off.Right(i) = [];
                break
            end
        end
    end
end

%% Generate a L & R plots with COP data
if strcmp(PPPlots.CloseUp, 'Yes')
    clc; disp('Generating close-up images');
%     fh = findobj( 'Type', 'Figure', 'Name', 'COPplots' );
%     % determine which trial they would like to use
%     if length(Selection) > 1
%         Trial2Cont = listdlg( 'PromptString','Select Trial to plot plantar pressures up close','ListString',ListStr);
%     else
%         Trial2Cont = Selection;
%     end
%     Trial2Cont(2) = Trial2Cont;
%     if DynamicPPTrials(Trial2Cont(1)).NumLeft == 0
%         Trial2Cont(1) = listdlg( 'PromptString','No eligible LEFTs, select another trial','ListString',ListStr);
%     end
%     if DynamicPPTrials(Trial2Cont(2)).NumRight == 0
%         Trial2Cont(2) = listdlg( 'PromptString','No eligible RIGHTs, select another trial','ListString',ListStr);
%     end
%     % define steps to pull out
%     % LEFT side
%     if DynamicPPTrials(Trial2Cont(1)).NumLeft == 1
%         Lcont = 1;
%     else
%         Question = 'Which LEFT foot strike would you like to plot?';
%         NumLeftCont = questdlg(Question, 'LeftFootPlot','1','2','3','1');
%         if NumLeftCont == '1'
%             Lcont = 1;
%         end
%         if NumLeftCont == '2'
%             Lcont = 2;
%         end
%         if NumLeftCont == '3'
%             Lcont = 3;
%         end
%     end
    Counter = 0;
    for i = 1:length(Regions{Trial2Cont(1)})
        if strcmp(Regions{Trial2Cont(1)}(i).Side,'Left')
            Counter = Counter + 1;
            if Counter == Lcont
                LcontReg = i;
                break
            end
        end
    end
    
    % RIGHT side
    if DynamicPPTrials(Trial2Cont(2)).NumRight == 1
        Rcont = 1;
    else
        Question = 'Which RIGHT foot strike would you like to plot?';
        NumRightCont = questdlg(Question, 'RightFootPlot','1','2','3','1');
        if NumRightCont == '1'
            Rcont = 1;
        end
        if NumRightCont == '2'
            Rcont = 2;
        end
        if NumRightCont == '3'
            Rcont = 3;
        end
    end
    Counter = 0;
    for i = 1:length(Regions{Trial2Cont(2)})
        if strcmp(Regions{Trial2Cont(2)}(i).Side,'Right')
            Counter = Counter + 1;
            if Counter == Rcont
                RcontReg = i;
                break
            end
        end
    end
    clearvars NumLeftCont NumRightCont question Question a i j q h1 h2 FPchange EditFPAngles Counter
    
    %%  Create COP figure
    COPplots = figure;
    DisplayDim(1) = DisplayDim(1)+10;
    DisplayDim(2) = DisplayDim(2)+10;
    set(COPplots,  'Position', DisplayDim);
    % Left Foot Plot
    subplot(1,2,1);
    contour(Regions{Trial2Cont(1)}(LcontReg).StepSum, 300, 'LineWidth', 1);
    colormap(jet);axis equal; hold on;
    plot(Regions{Trial2Cont(1)}(LcontReg).StepCoP(:,1),Regions{Trial2Cont(1)}(LcontReg).StepCoP(:,2), 'k.', 'MarkerSize',10);
    title('Left Foot Pressure Map');
    set(gca, 'YTick',[]);
    set(gca, 'XTick',[]);
    
    % Right Foot Plot
    subplot(1,2,2);
    contour(Regions{Trial2Cont(2)}(RcontReg).StepSum, 300, 'LineWidth', 1); 
    axis equal; hold on;
    plot(Regions{Trial2Cont(2)}(RcontReg).StepCoP(:,1),Regions{Trial2Cont(2)}(RcontReg).StepCoP(:,2), 'k.', 'MarkerSize',10);
    title('Right Foot Pressure Map');
    set(gca, 'YTick',[]);
    set(gca, 'XTick',[]);
    
    clearvars Del SumDel LogDel
    
    %     %% Load COP norms for plot presentation
    %     if strcmp(PPLots.COPNorms,'Yes')
    %         load('COPnormsData.mat');
    %         % Dynamic Pedobarography for Children: Use of the Center of Pressure Progression
    %         % Jameson et al. J Pediatry Orthop 2008; 28: 254-258.
    %
    %         % apply COP norm zones to pressure maps
    %         % LEFT
    %         subplot(121);
    %         % flip COP norms for left side only
    %         Flipx1 = [1 - COPnormsData.SD_1.x_med, 1 - COPnormsData.SD_1.x_lat];
    %         Flipy1 = [1 - COPnormsData.SD_1.y_med, 1 - COPnormsData.SD_1.y_lat];
    %         Flipx2 = [1 - COPnormsData.SD_2.x_med, 1 - COPnormsData.SD_2.x_lat];
    %         Flipy2 = [1 - COPnormsData.SD_2.y_med, 1 - COPnormsData.SD_2.y_lat];
    %         % define angle (a) to rotate COP norms
    %         a = DynamicPPTrials(Trial2Cont).LProg.AngActual{Lcont};
    %         % define scaling factor using foot length
    %         scale = 0.9*sqrt( (abs(LhtCont{Lcont}(1,1) - LhtCont{Lcont}(2,1)))^2 + (abs(Lht.GeneralCont{Lcont}(1,2) - LhtCont{Lcont}(2,2)))^2);
    %         % determine offsets for plotting norms using a combination of scaling,
    %         % centroids, and the box adjust. Y offset us just the height difference
    %         xOffSet = 0;
    %         yOffSet = min(LhtCont{Lcont}(:,2));
    %         % multiply COP norms by the scale and apply offsets
    %         x1 = [scale*Flipx1(:,1) + xOffSet, scale*Flipx1(:,2) + xOffSet];
    %         y1 = [scale*COPnormsData.SD_1.y_med + yOffSet, scale*COPnormsData.SD_1.y_lat+ yOffSet];
    %         x2 = [scale*Flipx2(:,1) + xOffSet, scale*Flipx2(:,2) + xOffSet];
    %         y2 = [scale*COPnormsData.SD_2.y_med + yOffSet, scale*COPnormsData.SD_2.y_lat + yOffSet];
    %         % plot and rotate the COP norms
    %         h1 =  plot(x1,y1,'k--'); hold on; % SD_1 lines
    %         rotate(h1, [0 0 1], a+10);
    %         h2 = plot(x2,y2,'r--'); axis equal; %SD_2 lines
    %         rotate(h2, [0 0 1], a+10);
    %         % RIGHT
    %         subplot(122);
    %         % define angle (a) to rotate COP norms
    %         a = DynamicPPTrials(Trial2Cont).RProg.AngActual{Rcont};
    %         % define scaling factor using foot length
    %         scale = 0.9*sqrt( (abs(RhtCont{Rcont}(1,1) - RhtCont{Rcont}(2,1)))^2 + (abs(RhtCont{Rcont}(1,2) - RhtCont{Rcont}(2,2)))^2);
    %         % determine offsets for plotting norms using a combination of scaling,
    %         % centroids, and the box adjust. Y offset us just the height difference
    %         xOffSet = ((scale*Flipx1(1,1) + scale*Flipx1(1,2)) /2) - (DynamicPPTrials(Trial2Cont).RightFeet(Rcont).Centroid(1) - RFootBox{1}(1,1));
    %         yOffSet = min(RhtCont{Rcont}(:,2)); % Y target
    %         % multiply COP norms by the scale and apply offsets
    %         x1 = [scale*COPnormsData.SD_1.x_med - xOffSet, scale*COPnormsData.SD_1.x_lat - xOffSet];
    %         y1 = [scale*COPnormsData.SD_1.y_med + yOffSet, scale*COPnormsData.SD_1.y_lat+ yOffSet];
    %         x2 = [scale*COPnormsData.SD_2.x_med - xOffSet, scale*COPnormsData.SD_2.x_lat - xOffSet];
    %         y2 = [scale*COPnormsData.SD_2.y_med + yOffSet, scale*COPnormsData.SD_2.y_lat + yOffSet];
    %         % plot and rotate the COP norms
    %         h1 =  plot(x1,y1,'k--'); hold on; % SD_1 lines
    %         rotate(h1, [0 0 1], a-10);
    %         h2 = plot(x2,y2,'r--'); axis equal; %SD_2 lines
    %         rotate(h2, [0 0 1], a-10);
    %     end
    subplotsqueeze(COPplots,1.1);
    saveas(COPplots,strcat(folder,'\','CloseUps.png'));
    % export figure to excel
    if strcmp(PPSettings.ExportReport, 'Yes') == 1
        xlsPasteTo('ImpressionsOutput.xlsx','Sheet1',700,580, 'U16');
    end
    clearvars x1 x2 y1  y2 Pt1 Pt2 Flipx1 Flipx2 Flipy1 Flipy2 LoadNew
end

%% Create foot zone masks
[m,n] = size(PPTrials(1).SumTM);
for j = 1:length(Selection)
    % LEFT Foot
    for i = 1:DynamicPPTrials(Selection(j)).NumLeft
        % Create Masks
        % Heel
        DynamicPPTrials(Selection(j)).Mask.L(i).LatHeel = poly2mask(DynamicPPTrials(Selection(j)).Zones.L(i).LatHeel(:,1), DynamicPPTrials(Selection(j)).Zones.L(i).LatHeel(:,2),m,n);
        DynamicPPTrials(Selection(j)).Mask.L(i).MedHeel = poly2mask(DynamicPPTrials(Selection(j)).Zones.L(i).MedHeel(:,1), DynamicPPTrials(Selection(j)).Zones.L(i).MedHeel(:,2),m,n);
        % Arch
        DynamicPPTrials(Selection(j)).Mask.L(i).LatArch = poly2mask(DynamicPPTrials(Selection(j)).Zones.L(i).LatArch(:,1), DynamicPPTrials(Selection(j)).Zones.L(i).LatArch(:,2),m,n);
        DynamicPPTrials(Selection(j)).Mask.L(i).MedArch = poly2mask(DynamicPPTrials(Selection(j)).Zones.L(i).MedArch(:,1), DynamicPPTrials(Selection(j)).Zones.L(i).MedArch(:,2),m,n);
        % Forefoot
        DynamicPPTrials(Selection(j)).Mask.L(i).LatFore = poly2mask(DynamicPPTrials(Selection(j)).Zones.L(i).LatFore(:,1), DynamicPPTrials(Selection(j)).Zones.L(i).LatFore(:,2),m,n);
        DynamicPPTrials(Selection(j)).Mask.L(i).MedFore = poly2mask(DynamicPPTrials(Selection(j)).Zones.L(i).MedFore(:,1), DynamicPPTrials(Selection(j)).Zones.L(i).MedFore(:,2),m,n);
        % Whole Foot
        DynamicPPTrials(Selection(j)).Mask.L(i).Whole = poly2mask(DynamicPPTrials(Selection(j)).Zones.L(i).Whole(:,1), DynamicPPTrials(Selection(j)).Zones.L(i).Whole(:,2),m,n);
    end
    % RIGHT Foot
    for i = 1:DynamicPPTrials(Selection(j)).NumRight
        % Heel
        DynamicPPTrials(Selection(j)).Mask.R(i).LatHeel = poly2mask(DynamicPPTrials(Selection(j)).Zones.R(i).LatHeel(:,1), DynamicPPTrials(Selection(j)).Zones.R(i).LatHeel(:,2),m,n);
        DynamicPPTrials(Selection(j)).Mask.R(i).MedHeel = poly2mask(DynamicPPTrials(Selection(j)).Zones.R(i).MedHeel(:,1), DynamicPPTrials(Selection(j)).Zones.R(i).MedHeel(:,2),m,n);
        % Arch
        DynamicPPTrials(Selection(j)).Mask.R(i).LatArch = poly2mask(DynamicPPTrials(Selection(j)).Zones.R(i).LatArch(:,1), DynamicPPTrials(Selection(j)).Zones.R(i).LatArch(:,2),m,n);
        DynamicPPTrials(Selection(j)).Mask.R(i).MedArch = poly2mask(DynamicPPTrials(Selection(j)).Zones.R(i).MedArch(:,1), DynamicPPTrials(Selection(j)).Zones.R(i).MedArch(:,2),m,n);
        % Forefoot
        DynamicPPTrials(Selection(j)).Mask.R(i).LatFore = poly2mask(DynamicPPTrials(Selection(j)).Zones.R(i).LatFore(:,1), DynamicPPTrials(Selection(j)).Zones.R(i).LatFore(:,2),m,n);
        DynamicPPTrials(Selection(j)).Mask.R(i).MedFore = poly2mask(DynamicPPTrials(Selection(j)).Zones.R(i).MedFore(:,1), DynamicPPTrials(Selection(j)).Zones.R(i).MedFore(:,2),m,n);
        % Whole Foot
        DynamicPPTrials(Selection(j)).Mask.R(i).Whole = poly2mask(DynamicPPTrials(Selection(j)).Zones.R(i).Whole(:,1), DynamicPPTrials(Selection(j)).Zones.R(i).Whole(:,2),m,n);
    end
end

%% Calculate Time Normalized Foot Pressures, Areas, and Forces
% clc; disp('Calculating foot/floor interactions (pressures, forces, & areas)');
% % Create time series Logical matrix for area calcualtions
% for k = 1:length(DynamicPPTrials)
%     DynamicPPTrials(k).TMLog = logical(DynamicPPTrials(k).TM > 0);
%     %DynamicPPTrials(k).TMLog(DynamicPPTrials(k).TMLog > 0 == 1);
% end
% % correct pressure, areas, and force units between the two types of inputs
% Subject.Newtons = Subject.Mass * 9.81;
% % double check that there are no 0s in the foot strike and foot off times
% for j = 1:length(Selection)
%     % delete zeros
%     if isfield(DynamicPPTrials(Selection(j)).Strike, 'Left')
%         for i = 1:length(DynamicPPTrials(Selection(j)).Strike.Left) % if an empty matrix is input, delete it
%             if DynamicPPTrials(Selection(j)).Strike.Left(i) == 0
%                 DynamicPPTrials(Selection(j)).Strike.Left(i) = [];
%                 DynamicPPTrials(Selection(j)).Off.Left(i) = [];
%                 break
%             end
%         end
%     end
%     if isfield(DynamicPPTrials(Selection(j)).Strike, 'Right')
%         for i = 1:length(DynamicPPTrials(Selection(j)).Strike.Right) % if an empty matrix is input, delete it
%             if DynamicPPTrials(Selection(j)).Strike.Right(i) == 0
%                 DynamicPPTrials(Selection(j)).Strike.Right(i) = [];
%                 DynamicPPTrials(Selection(j)).Off.Right(i) = [];
%                 break
%             end
%         end
%     end
% end
% % Time normalization is by percent of stance
% if strcmp(PPSettings.PPMatType,'RSScan') == 1 % trial data is from a RSScan mat
%     % RSScan exports in newtons (N)
%     for j = 1:length(Selection)
%         % LEFT
%         % stop loop from running through more times than there are foot
%         % offs. If there isn't a clean foot off, its not a valid full step
%         if DynamicPPTrials(Selection(j)).NumLeft > 0
%             if length(DynamicPPTrials(Selection(j)).Off.Left) < length(DynamicPPTrials(Selection(j)).Strike.Left)
%                 DynamicPPTrials(Selection(j)).NumFullStepsL = length(DynamicPPTrials(Selection(j)).Off.Left);
%             else
%                 DynamicPPTrials(Selection(j)).NumFullStepsL = DynamicPPTrials(Selection(j)).NumLeft;
%             end
%             for i = 1:DynamicPPTrials(Selection(j)).NumFullStepsL
%                 % loop through the stride and add up the sum of pressures per unit time
%                 for stride = DynamicPPTrials(Selection(j)).Strike.Left(i):DynamicPPTrials(Selection(j)).Off.Left(i)
%                     ind = stride - DynamicPPTrials(Selection(j)).Strike.Left(i) +1;
%                     % Heel
%                     DynamicPPTrials(Selection(j)).Force.Capture.L.LatHeel{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.L(i).LatHeel .* DynamicPPTrials(Selection(j)).TM(:,:,stride))); % force in N
%                     DynamicPPTrials(Selection(j)).Area.Capture.L.LatHeel{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.L(i).LatHeel .* DynamicPPTrials(Selection(j)).TMLog(:,:,stride))) * Adj.Area;  % area adjusted to cm^2
%                     DynamicPPTrials(Selection(j)).Pressure.Capture.L.LatHeel{i}(ind) = DynamicPPTrials(Selection(j)).Force.Capture.L.LatHeel{i}(ind) / DynamicPPTrials(Selection(j)).Area.Capture.L.LatHeel{i}(ind); % pressure in Pa
%                     DynamicPPTrials(Selection(j)).Force.Capture.L.MedHeel{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.L(i).MedHeel .* DynamicPPTrials(Selection(j)).TM(:,:,stride))); % force in N
%                     DynamicPPTrials(Selection(j)).Area.Capture.L.MedHeel{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.L(i).MedHeel .* DynamicPPTrials(Selection(j)).TMLog(:,:,stride))) * Adj.Area;  % area adjusted to cm^2
%                     DynamicPPTrials(Selection(j)).Pressure.Capture.L.MedHeel{i}(ind) = DynamicPPTrials(Selection(j)).Force.Capture.L.MedHeel{i}(ind) / DynamicPPTrials(Selection(j)).Area.Capture.L.MedHeel{i}(ind); % pressure in Pa
%                     % Arch
%                     DynamicPPTrials(Selection(j)).Force.Capture.L.LatArch{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.L(i).LatArch .* DynamicPPTrials(Selection(j)).TM(:,:,stride))); % force in N
%                     DynamicPPTrials(Selection(j)).Area.Capture.L.LatArch{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.L(i).LatArch .* DynamicPPTrials(Selection(j)).TMLog(:,:,stride))) * Adj.Area;  % area adjusted to cm^2
%                     DynamicPPTrials(Selection(j)).Pressure.Capture.L.LatArch{i}(ind) = DynamicPPTrials(Selection(j)).Force.Capture.L.LatArch{i}(ind) / DynamicPPTrials(Selection(j)).Area.Capture.L.LatArch{i}(ind); % pressure in Pa
%                     DynamicPPTrials(Selection(j)).Force.Capture.L.MedArch{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.L(i).MedArch .* DynamicPPTrials(Selection(j)).TM(:,:,stride))); % force in N
%                     DynamicPPTrials(Selection(j)).Area.Capture.L.MedArch{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.L(i).MedArch .* DynamicPPTrials(Selection(j)).TMLog(:,:,stride))) * Adj.Area;  % area adjusted to cm^2
%                     DynamicPPTrials(Selection(j)).Pressure.Capture.L.MedArch{i}(ind) = DynamicPPTrials(Selection(j)).Force.Capture.L.MedArch{i}(ind) / DynamicPPTrials(Selection(j)).Area.Capture.L.MedArch{i}(ind); % pressure in Pa
%                     % Forefoot
%                     DynamicPPTrials(Selection(j)).Force.Capture.L.LatFore{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.L(i).LatFore .* DynamicPPTrials(Selection(j)).TM(:,:,stride))); % force in N
%                     DynamicPPTrials(Selection(j)).Area.Capture.L.LatFore{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.L(i).LatFore .* DynamicPPTrials(Selection(j)).TMLog(:,:,stride))) * Adj.Area;  % area adjusted to cm^2
%                     DynamicPPTrials(Selection(j)).Pressure.Capture.L.LatFore{i}(ind) = DynamicPPTrials(Selection(j)).Force.Capture.L.LatFore{i}(ind) / DynamicPPTrials(Selection(j)).Area.Capture.L.LatFore{i}(ind); % pressure in Pa
%                     DynamicPPTrials(Selection(j)).Force.Capture.L.MedFore{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.L(i).MedFore .* DynamicPPTrials(Selection(j)).TM(:,:,stride))); % force in N
%                     DynamicPPTrials(Selection(j)).Area.Capture.L.MedFore{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.L(i).MedFore .* DynamicPPTrials(Selection(j)).TMLog(:,:,stride))) * Adj.Area;  % area adjusted to cm^2
%                     DynamicPPTrials(Selection(j)).Pressure.Capture.L.MedFore{i}(ind) = DynamicPPTrials(Selection(j)).Force.Capture.L.MedFore{i}(ind) / DynamicPPTrials(Selection(j)).Area.Capture.L.MedFore{i}(ind); % pressure in Pa
%                     % correct any NaNs
%                     if  isnan(DynamicPPTrials(Selection(j)).Pressure.Capture.L.LatHeel{i}(ind)) == 1
%                         DynamicPPTrials(Selection(j)).Pressure.Capture.L.LatHeel{i}(ind) = 0;
%                     end
%                     if  isnan(DynamicPPTrials(Selection(j)).Pressure.Capture.L.MedHeel{i}(ind)) == 1
%                         DynamicPPTrials(Selection(j)).Pressure.Capture.L.MedHeel{i}(ind) = 0;
%                     end
%                     if  isnan(DynamicPPTrials(Selection(j)).Pressure.Capture.L.LatArch{i}(ind)) == 1
%                         DynamicPPTrials(Selection(j)).Pressure.Capture.L.LatArch{i}(ind) = 0;
%                     end
%                     if  isnan(DynamicPPTrials(Selection(j)).Pressure.Capture.L.MedArch{i}(ind)) == 1
%                         DynamicPPTrials(Selection(j)).Pressure.Capture.L.MedArch{i}(ind) = 0;
%                     end
%                     if  isnan(DynamicPPTrials(Selection(j)).Pressure.Capture.L.LatFore{i}(ind)) == 1
%                         DynamicPPTrials(Selection(j)).Pressure.Capture.L.LatFore{i}(ind) = 0;
%                     end
%                     if  isnan(DynamicPPTrials(Selection(j)).Pressure.Capture.L.MedFore{i}(ind)) == 1
%                         DynamicPPTrials(Selection(j)).Pressure.Capture.L.MedFore{i}(ind) = 0;
%                     end
%                 end
%             end
%         end
%         % RIGHT Foot
%         % stop loop from running through more times than there are foot
%         % offs. If there isn't a clean foot off, its not a valid full step
%         if DynamicPPTrials(Selection(j)).NumRight > 0
%             if length(DynamicPPTrials(Selection(j)).Off.Right) < length(DynamicPPTrials(Selection(j)).Strike.Right)
%                 DynamicPPTrials(Selection(j)).NumFullStepsR = length(DynamicPPTrials(Selection(j)).Off.Right);
%             else
%                 DynamicPPTrials(Selection(j)).NumFullStepsR = DynamicPPTrials(Selection(j)).NumRight;
%             end
%             for i = 1:DynamicPPTrials(Selection(j)).NumFullStepsR
%                 for stride = DynamicPPTrials(Selection(j)).Strike.Right(i):DynamicPPTrials(Selection(j)).Off.Right(i)
%                     ind = stride - DynamicPPTrials(Selection(j)).Strike.Right(i) +1;
%                     % Heel
%                     DynamicPPTrials(Selection(j)).Force.Capture.R.LatHeel{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.R(i).LatHeel .* DynamicPPTrials(Selection(j)).TM(:,:,stride))); % force in N
%                     DynamicPPTrials(Selection(j)).Area.Capture.R.LatHeel{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.R(i).LatHeel .* DynamicPPTrials(Selection(j)).TMLog(:,:,stride))) * Adj.Area;  % area adjusted to cm^2
%                     DynamicPPTrials(Selection(j)).Pressure.Capture.R.LatHeel{i}(ind) = DynamicPPTrials(Selection(j)).Force.Capture.R.LatHeel{i}(ind) / DynamicPPTrials(Selection(j)).Area.Capture.R.LatHeel{i}(ind); % pressure in Pa
%                     DynamicPPTrials(Selection(j)).Force.Capture.R.MedHeel{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.R(i).MedHeel .* DynamicPPTrials(Selection(j)).TM(:,:,stride))); % force in N
%                     DynamicPPTrials(Selection(j)).Area.Capture.R.MedHeel{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.R(i).MedHeel .* DynamicPPTrials(Selection(j)).TMLog(:,:,stride))) * Adj.Area;  % area adjusted to cm^2
%                     DynamicPPTrials(Selection(j)).Pressure.Capture.R.MedHeel{i}(ind) = DynamicPPTrials(Selection(j)).Force.Capture.R.MedHeel{i}(ind) / DynamicPPTrials(Selection(j)).Area.Capture.R.MedHeel{i}(ind); % pressure in Pa
%                     % Arch
%                     DynamicPPTrials(Selection(j)).Force.Capture.R.LatArch{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.R(i).LatArch .* DynamicPPTrials(Selection(j)).TM(:,:,stride))); % force in N
%                     DynamicPPTrials(Selection(j)).Area.Capture.R.LatArch{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.R(i).LatArch .* DynamicPPTrials(Selection(j)).TMLog(:,:,stride))) * Adj.Area;  % area adjusted to cm^2
%                     DynamicPPTrials(Selection(j)).Pressure.Capture.R.LatArch{i}(ind) = DynamicPPTrials(Selection(j)).Force.Capture.R.LatArch{i}(ind) / DynamicPPTrials(Selection(j)).Area.Capture.R.LatArch{i}(ind); % pressure in Pa
%                     DynamicPPTrials(Selection(j)).Force.Capture.R.MedArch{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.R(i).MedArch .* DynamicPPTrials(Selection(j)).TM(:,:,stride))); % force in N
%                     DynamicPPTrials(Selection(j)).Area.Capture.R.MedArch{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.R(i).MedArch .* DynamicPPTrials(Selection(j)).TMLog(:,:,stride))) * Adj.Area;  % area adjusted to cm^2
%                     DynamicPPTrials(Selection(j)).Pressure.Capture.R.MedArch{i}(ind) = DynamicPPTrials(Selection(j)).Force.Capture.R.MedArch{i}(ind) / DynamicPPTrials(Selection(j)).Area.Capture.R.MedArch{i}(ind); % pressure in Pa
%                     % Forefoot
%                     DynamicPPTrials(Selection(j)).Force.Capture.R.LatFore{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.R(i).LatFore .* DynamicPPTrials(Selection(j)).TM(:,:,stride))); % force in N
%                     DynamicPPTrials(Selection(j)).Area.Capture.R.LatFore{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.R(i).LatFore .* DynamicPPTrials(Selection(j)).TMLog(:,:,stride))) * Adj.Area;  % area adjusted to cm^2
%                     DynamicPPTrials(Selection(j)).Pressure.Capture.R.LatFore{i}(ind) = DynamicPPTrials(Selection(j)).Force.Capture.R.LatFore{i}(ind) / DynamicPPTrials(Selection(j)).Area.Capture.R.LatFore{i}(ind); % pressure in Pa
%                     DynamicPPTrials(Selection(j)).Force.Capture.R.MedFore{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.R(i).MedFore .* DynamicPPTrials(Selection(j)).TM(:,:,stride))); % force in N
%                     DynamicPPTrials(Selection(j)).Area.Capture.R.MedFore{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.R(i).MedFore .* DynamicPPTrials(Selection(j)).TMLog(:,:,stride))) * Adj.Area;  % area adjusted to cm^2
%                     DynamicPPTrials(Selection(j)).Pressure.Capture.R.MedFore{i}(ind) = DynamicPPTrials(Selection(j)).Force.Capture.R.MedFore{i}(ind) / DynamicPPTrials(Selection(j)).Area.Capture.R.MedFore{i}(ind); % pressure in Pa
%                     % correct any NaNs
%                     if  isnan(DynamicPPTrials(Selection(j)).Pressure.Capture.R.LatHeel{i}(ind)) == 1
%                         DynamicPPTrials(Selection(j)).Pressure.Capture.R.LatHeel{i}(ind) = 0;
%                     end
%                     if isnan(DynamicPPTrials(Selection(j)).Pressure.Capture.R.MedHeel{i}(ind)) == 1
%                         DynamicPPTrials(Selection(j)).Pressure.Capture.R.MedHeel{i}(ind) = 0;
%                     end
%                     if  isnan(DynamicPPTrials(Selection(j)).Pressure.Capture.R.LatArch{i}(ind)) == 1
%                         DynamicPPTrials(Selection(j)).Pressure.Capture.R.LatArch{i}(ind) = 0;
%                     end
%                     if  isnan(DynamicPPTrials(Selection(j)).Pressure.Capture.R.MedArch{i}(ind)) == 1
%                         DynamicPPTrials(Selection(j)).Pressure.Capture.R.MedArch{i}(ind) = 0;
%                     end
%                     if  isnan(DynamicPPTrials(Selection(j)).Pressure.Capture.R.LatFore{i}(ind)) == 1
%                         DynamicPPTrials(Selection(j)).Pressure.Capture.R.LatFore{i}(ind) = 0;
%                     end
%                     if  isnan(DynamicPPTrials(Selection(j)).Pressure.Capture.R.MedFore{i}(ind)) == 1
%                         DynamicPPTrials(Selection(j)).Pressure.Capture.R.MedFore{i}(ind) = 0;
%                     end
%                 end
%             end
%         end
%     end
% elseif strcmp(PPSettings.PPMatType, 'Novel') % trial data is from Novel mat
%     % Novel exports the data in units of pressure
%     for j = 1:length(Selection)
%         % LEFT
%         % stop loop from running through more times than there are foot offs. If there isn't a clean foot off, its not a valid full step
%         if isfield(DynamicPPTrials(Selection(j)).Off, 'Left')
%             if length(DynamicPPTrials(Selection(j)).Off.Left) < length(DynamicPPTrials(Selection(j)).Strike.Left)
%                 DynamicPPTrials(Selection(j)).NumFullStepsL = length(DynamicPPTrials(Selection(j)).Off.Left);
%             else
%                 DynamicPPTrials(Selection(j)).NumFullStepsL = DynamicPPTrials(Selection(j)).NumLeft;
%             end
%         end
%         if isfield(DynamicPPTrials(Selection(j)), 'NumFullStepsL')
%             for i = 1:DynamicPPTrials(Selection(j)).NumFullStepsL
%                 if DynamicPPTrials(Selection(j)).Strike.Left(i) == 0 && DynamicPPTrials(Selection(j)).Off.Left(i) == 0
%                     DynamicPPTrials(Selection(j)).Off.Left(i) = [];
%                     DynamicPPTrials(Selection(j)).Strike.Left(i) = [];
%                 end
%                 % loop through the stride and add up the sum of pressures per unit time
%                 for stride = DynamicPPTrials(Selection(j)).Strike.Left(i):DynamicPPTrials(Selection(j)).Off.Left(i)
%                     ind = stride - DynamicPPTrials(Selection(j)).Strike.Left(i) +1;
%                     % Heel
%                     DynamicPPTrials(Selection(j)).Area.Capture.L.LatHeel{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.L(i).LatHeel .* DynamicPPTrials(Selection(j)).TMLog(:,:,stride))) .* Adj.Area;  % number of sensors
%                     DynamicPPTrials(Selection(j)).Pressure.Capture.L.LatHeel{i}(ind) =  sum(sum(DynamicPPTrials(Selection(j)).Mask.L(i).LatHeel .* DynamicPPTrials(Selection(j)).TM(:,:,stride))) ./ 1000; % pressure in kPa
%                     DynamicPPTrials(Selection(j)).Force.Capture.L.LatHeel{i}(ind) = DynamicPPTrials(Selection(j)).Pressure.Capture.L.LatHeel{i}(ind) .* DynamicPPTrials(Selection(j)).Area.Capture.L.LatHeel{i}(ind); % force in N
%                     DynamicPPTrials(Selection(j)).Area.Capture.L.MedHeel{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.L(i).MedHeel .* DynamicPPTrials(Selection(j)).TMLog(:,:,stride))) .* Adj.Area;  % number of sensors
%                     DynamicPPTrials(Selection(j)).Pressure.Capture.L.MedHeel{i}(ind) =  sum(sum(DynamicPPTrials(Selection(j)).Mask.L(i).MedHeel .* DynamicPPTrials(Selection(j)).TM(:,:,stride))) ./ 1000; % pressure in Pa
%                     DynamicPPTrials(Selection(j)).Force.Capture.L.MedHeel{i}(ind) = DynamicPPTrials(Selection(j)).Pressure.Capture.L.MedHeel{i}(ind) .* DynamicPPTrials(Selection(j)).Area.Capture.L.MedHeel{i}(ind); % force in N
%                     % Arch
%                     DynamicPPTrials(Selection(j)).Area.Capture.L.LatArch{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.L(i).LatArch .* DynamicPPTrials(Selection(j)).TMLog(:,:,stride))) .* Adj.Area;  % sensor area adjusted to cm^2
%                     DynamicPPTrials(Selection(j)).Pressure.Capture.L.LatArch{i}(ind) =  sum(sum(DynamicPPTrials(Selection(j)).Mask.L(i).LatArch .* DynamicPPTrials(Selection(j)).TM(:,:,stride))) ./ 1000; % pressure in kPa
%                     DynamicPPTrials(Selection(j)).Force.Capture.L.LatArch{i}(ind) = DynamicPPTrials(Selection(j)).Pressure.Capture.L.LatArch{i}(ind) .* DynamicPPTrials(Selection(j)).Area.Capture.L.LatArch{i}(ind); % force in N
%                     DynamicPPTrials(Selection(j)).Area.Capture.L.MedArch{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.L(i).MedArch .* DynamicPPTrials(Selection(j)).TMLog(:,:,stride))) .* Adj.Area;  % area adjusted to cm^2
%                     DynamicPPTrials(Selection(j)).Pressure.Capture.L.MedArch{i}(ind) =  sum(sum(DynamicPPTrials(Selection(j)).Mask.L(i).MedArch .* DynamicPPTrials(Selection(j)).TM(:,:,stride))) ./ 1000; % pressure in Pa
%                     DynamicPPTrials(Selection(j)).Force.Capture.L.MedArch{i}(ind) = DynamicPPTrials(Selection(j)).Pressure.Capture.L.MedArch{i}(ind) .* DynamicPPTrials(Selection(j)).Area.Capture.L.MedArch{i}(ind); % force in N
%                     % Forefoot
%                     DynamicPPTrials(Selection(j)).Area.Capture.L.LatFore{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.L(i).LatFore .* DynamicPPTrials(Selection(j)).TMLog(:,:,stride))) .* Adj.Area;  % sensor area adjusted to cm^2
%                     DynamicPPTrials(Selection(j)).Pressure.Capture.L.LatFore{i}(ind) =  sum(sum(DynamicPPTrials(Selection(j)).Mask.L(i).LatFore .* DynamicPPTrials(Selection(j)).TM(:,:,stride))) ./ 1000; % pressure in kPa
%                     DynamicPPTrials(Selection(j)).Force.Capture.L.LatFore{i}(ind) = DynamicPPTrials(Selection(j)).Pressure.Capture.L.LatFore{i}(ind) .* DynamicPPTrials(Selection(j)).Area.Capture.L.LatFore{i}(ind); % force in N
%                     DynamicPPTrials(Selection(j)).Area.Capture.L.MedFore{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.L(i).MedFore .* DynamicPPTrials(Selection(j)).TMLog(:,:,stride))) .* Adj.Area;  % area adjusted to cm^2
%                     DynamicPPTrials(Selection(j)).Pressure.Capture.L.MedFore{i}(ind) =  sum(sum(DynamicPPTrials(Selection(j)).Mask.L(i).MedFore .* DynamicPPTrials(Selection(j)).TM(:,:,stride))) ./ 1000; % pressure in Pa
%                     DynamicPPTrials(Selection(j)).Force.Capture.L.MedFore{i}(ind) = DynamicPPTrials(Selection(j)).Pressure.Capture.L.MedFore{i}(ind) .* DynamicPPTrials(Selection(j)).Area.Capture.L.MedFore{i}(ind); % force in N
%                 end
%             end
%         end
%         % RIGHT Foot
%         % stop loop from running through more times than there are foot offs. If there isn't a clean foot off, its not a valid full step
%         if isfield(DynamicPPTrials(Selection(j)).Off, 'Right')
%             if length(DynamicPPTrials(Selection(j)).Off.Right) < length(DynamicPPTrials(Selection(j)).Strike.Right)
%                 DynamicPPTrials(Selection(j)).NumFullStepsR = length(DynamicPPTrials(Selection(j)).Off.Right);
%             else
%                 DynamicPPTrials(Selection(j)).NumFullStepsR = DynamicPPTrials(Selection(j)).NumRight;
%             end
%         end
%         if isfield(DynamicPPTrials(Selection(j)), 'NumFullStepsR')
%             for i = 1:DynamicPPTrials(Selection(j)).NumFullStepsR
%                 if DynamicPPTrials(Selection(j)).Strike.Right(i) == 0 && DynamicPPTrials(Selection(j)).Off.Right(i)
%                     DynamicPPTrials(Selection(j)).Off.Right(i) = [];
%                     DynamicPPTrials(Selection(j)).Strike.Right(i) = [];
%                 end
%                 for stride = DynamicPPTrials(Selection(j)).Strike.Right(i):DynamicPPTrials(Selection(j)).Off.Right(i)
%                     ind = stride - DynamicPPTrials(Selection(j)).Strike.Right(i) +1; % sample counter
%                     % Heel
%                     DynamicPPTrials(Selection(j)).Area.Capture.R.LatHeel{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.R(i).LatHeel .* DynamicPPTrials(Selection(j)).TMLog(:,:,stride))) .* Adj.Area;  % sensor area adjusted to cm^2
%                     DynamicPPTrials(Selection(j)).Pressure.Capture.R.LatHeel{i}(ind) =  sum(sum(DynamicPPTrials(Selection(j)).Mask.R(i).LatHeel .* DynamicPPTrials(Selection(j)).TM(:,:,stride))) ./ 1000; % pressure in kPa
%                     DynamicPPTrials(Selection(j)).Force.Capture.R.LatHeel{i}(ind) = DynamicPPTrials(Selection(j)).Pressure.Capture.R.LatHeel{i}(ind) .* DynamicPPTrials(Selection(j)).Area.Capture.R.LatHeel{i}(ind); % force in N
%                     DynamicPPTrials(Selection(j)).Area.Capture.R.MedHeel{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.R(i).MedHeel .* DynamicPPTrials(Selection(j)).TMLog(:,:,stride))) .* Adj.Area;  % area adjusted to cm^2
%                     DynamicPPTrials(Selection(j)).Pressure.Capture.R.MedHeel{i}(ind) =  sum(sum(DynamicPPTrials(Selection(j)).Mask.R(i).MedHeel .* DynamicPPTrials(Selection(j)).TM(:,:,stride))) ./ 1000; % pressure in Pa
%                     DynamicPPTrials(Selection(j)).Force.Capture.R.MedHeel{i}(ind) = DynamicPPTrials(Selection(j)).Pressure.Capture.R.MedHeel{i}(ind) .* DynamicPPTrials(Selection(j)).Area.Capture.R.MedHeel{i}(ind); % force in N
%                     % Arch
%                     DynamicPPTrials(Selection(j)).Area.Capture.R.LatArch{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.R(i).LatArch .* DynamicPPTrials(Selection(j)).TMLog(:,:,stride))) .* Adj.Area;  % sensor area adjusted to cm^2
%                     DynamicPPTrials(Selection(j)).Pressure.Capture.R.LatArch{i}(ind) =  sum(sum(DynamicPPTrials(Selection(j)).Mask.R(i).LatArch .* DynamicPPTrials(Selection(j)).TM(:,:,stride))) ./ 1000; % pressure in kPa
%                     DynamicPPTrials(Selection(j)).Force.Capture.R.LatArch{i}(ind) = DynamicPPTrials(Selection(j)).Pressure.Capture.R.LatArch{i}(ind) .* DynamicPPTrials(Selection(j)).Area.Capture.R.LatArch{i}(ind); % force in N
%                     DynamicPPTrials(Selection(j)).Area.Capture.R.MedArch{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.R(i).MedArch .* DynamicPPTrials(Selection(j)).TMLog(:,:,stride))) .* Adj.Area;  % area adjusted to cm^2
%                     DynamicPPTrials(Selection(j)).Pressure.Capture.R.MedArch{i}(ind) =  sum(sum(DynamicPPTrials(Selection(j)).Mask.R(i).MedArch .* DynamicPPTrials(Selection(j)).TM(:,:,stride))) ./ 1000; % pressure in Pa
%                     DynamicPPTrials(Selection(j)).Force.Capture.R.MedArch{i}(ind) = DynamicPPTrials(Selection(j)).Pressure.Capture.R.MedArch{i}(ind) .* DynamicPPTrials(Selection(j)).Area.Capture.R.MedArch{i}(ind); % force in N
%                     % Forefoot
%                     DynamicPPTrials(Selection(j)).Area.Capture.R.LatFore{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.R(i).LatFore .* DynamicPPTrials(Selection(j)).TMLog(:,:,stride))) .* Adj.Area;  % sensor area adjusted to cm^2
%                     DynamicPPTrials(Selection(j)).Pressure.Capture.R.LatFore{i}(ind) =  sum(sum(DynamicPPTrials(Selection(j)).Mask.R(i).LatFore .* DynamicPPTrials(Selection(j)).TM(:,:,stride))) ./ 1000; % pressure in kPa
%                     DynamicPPTrials(Selection(j)).Force.Capture.R.LatFore{i}(ind) = DynamicPPTrials(Selection(j)).Pressure.Capture.R.LatFore{i}(ind) .* DynamicPPTrials(Selection(j)).Area.Capture.R.LatFore{i}(ind); % force in N
%                     DynamicPPTrials(Selection(j)).Area.Capture.R.MedFore{i}(ind) = sum(sum(DynamicPPTrials(Selection(j)).Mask.R(i).MedFore .* DynamicPPTrials(Selection(j)).TMLog(:,:,stride))) .* Adj.Area;  % area adjusted to cm^2
%                     DynamicPPTrials(Selection(j)).Pressure.Capture.R.MedFore{i}(ind) =  sum(sum(DynamicPPTrials(Selection(j)).Mask.R(i).MedFore .* DynamicPPTrials(Selection(j)).TM(:,:,stride))) ./ 1000; % pressure in Pa
%                     DynamicPPTrials(Selection(j)).Force.Capture.R.MedFore{i}(ind) = DynamicPPTrials(Selection(j)).Pressure.Capture.R.MedFore{i}(ind) .* DynamicPPTrials(Selection(j)).Area.Capture.R.MedFore{i}(ind); % force in N
%                 end
%             end
%         end
%     end
% end
% 
% %% Interpolate stance duration to % of stance period
% % LEFT
% for j = 1:length(Selection)
%     for i = 1:DynamicPPTrials(Selection(j)).NumFullStepsL
%         a = length(DynamicPPTrials(Selection(j)).Pressure.Capture.L.LatFore{i});
%         % Forefoot
%         DynamicPPTrials(Selection(j)).Pressure.ReSam.L.LatFore{i} = resample(DynamicPPTrials(Selection(j)).Pressure.Capture.L.LatFore{i}(:),100, a);
%         DynamicPPTrials(Selection(j)).Pressure.ReSam.L.MedFore{i} = resample(DynamicPPTrials(Selection(j)).Pressure.Capture.L.MedFore{i}(:),100, a);
%         DynamicPPTrials(Selection(j)).Force.ReSam.L.LatFore{i} = resample(DynamicPPTrials(Selection(j)).Force.Capture.L.LatFore{i}(:),100, a);
%         DynamicPPTrials(Selection(j)).Force.ReSam.L.MedFore{i} = resample(DynamicPPTrials(Selection(j)).Force.Capture.L.MedFore{i}(:),100, a);
%         DynamicPPTrials(Selection(j)).Area.ReSam.L.LatFore{i} = resample(DynamicPPTrials(Selection(j)).Area.Capture.L.LatFore{i}(:),100, a);
%         DynamicPPTrials(Selection(j)).Area.ReSam.L.MedFore{i} = resample(DynamicPPTrials(Selection(j)).Area.Capture.L.MedFore{i}(:),100, a);
%         % Arch
%         DynamicPPTrials(Selection(j)).Pressure.ReSam.L.LatArch{i} = resample(DynamicPPTrials(Selection(j)).Pressure.Capture.L.LatArch{i}(:),100, a);
%         DynamicPPTrials(Selection(j)).Pressure.ReSam.L.MedArch{i} = resample(DynamicPPTrials(Selection(j)).Pressure.Capture.L.MedArch{i}(:),100, a);
%         DynamicPPTrials(Selection(j)).Force.ReSam.L.LatArch{i} = resample(DynamicPPTrials(Selection(j)).Force.Capture.L.LatArch{i}(:),100, a);
%         DynamicPPTrials(Selection(j)).Force.ReSam.L.MedArch{i} = resample(DynamicPPTrials(Selection(j)).Force.Capture.L.MedArch{i}(:),100, a);
%         DynamicPPTrials(Selection(j)).Area.ReSam.L.LatArch{i} = resample(DynamicPPTrials(Selection(j)).Area.Capture.L.LatArch{i}(:),100, a);
%         DynamicPPTrials(Selection(j)).Area.ReSam.L.MedArch{i} = resample(DynamicPPTrials(Selection(j)).Area.Capture.L.MedArch{i}(:),100, a);
%         % Heel
%         DynamicPPTrials(Selection(j)).Pressure.ReSam.L.LatHeel{i} = resample(DynamicPPTrials(Selection(j)).Pressure.Capture.L.LatHeel{i}(:),100, a);
%         DynamicPPTrials(Selection(j)).Pressure.ReSam.L.MedHeel{i} = resample(DynamicPPTrials(Selection(j)).Pressure.Capture.L.MedHeel{i}(:),100, a);
%         DynamicPPTrials(Selection(j)).Force.ReSam.L.LatHeel{i} = resample(DynamicPPTrials(Selection(j)).Force.Capture.L.LatHeel{i}(:),100, a);
%         DynamicPPTrials(Selection(j)).Force.ReSam.L.MedHeel{i} = resample(DynamicPPTrials(Selection(j)).Force.Capture.L.MedHeel{i}(:),100, a);
%         DynamicPPTrials(Selection(j)).Area.ReSam.L.LatHeel{i} = resample(DynamicPPTrials(Selection(j)).Area.Capture.L.LatHeel{i}(:),100, a);
%         DynamicPPTrials(Selection(j)).Area.ReSam.L.MedHeel{i} = resample(DynamicPPTrials(Selection(j)).Area.Capture.L.MedHeel{i}(:),100, a);
%         % Whole Foot
%         DynamicPPTrials(Selection(j)).Pressure.ReSam.L.WholeFoot{i} = DynamicPPTrials(Selection(j)).Pressure.ReSam.L.LatHeel{i} + DynamicPPTrials(Selection(j)).Pressure.ReSam.L.MedHeel{i} + ...
%             DynamicPPTrials(Selection(j)).Pressure.ReSam.L.LatArch{i} + DynamicPPTrials(Selection(j)).Pressure.ReSam.L.MedHeel{i} + ...
%             DynamicPPTrials(Selection(j)).Pressure.ReSam.L.LatFore{i} + DynamicPPTrials(Selection(j)).Pressure.ReSam.L.MedFore{i};
%         DynamicPPTrials(Selection(j)).Area.ReSam.L.WholeFoot{i} = DynamicPPTrials(Selection(j)).Area.ReSam.L.LatHeel{i} + DynamicPPTrials(Selection(j)).Area.ReSam.L.MedHeel{i} + ...
%             DynamicPPTrials(Selection(j)).Area.ReSam.L.LatArch{i} + DynamicPPTrials(Selection(j)).Area.ReSam.L.MedHeel{i} + ...
%             DynamicPPTrials(Selection(j)).Area.ReSam.L.LatFore{i} + DynamicPPTrials(Selection(j)).Area.ReSam.L.MedFore{i};
%         DynamicPPTrials(Selection(j)).Force.ReSam.L.WholeFoot{i} = DynamicPPTrials(Selection(j)).Force.ReSam.L.LatHeel{i} + DynamicPPTrials(Selection(j)).Force.ReSam.L.MedHeel{i} + ...
%             DynamicPPTrials(Selection(j)).Force.ReSam.L.LatArch{i} + DynamicPPTrials(Selection(j)).Force.ReSam.L.MedArch{i} + ...
%             DynamicPPTrials(Selection(j)).Force.ReSam.L.LatFore{i} + DynamicPPTrials(Selection(j)).Force.ReSam.L.MedFore{i};
%     end
%     % RIGHT
%     for i = 1:DynamicPPTrials(Selection(j)).NumFullStepsR
%         a = length(DynamicPPTrials(Selection(j)).Pressure.Capture.R.LatFore{i});
%         % Forefoot
%         DynamicPPTrials(Selection(j)).Pressure.ReSam.R.LatFore{i} = resample(DynamicPPTrials(Selection(j)).Pressure.Capture.R.LatFore{i}(:),100, a);
%         DynamicPPTrials(Selection(j)).Pressure.ReSam.R.MedFore{i} = resample(DynamicPPTrials(Selection(j)).Pressure.Capture.R.MedFore{i}(:),100, a);
%         DynamicPPTrials(Selection(j)).Force.ReSam.R.LatFore{i} = resample(DynamicPPTrials(Selection(j)).Force.Capture.R.LatFore{i}(:),100, a);
%         DynamicPPTrials(Selection(j)).Force.ReSam.R.MedFore{i} = resample(DynamicPPTrials(Selection(j)).Force.Capture.R.MedFore{i}(:),100, a);
%         DynamicPPTrials(Selection(j)).Area.ReSam.R.LatFore{i} = resample(DynamicPPTrials(Selection(j)).Area.Capture.R.LatFore{i}(:),100, a);
%         DynamicPPTrials(Selection(j)).Area.ReSam.R.MedFore{i} = resample(DynamicPPTrials(Selection(j)).Area.Capture.R.MedFore{i}(:),100, a);
%         % Arch
%         DynamicPPTrials(Selection(j)).Pressure.ReSam.R.LatArch{i} = resample(DynamicPPTrials(Selection(j)).Pressure.Capture.R.LatArch{i}(:),100, a);
%         DynamicPPTrials(Selection(j)).Pressure.ReSam.R.MedArch{i} = resample(DynamicPPTrials(Selection(j)).Pressure.Capture.R.MedArch{i}(:),100, a);
%         DynamicPPTrials(Selection(j)).Force.ReSam.R.LatArch{i} = resample(DynamicPPTrials(Selection(j)).Force.Capture.R.LatArch{i}(:),100, a);
%         DynamicPPTrials(Selection(j)).Force.ReSam.R.MedArch{i} = resample(DynamicPPTrials(Selection(j)).Force.Capture.R.MedArch{i}(:),100, a);
%         DynamicPPTrials(Selection(j)).Area.ReSam.R.LatArch{i} = resample(DynamicPPTrials(Selection(j)).Area.Capture.R.LatArch{i}(:),100, a);
%         DynamicPPTrials(Selection(j)).Area.ReSam.R.MedArch{i} = resample(DynamicPPTrials(Selection(j)).Area.Capture.R.MedArch{i}(:),100, a);
%         % Heel
%         DynamicPPTrials(Selection(j)).Pressure.ReSam.R.LatHeel{i} = resample(DynamicPPTrials(Selection(j)).Pressure.Capture.R.LatHeel{i}(:),100, a);
%         DynamicPPTrials(Selection(j)).Pressure.ReSam.R.MedHeel{i} = resample(DynamicPPTrials(Selection(j)).Pressure.Capture.R.MedHeel{i}(:),100, a);
%         DynamicPPTrials(Selection(j)).Force.ReSam.R.LatHeel{i} = resample(DynamicPPTrials(Selection(j)).Force.Capture.R.LatHeel{i}(:),100, a);
%         DynamicPPTrials(Selection(j)).Force.ReSam.R.MedHeel{i} = resample(DynamicPPTrials(Selection(j)).Force.Capture.R.MedHeel{i}(:),100, a);
%         DynamicPPTrials(Selection(j)).Area.ReSam.R.LatHeel{i} = resample(DynamicPPTrials(Selection(j)).Area.Capture.R.LatHeel{i}(:),100, a);
%         DynamicPPTrials(Selection(j)).Area.ReSam.R.MedHeel{i} = resample(DynamicPPTrials(Selection(j)).Area.Capture.R.MedHeel{i}(:),100, a);
%         % Whole Foot
%         DynamicPPTrials(Selection(j)).Pressure.ReSam.R.WholeFoot{i} = DynamicPPTrials(Selection(j)).Pressure.ReSam.R.LatHeel{i} + DynamicPPTrials(Selection(j)).Pressure.ReSam.R.MedHeel{i} + ...
%             DynamicPPTrials(Selection(j)).Pressure.ReSam.R.LatArch{i} + DynamicPPTrials(Selection(j)).Pressure.ReSam.R.MedHeel{i} + ...
%             DynamicPPTrials(Selection(j)).Pressure.ReSam.R.LatFore{i} + DynamicPPTrials(Selection(j)).Pressure.ReSam.R.MedFore{i};
%         DynamicPPTrials(Selection(j)).Area.ReSam.R.WholeFoot{i} = DynamicPPTrials(Selection(j)).Area.ReSam.R.LatHeel{i} + DynamicPPTrials(Selection(j)).Area.ReSam.R.MedHeel{i} + ...
%             DynamicPPTrials(Selection(j)).Area.ReSam.R.LatArch{i} + DynamicPPTrials(Selection(j)).Area.ReSam.R.MedHeel{i} + ...
%             DynamicPPTrials(Selection(j)).Area.ReSam.R.LatFore{i} + DynamicPPTrials(Selection(j)).Area.ReSam.R.MedFore{i};
%         DynamicPPTrials(Selection(j)).Force.ReSam.R.WholeFoot{i} = DynamicPPTrials(Selection(j)).Force.ReSam.R.LatHeel{i} + DynamicPPTrials(Selection(j)).Force.ReSam.R.MedHeel{i} + ...
%             DynamicPPTrials(Selection(j)).Force.ReSam.R.LatArch{i} + DynamicPPTrials(Selection(j)).Force.ReSam.R.MedArch{i} + ...
%             DynamicPPTrials(Selection(j)).Force.ReSam.R.LatFore{i} + DynamicPPTrials(Selection(j)).Force.ReSam.R.MedFore{i};
%     end
% end
% 
% %% Normalize pressures and forces to % weight
% % normalize to shoe area
% % LEFT
% for j = 1:length(Selection)
%     for i = 1:DynamicPPTrials(Selection(j)).NumFullStepsL
%         % Forefoot
%         DynamicPPTrials(Selection(j)).Pressure.Norm.L.LatFore{i} = DynamicPPTrials(Selection(j)).Pressure.ReSam.L.LatFore{i} ./ Subject.Newtons;
%         DynamicPPTrials(Selection(j)).Pressure.Norm.L.MedFore{i} =  DynamicPPTrials(Selection(j)).Pressure.ReSam.L.MedFore{i} ./ Subject.Newtons;
%         DynamicPPTrials(Selection(j)).Force.Norm.L.LatFore{i} = DynamicPPTrials(Selection(j)).Force.ReSam.L.LatFore{i} ./ Subject.Newtons;
%         DynamicPPTrials(Selection(j)).Force.Norm.L.MedFore{i} = DynamicPPTrials(Selection(j)).Force.ReSam.L.MedFore{i} ./ Subject.Newtons;
%         % Arch
%         DynamicPPTrials(Selection(j)).Pressure.Norm.L.LatArch{i} = DynamicPPTrials(Selection(j)).Pressure.ReSam.L.LatArch{i} ./ Subject.Newtons;
%         DynamicPPTrials(Selection(j)).Pressure.Norm.L.MedArch{i} = DynamicPPTrials(Selection(j)).Pressure.ReSam.L.MedArch{i} ./ Subject.Newtons;
%         DynamicPPTrials(Selection(j)).Force.Norm.L.LatArch{i} = DynamicPPTrials(Selection(j)).Force.ReSam.L.LatArch{i} ./ Subject.Newtons;
%         DynamicPPTrials(Selection(j)).Force.Norm.L.MedArch{i} = DynamicPPTrials(Selection(j)).Force.ReSam.L.MedArch{i} ./ Subject.Newtons;
%         % Heel
%         DynamicPPTrials(Selection(j)).Pressure.Norm.L.LatHeel{i} = DynamicPPTrials(Selection(j)).Pressure.ReSam.L.LatHeel{i} ./ Subject.Newtons;
%         DynamicPPTrials(Selection(j)).Pressure.Norm.L.MedHeel{i} = DynamicPPTrials(Selection(j)).Pressure.ReSam.L.MedHeel{i} ./ Subject.Newtons;
%         DynamicPPTrials(Selection(j)).Force.Norm.L.LatHeel{i} = DynamicPPTrials(Selection(j)).Force.ReSam.L.LatHeel{i} ./ Subject.Newtons;
%         DynamicPPTrials(Selection(j)).Force.Norm.L.MedHeel{i} = DynamicPPTrials(Selection(j)).Force.ReSam.L.MedHeel{i} ./ Subject.Newtons;
%         % Whole Foot
%         DynamicPPTrials(Selection(j)).Pressure.Norm.L.WholeFoot{i} = DynamicPPTrials(Selection(j)).Pressure.Norm.L.LatHeel{i} + DynamicPPTrials(Selection(j)).Pressure.Norm.L.MedHeel{i} + ...
%             DynamicPPTrials(Selection(j)).Pressure.Norm.L.LatArch{i} + DynamicPPTrials(Selection(j)).Pressure.Norm.L.MedHeel{i} + ...
%             DynamicPPTrials(Selection(j)).Pressure.Norm.L.LatFore{i} + DynamicPPTrials(Selection(j)).Pressure.Norm.L.MedFore{i};
%         DynamicPPTrials(Selection(j)).Force.Norm.L.WholeFoot{i} = DynamicPPTrials(Selection(j)).Force.Norm.L.LatHeel{i} + DynamicPPTrials(Selection(j)).Force.Norm.L.MedHeel{i} + ...
%             DynamicPPTrials(Selection(j)).Force.Norm.L.LatArch{i} + DynamicPPTrials(Selection(j)).Force.Norm.L.MedArch{i} + ...
%             DynamicPPTrials(Selection(j)).Force.Norm.L.LatFore{i} + DynamicPPTrials(Selection(j)).Force.Norm.L.MedFore{i};
%     end
%     % RIGHT
%     for i = 1:DynamicPPTrials(Selection(j)).NumFullStepsR
%         % Forefoot
%         DynamicPPTrials(Selection(j)).Pressure.Norm.R.LatFore{i} = DynamicPPTrials(Selection(j)).Pressure.ReSam.R.LatFore{i} ./ Subject.Newtons;
%         DynamicPPTrials(Selection(j)).Pressure.Norm.R.MedFore{i} = DynamicPPTrials(Selection(j)).Pressure.ReSam.R.MedFore{i} ./ Subject.Newtons;
%         DynamicPPTrials(Selection(j)).Force.Norm.R.LatFore{i} = DynamicPPTrials(Selection(j)).Force.ReSam.R.LatFore{i} ./ Subject.Newtons;
%         DynamicPPTrials(Selection(j)).Force.Norm.R.MedFore{i} = DynamicPPTrials(Selection(j)).Force.ReSam.R.MedFore{i} ./ Subject.Newtons;
%         % Arch
%         DynamicPPTrials(Selection(j)).Pressure.Norm.R.LatArch{i} = DynamicPPTrials(Selection(j)).Pressure.ReSam.R.LatArch{i} ./ Subject.Newtons;
%         DynamicPPTrials(Selection(j)).Pressure.Norm.R.MedArch{i} = DynamicPPTrials(Selection(j)).Pressure.ReSam.R.MedArch{i} ./ Subject.Newtons;
%         DynamicPPTrials(Selection(j)).Force.Norm.R.LatArch{i} = DynamicPPTrials(Selection(j)).Force.ReSam.R.LatArch{i} ./ Subject.Newtons;
%         DynamicPPTrials(Selection(j)).Force.Norm.R.MedArch{i} = DynamicPPTrials(Selection(j)).Force.ReSam.R.MedArch{i} ./ Subject.Newtons;
%         % Heel
%         DynamicPPTrials(Selection(j)).Pressure.Norm.R.LatHeel{i} = DynamicPPTrials(Selection(j)).Pressure.ReSam.R.LatHeel{i} ./ Subject.Newtons;
%         DynamicPPTrials(Selection(j)).Pressure.Norm.R.MedHeel{i} = DynamicPPTrials(Selection(j)).Pressure.ReSam.R.MedHeel{i} ./ Subject.Newtons;
%         DynamicPPTrials(Selection(j)).Force.Norm.R.LatHeel{i} = DynamicPPTrials(Selection(j)).Force.ReSam.R.LatHeel{i} ./ Subject.Newtons;
%         DynamicPPTrials(Selection(j)).Force.Norm.R.MedHeel{i} = DynamicPPTrials(Selection(j)).Force.ReSam.R.MedHeel{i} ./ Subject.Newtons;
%         % Whole Foot
%         DynamicPPTrials(Selection(j)).Pressure.Norm.R.WholeFoot{i} = DynamicPPTrials(Selection(j)).Pressure.Norm.R.LatHeel{i} + DynamicPPTrials(Selection(j)).Pressure.Norm.R.MedHeel{i} + ...
%             DynamicPPTrials(Selection(j)).Pressure.Norm.R.LatArch{i} + DynamicPPTrials(Selection(j)).Pressure.Norm.R.MedHeel{i} + ...
%             DynamicPPTrials(Selection(j)).Pressure.Norm.R.LatFore{i} + DynamicPPTrials(Selection(j)).Pressure.Norm.R.MedFore{i};
%         DynamicPPTrials(Selection(j)).Force.Norm.R.WholeFoot{i} = DynamicPPTrials(Selection(j)).Force.Norm.R.LatHeel{i} + DynamicPPTrials(Selection(j)).Force.Norm.R.MedHeel{i} + ...
%             DynamicPPTrials(Selection(j)).Force.Norm.R.LatArch{i} + DynamicPPTrials(Selection(j)).Force.Norm.R.MedArch{i} + ...
%             DynamicPPTrials(Selection(j)).Force.Norm.R.LatFore{i} + DynamicPPTrials(Selection(j)).Force.Norm.R.MedFore{i};
%     end
% end
% 
% 
% %% Create forces figure
% if strcmp(PPPlots.AllForces, 'Yes')
%     % create pressure figure
%     TimeSeriesForce = figure('Position', [50, 50, 900, 700]);
%     x = 1:100;
%     % LEFT Foot
%     subplot(221);
%     hold on;
%     for j = 1:length(Selection)
%         for k = 1:DynamicPPTrials(Selection(j)).NumFullStepsL
%             LF = plot(x,DynamicPPTrials(Selection(j)).Force.Norm.L.LatFore{k},'b--');
%             plot(x,DynamicPPTrials(Selection(j)).Force.Norm.L.MedFore{k},'b.');
%             LA = plot(x,DynamicPPTrials(Selection(j)).Force.Norm.L.LatArch{k},'g--');
%             plot(x,DynamicPPTrials(Selection(j)).Force.Norm.L.MedArch{k},'g.');
%             LH = plot(x,DynamicPPTrials(Selection(j)).Force.Norm.L.LatHeel{k},'r--');
%             plot(x,DynamicPPTrials(Selection(j)).Force.Norm.L.MedHeel{k},'r.');
%         end
%     end
%     title('Left Foot Forces Med/Lat');
%     ylabel('% of Weight');
%     ylim([0 1.5]);
%     ax = gca;
%     ax.YTick = ([0 .25 .50 .75 1 1.25 1.5]);
%     ax.YTickLabel = ({'0%','25%','50%','75%','100%','125%','150%'});
%     legend([LF LA LH],{'Lat Fore' 'Lat Mid' 'Lat Hind'}, 'Location','North');
%     ax.XTick = ([0 25 50 75 100]);
%     ax.XTickLabel = [];
%     grid on;
%     
%     subplot(223);
%     hold on;
%     % plot general heel, arch, fore, and whole foot pressures
%     for j = 1:length(Selection)
%         for k = 1:DynamicPPTrials(Selection(j)).NumFullStepsL
%             Lfore = (DynamicPPTrials(Selection(j)).Force.Norm.L.LatFore{k} + DynamicPPTrials(Selection(j)).Force.Norm.L.MedFore{k});
%             Larch = (DynamicPPTrials(Selection(j)).Force.Norm.L.LatArch{k} + DynamicPPTrials(Selection(j)).Force.Norm.L.MedArch{k});
%             Lheel = (DynamicPPTrials(Selection(j)).Force.Norm.L.LatHeel{k} + DynamicPPTrials(Selection(j)).Force.Norm.L.MedHeel{k});
%             plot(x, Lfore, 'b-', 'LineWidth', 1);
%             plot(x, Larch, 'g-', 'LineWidth', 1);
%             plot(x, Lheel, 'r-', 'LineWidth', 1);
%             plot(x, Lfore + Larch + Lheel, 'k-','LineWidth',2);
%         end
%     end
%     % create bar series along the top to display timing norms
%     barh(2.8,100,'b');
%     barh(2.8,53,'g');
%     barh(2.8,24,'r');
%     title('Left Foot Forces');
%     ylabel('% of Weight');
%     ylim([0 2.5]);
%     xlabel('% of Stance Period');
%     xlim([0 100]);
%     ax = gca;
%     %ax.XTickLabel = [];
%     ax.YTick = ([0 .25 .50 .75 1 1.25 1.5 1.75 2 2.25 2.5]);
%     ax.YTickLabel = ({'0%','','50%','','100%','','150%','','200%','','250%'});
%     ax.XTick = ([0 25 50 75 100]);
%     ax.XTickLabel = ({'0%','25%','50%','75%','100%'});
%     grid on;
%     
%     subplot(222);
%     hold on;
%     for j = 1:length(Selection)
%         for k = 1:DynamicPPTrials(Selection(j)).NumFullStepsR
%             plot(x,DynamicPPTrials(Selection(j)).Force.Norm.R.LatFore{k},'b--');
%             MF = plot(x,DynamicPPTrials(Selection(j)).Force.Norm.R.MedFore{k},'b.');
%             plot(x,DynamicPPTrials(Selection(j)).Force.Norm.R.LatArch{k},'g--');
%             MA = plot(x,DynamicPPTrials(Selection(j)).Force.Norm.R.MedArch{k},'g.');
%             plot(x,DynamicPPTrials(Selection(j)).Force.Norm.R.LatHeel{k},'r--');
%             MH = plot(x,DynamicPPTrials(Selection(j)).Force.Norm.R.MedHeel{k},'r.');
%         end
%     end
%     title('Right Foot Forces Med/Lat');
%     ylim([0 1.5]);
%     ax = gca;
%     ax.YTick = ([0 .25 .50 .75 1 1.25 1.5]);
%     ax.YTickLabel = [];
%     ax.XTick = ([0 25 50 75 100]);
%     ax.XTickLabel = [];
%     grid on;
%     
%     subplot(224);
%     hold on;
%     for j = 1:length(Selection)
%         for k = 1:DynamicPPTrials(Selection(j)).NumFullStepsR
%             % plot general heel, arch, fore, and whole foot pressures
%             Rfore = (DynamicPPTrials(Selection(j)).Force.Norm.R.LatFore{k} + DynamicPPTrials(Selection(j)).Force.Norm.R.MedFore{k});
%             Rarch = (DynamicPPTrials(Selection(j)).Force.Norm.R.LatArch{k} + DynamicPPTrials(Selection(j)).Force.Norm.R.MedArch{k});
%             Rheel = (DynamicPPTrials(Selection(j)).Force.Norm.R.LatHeel{k} + DynamicPPTrials(Selection(j)).Force.Norm.R.MedHeel{k});
%             plot(x, Rfore, 'b-', 'LineWidth', 1);
%             plot(x, Rarch, 'g-', 'LineWidth', 1);
%             plot(x, Rheel, 'r-', 'LineWidth', 1);
%             plot(x, Rfore + Rarch + Rheel, 'k-','LineWidth',2);
%         end
%     end
%     % create bar series along the top to display timing norms
%     barh(2.8,100,'b');
%     barh(2.8,53,'g');
%     barh(2.8,24,'r');
%     
%     title('Right Foot Forces');
%     xlabel('% of Stance Period');
%     xlim([0 100]);
%     ylim([0 2.5]);
%     legend([MF MA MH],{'Med Fore' 'Med Mid' 'Med Hind'}, 'Location','North');
%     ax = gca;
%     ax.XTick = ([0 25 50 75 100]);
%     ax.XTickLabel = ({'0%','25%','50%','75%','100%'});
%     ax.YTick = ([0 .25 .50 .75 1 1.25 1.5 1.75 2 2.25 2.5]);
%     ax.YTickLabel = [];
%     grid on;
%     
%     subplotsqueeze(TimeSeriesForce, 1.2);
%     % save and export figure to excel
%     saveas(TimeSeriesForce,strcat(folder,'\','Force Distribution.png'));
%     if strcmp(PPSettings.ExportReport, 'Yes') == 1
%         xlsPasteTo('ImpressionsOutput.xlsx','Sheet1',610,600, 'AE17');
%     end
%     clearvars x LF MF LA MA LH MH Rfore Rarch Rheel Lfore Larch Lheel
% end
% 
% %% Create Pressure figure
% if strcmp(PPPlots.AllPressures, 'Yes')
%     TimeSeriesPress = figure('Position', [50, 50, 900, 700]);
%     x = 1:100;
%     % LEFT Foot
%     subplot(221);
%     hold on;
%     for j = 1:length(Selection)
%         for k = 1:DynamicPPTrials(Selection(j)).NumFullStepsL
%             LF = plot(x,DynamicPPTrials(Selection(j)).Pressure.Norm.L.LatFore{k},'b--');
%             plot(x,DynamicPPTrials(Selection(j)).Pressure.Norm.L.MedFore{k},'b.');
%             LA = plot(x,DynamicPPTrials(Selection(j)).Pressure.Norm.L.LatArch{k},'g--');
%             plot(x,DynamicPPTrials(Selection(j)).Pressure.Norm.L.MedArch{k},'g.');
%             LH = plot(x,DynamicPPTrials(Selection(j)).Pressure.Norm.L.LatHeel{k},'r--');
%             plot(x,DynamicPPTrials(Selection(j)).Pressure.Norm.L.MedHeel{k},'r.');
%         end
%     end
%     title('Left Foot Pressures Med/Lat');
%     ylabel('Normalized Pressure (N / cm^2 / kg)');
%     % ylim([0 1]);
%     ax = gca;
%     % ax.YTick = ([0 .25 .50 .75 1]);
%     % ax.YTickLabel = ({'0%','25%','50%','75%','100%'});
%     legend([LF LA LH],{'Lat Fore' 'Lat Arch' 'Lat Heel'}, 'Location','North');
%     ax.XTick = ([0 25 50 75 100]);
%     ax.XTickLabel = [];
%     grid on;
%     
%     subplot(223);
%     hold on;
%     % plot general heel, arch, fore, and whole foot pressures
%     for j = 1:length(Selection)
%         for k = 1:DynamicPPTrials(Selection(j)).NumFullStepsL
%             Lfore = (DynamicPPTrials(Selection(j)).Pressure.Norm.L.LatFore{k} + DynamicPPTrials(Selection(j)).Pressure.Norm.L.MedFore{k});
%             Larch = (DynamicPPTrials(Selection(j)).Pressure.Norm.L.LatArch{k} + DynamicPPTrials(Selection(j)).Pressure.Norm.L.MedArch{k});
%             Lheel = (DynamicPPTrials(Selection(j)).Pressure.Norm.L.LatHeel{k} + DynamicPPTrials(Selection(j)).Pressure.Norm.L.MedHeel{k});
%             plot(x, Lfore, 'b-', 'LineWidth', 1);
%             plot(x, Larch, 'g-', 'LineWidth', 1);
%             plot(x, Lheel, 'r-', 'LineWidth', 1);
%             plot(x, Lfore + Larch + Lheel, 'k-','LineWidth',2);
%         end
%     end
%     % create bar series along the top to display timing norms
%     % barh(2.3,100,'b');
%     % barh(2.3,53,'g');
%     % barh(2.3,24,'r');
%     title('Left Foot Pressures');
%     ylabel('Normalized Pressure (N / cm^2 / kg)');
%     % ylim([0 2]);
%     xlabel('% of Stance Period');
%     xlim([0 100]);
%     ax = gca;
%     %ax.XTickLabel = [];
%     % ax.YTick = ([0 .25 .50 .75 1 1.25 1.5 1.75 2]);
%     % ax.YTickLabel = ({'0%','25%','50%','75%','100%','125%','150%','175%','200%'});
%     ax.XTick = ([0 25 50 75 100]);
%     ax.XTickLabel = ({'0%','25%','50%','75%','100%'});
%     grid on;
%     
%     subplot(222);
%     hold on;
%     for j = 1:length(Selection)
%         for k = 1:DynamicPPTrials(Selection(j)).NumFullStepsR
%             plot(x,DynamicPPTrials(Selection(j)).Pressure.Norm.R.LatFore{k},'b--');
%             MF = plot(x,DynamicPPTrials(Selection(j)).Pressure.Norm.R.MedFore{k},'b.');
%             plot(x,DynamicPPTrials(Selection(j)).Pressure.Norm.R.LatArch{k},'g--');
%             MA = plot(x,DynamicPPTrials(Selection(j)).Pressure.Norm.R.MedArch{k},'g.');
%             plot(x,DynamicPPTrials(Selection(j)).Pressure.Norm.R.LatHeel{k},'r--');
%             MH = plot(x,DynamicPPTrials(Selection(j)).Pressure.Norm.R.MedHeel{k},'r.');
%         end
%     end
%     title('Right Foot Pressures Med/Lat');
%     % ylim([0 1]);
%     ax = gca;
%     % ax.YTick = ([0 .25 .50 .75 1]);
%     % ax.YTickLabel = [];
%     ax.XTick = ([0 25 50 75 100]);
%     ax.XTickLabel = [];
%     grid on;
%     
%     subplot(224);
%     hold on;
%     for j = 1:length(Selection)
%         for k = 1:DynamicPPTrials(Selection(j)).NumFullStepsR
%             % plot general heel, arch, fore, and whole foot pressures
%             Rfore = (DynamicPPTrials(Selection(j)).Pressure.Norm.R.LatFore{k} + DynamicPPTrials(Selection(j)).Pressure.Norm.R.MedFore{k});
%             Rarch = (DynamicPPTrials(Selection(j)).Pressure.Norm.R.LatArch{k} + DynamicPPTrials(Selection(j)).Pressure.Norm.R.MedArch{k});
%             Rheel = (DynamicPPTrials(Selection(j)).Pressure.Norm.R.LatHeel{k} + DynamicPPTrials(Selection(j)).Pressure.Norm.R.MedHeel{k});
%             plot(x, Rfore, 'b-', 'LineWidth', 1);
%             plot(x, Rarch, 'g-', 'LineWidth', 1);
%             plot(x, Rheel, 'r-', 'LineWidth', 1);
%             plot(x, Rfore + Rarch + Rheel, 'k-','LineWidth',2);
%         end
%     end
%     % create bar series along the top to display timing norms
%     % barh(2.3,100,'b');
%     % barh(2.3,53,'g');
%     % barh(2.3,24,'r');
%     
%     title('Right Foot Pressures');
%     xlabel('% of Stance Period');
%     xlim([0 100]);
%     % ylim([0 2]);
%     legend([MF MA MH],{'Med Fore' 'Med Arch' 'Med Heel'}, 'Location','North');
%     ax = gca;
%     ax.XTick = ([0 25 50 75 100]);
%     ax.XTickLabel = ({'0%','25%','50%','75%','100%'});
%     % ax.YTick = ([0 .25 .50 .75 1 1.25 1.5 1.75 2]);
%     % ax.YTickLabel = [];
%     grid on;
%     
%     subplotsqueeze(TimeSeriesPress, 1.2);
%     % export figure to excel
%     % if strcmp(PPSettings.ExportReport, 'Yes') == 1
%     %     xlsPasteTo('ImpressionsOutput.xlsx','Sheet1',600,600, 'AA17');
%     % end
%     clearvars x LF MF LA MA LH MH Rfore Rarch Rheel Lfore Larch Lheel
% end
% 
% %% Create Area figure
% if strcmp(PPPlots.AllAreas, 'Yes')
%     TimeSeriesArea = figure('Position', [50, 50, 900, 700]);
%     x = 1:100;
%     % LEFT Foot
%     subplot(221);
%     hold on;
%     for j = 1:length(Selection)
%         for k = 1:DynamicPPTrials(Selection(j)).NumFullStepsL
%             LF = plot(x,DynamicPPTrials(Selection(j)).Area.ReSam.L.LatFore{k} ,'b--');
%             plot(x,DynamicPPTrials(Selection(j)).Area.ReSam.L.MedFore{k} ,'b.');
%             LA = plot(x,DynamicPPTrials(Selection(j)).Area.ReSam.L.LatArch{k} ,'g--');
%             plot(x,DynamicPPTrials(Selection(j)).Area.ReSam.L.MedArch{k} ,'g.');
%             LH = plot(x,DynamicPPTrials(Selection(j)).Area.ReSam.L.LatHeel{k} ,'r--');
%             plot(x,DynamicPPTrials(Selection(j)).Area.ReSam.L.MedHeel{k} ,'r.');
%         end
%     end
%     title('Left Foot Areas Med/Lat');
%     ylabel('Area (cm^2)');
%     % ylim([0 1]);
%     ax = gca;
%     % ax.YTick = ([0 .25 .50 .75 1]);
%     % ax.YTickLabel = ({'0%','25%','50%','75%','100%'});
%     legend([LF LA LH],{'Lat Fore' 'Lat Arch' 'Lat Heel'}, 'Location','North');
%     ax.XTick = ([0 25 50 75 100]);
%     ax.XTickLabel = [];
%     grid on;
%     
%     subplot(223);
%     hold on;
%     % plot general heel, arch, fore, and whole foot areas
%     for j = 1:length(Selection)
%         for k = 1:DynamicPPTrials(Selection(j)).NumFullStepsL
%             Lfore = (DynamicPPTrials(Selection(j)).Area.ReSam.L.LatFore{k} + DynamicPPTrials(Selection(j)).Area.ReSam.L.MedFore{k});
%             Larch = (DynamicPPTrials(Selection(j)).Area.ReSam.L.LatArch{k} + DynamicPPTrials(Selection(j)).Area.ReSam.L.MedArch{k});
%             Lheel = (DynamicPPTrials(Selection(j)).Area.ReSam.L.LatHeel{k} + DynamicPPTrials(Selection(j)).Area.ReSam.L.MedHeel{k});
%             plot(x, Lfore, 'b-', 'LineWidth', 1);
%             plot(x, Larch, 'g-', 'LineWidth', 1);
%             plot(x, Lheel, 'r-', 'LineWidth', 1);
%             plot(x, Lfore + Larch + Lheel, 'k-','LineWidth',2);
%         end
%     end
%     % create bar series along the top to display timing norms
%     % barh(2.3,100,'b');
%     % barh(2.3,53,'g');
%     % barh(2.3,24,'r');
%     title('Left Foot Areas');
%     ylabel('Area (cm^2)');
%     % ylim([0 2]);
%     xlabel('% of Stance Period');
%     xlim([0 100]);
%     ax = gca;
%     %ax.XTickLabel = [];
%     % ax.YTick = ([0 .25 .50 .75 1 1.25 1.5 1.75 2]);
%     % ax.YTickLabel = ({'0%','25%','50%','75%','100%','125%','150%','175%','200%'});
%     ax.XTick = ([0 25 50 75 100]);
%     ax.XTickLabel = ({'0%','25%','50%','75%','100%'});
%     grid on;
%     
%     subplot(222);
%     hold on;
%     for j = 1:length(Selection)
%         for k = 1:DynamicPPTrials(Selection(j)).NumFullStepsR
%             plot(x,DynamicPPTrials(Selection(j)).Area.ReSam.R.LatFore{k} ,'b--');
%             MF = plot(x,DynamicPPTrials(Selection(j)).Area.ReSam.R.MedFore{k} ,'b.');
%             plot(x,DynamicPPTrials(Selection(j)).Area.ReSam.R.LatArch{k} ,'g--');
%             MA = plot(x,DynamicPPTrials(Selection(j)).Area.ReSam.R.MedArch{k} ,'g.');
%             plot(x,DynamicPPTrials(Selection(j)).Area.ReSam.R.LatHeel{k} ,'r--');
%             MH = plot(x,DynamicPPTrials(Selection(j)).Area.ReSam.R.MedHeel{k} ,'r.');
%         end
%     end
%     title('Right Foot Areas Med/Lat');
%     % ylim([0 1]);
%     ax = gca;
%     % ax.YTick = ([0 .25 .50 .75 1]);
%     % ax.YTickLabel = [];
%     ax.XTick = ([0 25 50 75 100]);
%     ax.XTickLabel = [];
%     grid on;
%     
%     subplot(224);
%     hold on;
%     for j = 1:length(Selection)
%         for k = 1:DynamicPPTrials(Selection(j)).NumFullStepsR
%             % plot general heel, arch, fore, and whole foot pressures
%             Rfore = (DynamicPPTrials(Selection(j)).Area.ReSam.R.LatFore{k} + DynamicPPTrials(Selection(j)).Area.ReSam.R.MedFore{k});
%             Rarch = (DynamicPPTrials(Selection(j)).Area.ReSam.R.LatArch{k} + DynamicPPTrials(Selection(j)).Area.ReSam.R.MedArch{k});
%             Rheel = (DynamicPPTrials(Selection(j)).Area.ReSam.R.LatHeel{k} + DynamicPPTrials(Selection(j)).Area.ReSam.R.MedHeel{k});
%             plot(x, Rfore, 'b-', 'LineWidth', 1);
%             plot(x, Rarch, 'g-', 'LineWidth', 1);
%             plot(x, Rheel, 'r-', 'LineWidth', 1);
%             plot(x, Rfore + Rarch + Rheel, 'k-','LineWidth',2);
%         end
%     end
%     % create bar series along the top to display timing norms
%     % barh(2.3,100,'b');
%     % barh(2.3,53,'g');
%     % barh(2.3,24,'r');
%     
%     title('Right Foot Areas');
%     xlabel('% of Stance Period');
%     xlim([0 100]);
%     % ylim([0 2]);
%     legend([MF MA MH],{'Med Fore' 'Med Arch' 'Med Heel'}, 'Location','North');
%     ax = gca;
%     ax.XTick = ([0 25 50 75 100]);
%     ax.XTickLabel = ({'0%','25%','50%','75%','100%'});
%     % ax.YTick = ([0 .25 .50 .75 1 1.25 1.5 1.75 2]);
%     % ax.YTickLabel = [];
%     grid on;
%     
%     subplotsqueeze(TimeSeriesArea, 1.2);
%     % export figure to excel
%     % if strcmp(PPSettings.ExportReport, 'Yes') == 1
%     %     xlsPasteTo('ImpressionsOutput.xlsx','Sheet1',600,600, 'AA17');
%     % end
%     clearvars x LF MF LA MA LH MH Rfore Rarch Rheel Lfore Larch Lheel
% end
% 
% %% Compute average and standard deviations of pressures and forces across the various trials
% k = 1; % save all data into side by side columns
% for j = 1:length(Selection)
%     for i = 1:DynamicPPTrials(Selection(j)).NumFullStepsL
%         % Forefoot
%         TrialMean.Pressure.L.LatFore(1:100,k) = DynamicPPTrials(Selection(j)).Pressure.Norm.L.LatFore{i};
%         TrialMean.Pressure.L.MedFore(1:100,k) =  DynamicPPTrials(Selection(j)).Pressure.Norm.L.MedFore{i};
%         TrialMean.Force.L.LatFore(1:100,k) = DynamicPPTrials(Selection(j)).Force.Norm.L.LatFore{i};
%         TrialMean.Force.L.MedFore(1:100,k) = DynamicPPTrials(Selection(j)).Force.Norm.L.MedFore{i};
%         % Arch
%         TrialMean.Pressure.L.LatArch(1:100,k) = DynamicPPTrials(Selection(j)).Pressure.Norm.L.LatArch{i};
%         TrialMean.Pressure.L.MedArch(1:100,k) = DynamicPPTrials(Selection(j)).Pressure.Norm.L.MedArch{i};
%         TrialMean.Force.L.LatArch(1:100,k) = DynamicPPTrials(Selection(j)).Force.Norm.L.LatArch{i};
%         TrialMean.Force.L.MedArch(1:100,k) = DynamicPPTrials(Selection(j)).Force.Norm.L.MedArch{i};
%         % Heel
%         TrialMean.Pressure.L.LatHeel(1:100,k) = DynamicPPTrials(Selection(j)).Pressure.Norm.L.LatHeel{i};
%         TrialMean.Pressure.L.MedHeel(1:100,k) = DynamicPPTrials(Selection(j)).Pressure.Norm.L.MedHeel{i};
%         TrialMean.Force.L.LatHeel(1:100,k) = DynamicPPTrials(Selection(j)).Force.Norm.L.LatHeel{i};
%         TrialMean.Force.L.MedHeel(1:100,k) = DynamicPPTrials(Selection(j)).Force.Norm.L.MedHeel{i};
%         % Whole
%         TrialMean.Pressure.L.WholeFoot(1:100,k) = DynamicPPTrials(Selection(j)).Pressure.Norm.L.WholeFoot{i};
%         TrialMean.Force.L.WholeFoot(1:100,k) = DynamicPPTrials(Selection(j)).Force.Norm.L.WholeFoot{i};
%         k = k+1;
%         if j == Trial2Cont(1) && i == Lcont
%             % Forefoot
%             TrialSelection.Pressure.L.LatFore.Avg = nanmean(DynamicPPTrials(Selection(j)).Pressure.Norm.L.LatFore{i});
%             TrialSelection.Pressure.L.LatFore.Std = nanstd(DynamicPPTrials(Selection(j)).Pressure.Norm.L.LatFore{i});
%             TrialSelection.Pressure.L.MedFore.Avg = nanmean(DynamicPPTrials(Selection(j)).Pressure.Norm.L.MedFore{i});
%             TrialSelection.Pressure.L.MedFore.Std = nanstd(DynamicPPTrials(Selection(j)).Pressure.Norm.L.MedFore{i});
%             TrialSelection.Force.L.LatFore.Avg = nanmean(DynamicPPTrials(Selection(j)).Force.Norm.L.LatFore{i});
%             TrialSelection.Force.L.LatFore.Std = nanstd(DynamicPPTrials(Selection(j)).Force.Norm.L.LatFore{i});
%             TrialSelection.Force.L.MedFore.Avg = nanmean(DynamicPPTrials(Selection(j)).Force.Norm.L.MedFore{i});
%             TrialSelection.Force.L.MedFore.Std = nanstd(DynamicPPTrials(Selection(j)).Force.Norm.L.MedFore{i});
%             % arch
%             TrialSelection.Pressure.L.LatArch.Avg = nanmean(DynamicPPTrials(Selection(j)).Pressure.Norm.L.LatArch{i});
%             TrialSelection.Pressure.L.LatArch.Std = nanstd(DynamicPPTrials(Selection(j)).Pressure.Norm.L.LatArch{i});
%             TrialSelection.Pressure.L.MedArch.Avg = nanmean(DynamicPPTrials(Selection(j)).Pressure.Norm.L.MedArch{i});
%             TrialSelection.Pressure.L.MedArch.Std = nanstd(DynamicPPTrials(Selection(j)).Pressure.Norm.L.MedArch{i});
%             TrialSelection.Force.L.LatArch.Avg = nanmean(DynamicPPTrials(Selection(j)).Force.Norm.L.LatArch{i});
%             TrialSelection.Force.L.LatArch.Std = nanstd(DynamicPPTrials(Selection(j)).Force.Norm.L.LatArch{i});
%             TrialSelection.Force.L.MedArch.Avg = nanmean(DynamicPPTrials(Selection(j)).Force.Norm.L.MedArch{i});
%             TrialSelection.Force.L.MedArch.Std = nanstd(DynamicPPTrials(Selection(j)).Force.Norm.L.MedArch{i});
%             % heel
%             TrialSelection.Pressure.L.LatHeel.Avg = nanmean(DynamicPPTrials(Selection(j)).Pressure.Norm.L.LatHeel{i});
%             TrialSelection.Pressure.L.LatHeel.Std = nanstd(DynamicPPTrials(Selection(j)).Pressure.Norm.L.LatHeel{i});
%             TrialSelection.Pressure.L.MedHeel.Avg = nanmean(DynamicPPTrials(Selection(j)).Pressure.Norm.L.MedHeel{i});
%             TrialSelection.Pressure.L.MedHeel.Std = nanstd(DynamicPPTrials(Selection(j)).Pressure.Norm.L.MedHeel{i});
%             TrialSelection.Force.L.LatHeel.Avg = nanmean(DynamicPPTrials(Selection(j)).Force.Norm.L.LatHeel{i});
%             TrialSelection.Force.L.LatHeel.Std = nanstd(DynamicPPTrials(Selection(j)).Force.Norm.L.LatHeel{i});
%             TrialSelection.Force.L.MedHeel.Avg = nanmean(DynamicPPTrials(Selection(j)).Force.Norm.L.MedHeel{i});
%             TrialSelection.Force.L.MedHeel.Std = nanstd(DynamicPPTrials(Selection(j)).Force.Norm.L.MedHeel{i});
%             % WholeFoot
%             TrialSelection.Pressure.L.WholeFoot.Avg = nanmean(DynamicPPTrials(Selection(j)).Pressure.Norm.L.WholeFoot{i});
%             TrialSelection.Pressure.L.WholeFoot.Std = nanstd(DynamicPPTrials(Selection(j)).Pressure.Norm.L.WholeFoot{i});
%             TrialSelection.Force.L.WholeFoot.Avg = nanmean(DynamicPPTrials(Selection(j)).Force.Norm.L.WholeFoot{i});
%             TrialSelection.Force.L.WholeFoot.Std = nanstd(DynamicPPTrials(Selection(j)).Force.Norm.L.WholeFoot{i});
%         end
%     end
% end
% k = 1;
% for j = 1:length(Selection)
%     % RIGHT
%     for i = 1:DynamicPPTrials(Selection(j)).NumFullStepsR
%         % Forefoot
%         TrialMean.Pressure.R.LatFore(1:100,k) = DynamicPPTrials(Selection(j)).Pressure.Norm.R.LatFore{i};
%         TrialMean.Pressure.R.MedFore(1:100,k) = DynamicPPTrials(Selection(j)).Pressure.Norm.R.MedFore{i};
%         TrialMean.Force.R.LatFore(1:100,k) = DynamicPPTrials(Selection(j)).Force.Norm.R.LatFore{i};
%         TrialMean.Force.R.MedFore(1:100,k) = DynamicPPTrials(Selection(j)).Force.Norm.R.MedFore{i};
%         % Arch
%         TrialMean.Pressure.R.LatArch(1:100,k) = DynamicPPTrials(Selection(j)).Pressure.Norm.R.LatArch{i};
%         TrialMean.Pressure.R.MedArch(1:100,k) = DynamicPPTrials(Selection(j)).Pressure.Norm.R.MedArch{i};
%         TrialMean.Force.R.LatArch(1:100,k) = DynamicPPTrials(Selection(j)).Force.Norm.R.LatArch{i};
%         TrialMean.Force.R.MedArch(1:100,k) = DynamicPPTrials(Selection(j)).Force.Norm.R.MedArch{i};
%         % Heel
%         TrialMean.Pressure.R.LatHeel(1:100,k) = DynamicPPTrials(Selection(j)).Pressure.Norm.R.LatHeel{i};
%         TrialMean.Pressure.R.MedHeel(1:100,k) = DynamicPPTrials(Selection(j)).Pressure.Norm.R.MedHeel{i};
%         TrialMean.Force.R.LatHeel(1:100,k) = DynamicPPTrials(Selection(j)).Force.Norm.R.LatHeel{i};
%         TrialMean.Force.R.MedHeel(1:100,k) = DynamicPPTrials(Selection(j)).Force.Norm.R.MedHeel{i};
%         % Whole
%         TrialMean.Pressure.R.WholeFoot(1:100,k) = DynamicPPTrials(Selection(j)).Pressure.Norm.R.WholeFoot{i};
%         TrialMean.Force.R.WholeFoot(1:100,k) = DynamicPPTrials(Selection(j)).Force.Norm.R.WholeFoot{i};
%         k = k+1;
%         if j == Trial2Cont(2) && i == Rcont
%             % Forefoot
%             TrialSelection.Pressure.R.LatFore.Avg = nanmean(DynamicPPTrials(Selection(j)).Pressure.Norm.R.LatFore{i});
%             TrialSelection.Pressure.R.LatFore.Std = nanstd(DynamicPPTrials(Selection(j)).Pressure.Norm.R.LatFore{i});
%             TrialSelection.Pressure.R.MedFore.Avg = nanmean(DynamicPPTrials(Selection(j)).Pressure.Norm.R.MedFore{i});
%             TrialSelection.Pressure.R.MedFore.Std = nanstd(DynamicPPTrials(Selection(j)).Pressure.Norm.R.MedFore{i});
%             TrialSelection.Force.R.LatFore.Avg = nanmean(DynamicPPTrials(Selection(j)).Force.Norm.R.LatFore{i});
%             TrialSelection.Force.R.LatFore.Std = nanstd(DynamicPPTrials(Selection(j)).Force.Norm.R.LatFore{i});
%             TrialSelection.Force.R.MedFore.Avg = nanmean(DynamicPPTrials(Selection(j)).Force.Norm.R.MedFore{i});
%             TrialSelection.Force.R.MedFore.Std = nanstd(DynamicPPTrials(Selection(j)).Force.Norm.R.MedFore{i});
%             % arch
%             TrialSelection.Pressure.R.LatArch.Avg = nanmean(DynamicPPTrials(Selection(j)).Pressure.Norm.R.LatArch{i});
%             TrialSelection.Pressure.R.LatArch.Std = nanstd(DynamicPPTrials(Selection(j)).Pressure.Norm.R.LatArch{i});
%             TrialSelection.Pressure.R.MedArch.Avg = nanmean(DynamicPPTrials(Selection(j)).Pressure.Norm.R.MedArch{i});
%             TrialSelection.Pressure.R.MedArch.Std = nanstd(DynamicPPTrials(Selection(j)).Pressure.Norm.R.MedArch{i});
%             TrialSelection.Force.R.LatArch.Avg = nanmean(DynamicPPTrials(Selection(j)).Force.Norm.R.LatArch{i});
%             TrialSelection.Force.R.LatArch.Std = nanstd(DynamicPPTrials(Selection(j)).Force.Norm.R.LatArch{i});
%             TrialSelection.Force.R.MedArch.Avg = nanmean(DynamicPPTrials(Selection(j)).Force.Norm.R.MedArch{i});
%             TrialSelection.Force.R.MedArch.Std = nanstd(DynamicPPTrials(Selection(j)).Force.Norm.R.MedArch{i});
%             % heel
%             TrialSelection.Pressure.R.LatHeel.Avg = nanmean(DynamicPPTrials(Selection(j)).Pressure.Norm.R.LatHeel{i});
%             TrialSelection.Pressure.R.LatHeel.Std = nanstd(DynamicPPTrials(Selection(j)).Pressure.Norm.R.LatHeel{i});
%             TrialSelection.Pressure.R.MedHeel.Avg = nanmean(DynamicPPTrials(Selection(j)).Pressure.Norm.R.MedHeel{i});
%             TrialSelection.Pressure.R.MedHeel.Std = nanstd(DynamicPPTrials(Selection(j)).Pressure.Norm.R.MedHeel{i});
%             TrialSelection.Force.R.LatHeel.Avg = nanmean(DynamicPPTrials(Selection(j)).Force.Norm.R.LatHeel{i});
%             TrialSelection.Force.R.LatHeel.Std = nanstd(DynamicPPTrials(Selection(j)).Force.Norm.R.LatHeel{i});
%             TrialSelection.Force.R.MedHeel.Avg = nanmean(DynamicPPTrials(Selection(j)).Force.Norm.R.MedHeel{i});
%             TrialSelection.Force.R.MedHeel.Std = nanstd(DynamicPPTrials(Selection(j)).Force.Norm.R.MedHeel{i});
%             % WholeFoot
%             TrialSelection.Pressure.R.WholeFoot.Avg = nanmean(DynamicPPTrials(Selection(j)).Pressure.Norm.R.WholeFoot{i});
%             TrialSelection.Pressure.R.WholeFoot.Std = nanstd(DynamicPPTrials(Selection(j)).Pressure.Norm.R.WholeFoot{i});
%             TrialSelection.Force.R.WholeFoot.Avg = nanmean(DynamicPPTrials(Selection(j)).Force.Norm.R.WholeFoot{i});
%             TrialSelection.Force.R.WholeFoot.Std = nanstd(DynamicPPTrials(Selection(j)).Force.Norm.R.WholeFoot{i});
%         end
%     end
% end
% % compute means and STDs
% % LEFT
% % pressure
% Trial_Mean.Pressure.L.LatFore.Avg = mean(TrialMean.Pressure.L.LatFore, 2);  % Fore
% Trial_Mean.Pressure.L.LatFore.STD = std(TrialMean.Pressure.L.LatFore,0, 2);
% Trial_Mean.Pressure.L.MedFore.Avg = mean(TrialMean.Pressure.L.MedFore, 2);
% Trial_Mean.Pressure.L.MedFore.STD = std(TrialMean.Pressure.L.MedFore,0, 2);
% Trial_Mean.Pressure.L.LatArch.Avg = mean(TrialMean.Pressure.L.LatArch, 2);  % Arch
% Trial_Mean.Pressure.L.LatArch.STD = std(TrialMean.Pressure.L.LatArch,0, 2);
% Trial_Mean.Pressure.L.MedArch.Avg = mean(TrialMean.Pressure.L.MedArch, 2);
% Trial_Mean.Pressure.L.MedArch.STD = std(TrialMean.Pressure.L.MedArch,0, 2);
% Trial_Mean.Pressure.L.LatHeel.Avg = mean(TrialMean.Pressure.L.LatHeel, 2);  % Heel
% Trial_Mean.Pressure.L.LatHeel.STD = std(TrialMean.Pressure.L.LatHeel,0, 2);
% Trial_Mean.Pressure.L.MedHeel.Avg = mean(TrialMean.Pressure.L.MedHeel, 2);
% Trial_Mean.Pressure.L.MedHeel.STD = std(TrialMean.Pressure.L.MedHeel,0, 2);
% Trial_Mean.Pressure.L.WholeFoot.Avg = mean(TrialMean.Pressure.L.WholeFoot, 2); % whole
% Trial_Mean.Pressure.L.WholeFoot.STD = std(TrialMean.Pressure.L.WholeFoot,0, 2);
% % force
% Trial_Mean.Force.L.LatFore.Avg = mean(TrialMean.Force.L.LatFore, 2);  % Fore
% Trial_Mean.Force.L.LatFore.STD = std(TrialMean.Force.L.LatFore,0, 2);
% Trial_Mean.Force.L.MedFore.Avg = mean(TrialMean.Force.L.MedFore, 2);
% Trial_Mean.Force.L.MedFore.STD = std(TrialMean.Force.L.MedFore,0, 2);
% Trial_Mean.Force.L.LatArch.Avg = mean(TrialMean.Force.L.LatArch, 2);  % Arch
% Trial_Mean.Force.L.LatArch.STD = std(TrialMean.Force.L.LatArch,0, 2);
% Trial_Mean.Force.L.MedArch.Avg = mean(TrialMean.Force.L.MedArch, 2);
% Trial_Mean.Force.L.MedArch.STD = std(TrialMean.Force.L.MedArch,0, 2);
% Trial_Mean.Force.L.LatHeel.Avg = mean(TrialMean.Force.L.LatHeel, 2);  % Heel
% Trial_Mean.Force.L.LatHeel.STD = std(TrialMean.Force.L.LatHeel,0, 2);
% Trial_Mean.Force.L.MedHeel.Avg = mean(TrialMean.Force.L.MedHeel, 2);
% Trial_Mean.Force.L.MedHeel.STD = std(TrialMean.Force.L.MedHeel,0, 2);
% Trial_Mean.Force.L.WholeFoot.Avg = mean(TrialMean.Force.L.WholeFoot, 2); % whole
% Trial_Mean.Force.L.WholeFoot.STD = std(TrialMean.Force.L.WholeFoot,0, 2);
% % Right
% % pressure
% Trial_Mean.Pressure.R.LatFore.Avg = mean(TrialMean.Pressure.R.LatFore, 2);  % Fore
% Trial_Mean.Pressure.R.LatFore.STD = std(TrialMean.Pressure.R.LatFore,0, 2);
% Trial_Mean.Pressure.R.MedFore.Avg = mean(TrialMean.Pressure.R.MedFore, 2);
% Trial_Mean.Pressure.R.MedFore.STD = std(TrialMean.Pressure.R.MedFore,0, 2);
% Trial_Mean.Pressure.R.LatArch.Avg = mean(TrialMean.Pressure.R.LatArch, 2);  % Arch
% Trial_Mean.Pressure.R.LatArch.STD = std(TrialMean.Pressure.R.LatArch,0, 2);
% Trial_Mean.Pressure.R.MedArch.Avg = mean(TrialMean.Pressure.R.MedArch, 2);
% Trial_Mean.Pressure.R.MedArch.STD = std(TrialMean.Pressure.R.MedArch,0, 2);
% Trial_Mean.Pressure.R.LatHeel.Avg = mean(TrialMean.Pressure.R.LatHeel, 2);  % Heel
% Trial_Mean.Pressure.R.LatHeel.STD = std(TrialMean.Pressure.R.LatHeel,0, 2);
% Trial_Mean.Pressure.R.MedHeel.Avg = mean(TrialMean.Pressure.R.MedHeel, 2);
% Trial_Mean.Pressure.R.MedHeel.STD = std(TrialMean.Pressure.R.MedHeel,0, 2);
% Trial_Mean.Pressure.R.WholeFoot.Avg = mean(TrialMean.Pressure.R.WholeFoot, 2); % whole
% Trial_Mean.Pressure.R.WholeFoot.STD = std(TrialMean.Pressure.R.WholeFoot,0, 2);
% % force
% Trial_Mean.Force.R.LatFore.Avg = mean(TrialMean.Force.R.LatFore, 2);  % Fore
% Trial_Mean.Force.R.LatFore.STD = std(TrialMean.Force.R.LatFore,0, 2);
% Trial_Mean.Force.R.MedFore.Avg = mean(TrialMean.Force.R.MedFore, 2);
% Trial_Mean.Force.R.MedFore.STD = std(TrialMean.Force.R.MedFore,0, 2);
% Trial_Mean.Force.R.LatArch.Avg = mean(TrialMean.Force.R.LatArch, 2);  % Arch
% Trial_Mean.Force.R.LatArch.STD = std(TrialMean.Force.R.LatArch,0, 2);
% Trial_Mean.Force.R.MedArch.Avg = mean(TrialMean.Force.R.MedArch, 2);
% Trial_Mean.Force.R.MedArch.STD = std(TrialMean.Force.R.MedArch,0, 2);
% Trial_Mean.Force.R.LatHeel.Avg = mean(TrialMean.Force.R.LatHeel, 2);  % Heel
% Trial_Mean.Force.R.LatHeel.STD = std(TrialMean.Force.R.LatHeel,0, 2);
% Trial_Mean.Force.R.MedHeel.Avg = mean(TrialMean.Force.R.MedHeel, 2);
% Trial_Mean.Force.R.MedHeel.STD = std(TrialMean.Force.R.MedHeel,0, 2);
% Trial_Mean.Force.R.WholeFoot.Avg = mean(TrialMean.Force.R.WholeFoot, 2); % whole
% Trial_Mean.Force.R.WholeFoot.STD = std(TrialMean.Force.R.WholeFoot,0, 2);
% 
% %% plot means and STDs
% x = 1:100;
% if strcmp(PPPlots.AvgPressures, 'Yes')
%     figure;  % Pressures
%     for j = 1: length(Selection)
%         for i = 1:DynamicPPTrials(Selection(j)).NumFullStepsL
%             subplot(421); hold on; title('Average Forefoot Medial Pressures'); hline(0,'k');
%             plot(x, Trial_Mean.Pressure.L.MedFore.Avg ,'r-', 'LineWidth',2');
%             plot(x, Trial_Mean.Pressure.L.MedFore.Avg + Trial_Mean.Pressure.L.MedFore.STD ,'r--', 'LineWidth',1');
%             plot(x, Trial_Mean.Pressure.L.MedFore.Avg - Trial_Mean.Pressure.L.MedFore.STD,'r--', 'LineWidth',1');
%             subplot(422); hold on; title('Average Forefoot Lateral Pressures'); hline(0,'k');
%             plot(x, Trial_Mean.Pressure.L.LatFore.Avg ,'r-', 'LineWidth',2');
%             plot(x, Trial_Mean.Pressure.L.LatFore.Avg + Trial_Mean.Pressure.L.LatFore.STD ,'r--', 'LineWidth',1');
%             plot(x, Trial_Mean.Pressure.L.LatFore.Avg - Trial_Mean.Pressure.L.LatFore.STD,'r--', 'LineWidth',1');
%             subplot(423); hold on; title('Average Midfoot Medial Pressures'); hline(0,'k');
%             plot(x, Trial_Mean.Pressure.L.MedArch.Avg ,'r-', 'LineWidth',2');
%             plot(x, Trial_Mean.Pressure.L.MedArch.Avg + Trial_Mean.Pressure.L.MedArch.STD ,'r--', 'LineWidth',1');
%             plot(x, Trial_Mean.Pressure.L.MedArch.Avg - Trial_Mean.Pressure.L.MedArch.STD,'r--', 'LineWidth',1');
%             subplot(424); hold on; title('Average Midfoot Lateral Pressures'); hline(0,'k');
%             plot(x, Trial_Mean.Pressure.L.LatArch.Avg ,'r-', 'LineWidth',2');
%             plot(x, Trial_Mean.Pressure.L.LatArch.Avg + Trial_Mean.Pressure.L.LatArch.STD ,'r--', 'LineWidth',1');
%             plot(x, Trial_Mean.Pressure.L.LatArch.Avg - Trial_Mean.Pressure.L.LatArch.STD,'r--', 'LineWidth',1');
%             subplot(425); hold on; title('Average Hindfoot Medial Pressures'); hline(0,'k');
%             plot(x, Trial_Mean.Pressure.L.MedHeel.Avg ,'r-', 'LineWidth',2');
%             plot(x, Trial_Mean.Pressure.L.MedHeel.Avg + Trial_Mean.Pressure.L.MedHeel.STD ,'r--', 'LineWidth',1');
%             plot(x, Trial_Mean.Pressure.L.MedHeel.Avg - Trial_Mean.Pressure.L.MedHeel.STD,'r--', 'LineWidth',1');
%             subplot(426); hold on; title('Average Hindfoot Lateral Pressures'); hline(0,'k');
%             plot(x, Trial_Mean.Pressure.L.LatHeel.Avg ,'r-', 'LineWidth',2');
%             plot(x, Trial_Mean.Pressure.L.LatHeel.Avg + Trial_Mean.Pressure.L.LatHeel.STD ,'r--', 'LineWidth',1');
%             plot(x, Trial_Mean.Pressure.L.LatHeel.Avg - Trial_Mean.Pressure.L.LatHeel.STD,'r--', 'LineWidth',1');
%             subplot(427); hold on; title('Average Whole Foot Pressures'); hline(0,'k');
%             plot(x, Trial_Mean.Pressure.L.WholeFoot.Avg ,'r-', 'LineWidth',2');
%             plot(x, Trial_Mean.Pressure.L.WholeFoot.Avg + Trial_Mean.Pressure.L.WholeFoot.STD ,'r--', 'LineWidth',1');
%             plot(x, Trial_Mean.Pressure.L.WholeFoot.Avg - Trial_Mean.Pressure.L.WholeFoot.STD,'r--', 'LineWidth',1');
%         end
%         for i = 1:DynamicPPTrials(Selection(j)).NumFullStepsR
%             subplot(421); hold on; title('Average Forefoot Medial Pressures'); hline(0,'k');grid on;
%             plot(x, Trial_Mean.Pressure.R.MedFore.Avg ,'g-', 'LineWidth',2');
%             plot(x, Trial_Mean.Pressure.R.MedFore.Avg + Trial_Mean.Pressure.R.MedFore.STD ,'g--', 'LineWidth',1');
%             plot(x, Trial_Mean.Pressure.R.MedFore.Avg - Trial_Mean.Pressure.R.MedFore.STD,'g--', 'LineWidth',1');
%             subplot(422); hold on; title('Average Forefoot Lateral Pressures'); hline(0,'k');grid on;
%             plot(x, Trial_Mean.Pressure.R.LatFore.Avg ,'g-', 'LineWidth',2');
%             plot(x, Trial_Mean.Pressure.R.LatFore.Avg + Trial_Mean.Pressure.R.LatFore.STD ,'g--', 'LineWidth',1');
%             plot(x, Trial_Mean.Pressure.R.LatFore.Avg - Trial_Mean.Pressure.R.LatFore.STD,'g--', 'LineWidth',1');
%             subplot(423); hold on; title('Average Midfoot Medial Pressures'); hline(0,'k');grid on;
%             plot(x, Trial_Mean.Pressure.R.MedArch.Avg ,'g-', 'LineWidth',2');
%             plot(x, Trial_Mean.Pressure.R.MedArch.Avg + Trial_Mean.Pressure.R.MedArch.STD ,'g--', 'LineWidth',1');
%             plot(x, Trial_Mean.Pressure.R.MedArch.Avg - Trial_Mean.Pressure.R.MedArch.STD,'g--', 'LineWidth',1');
%             subplot(424); hold on; title('Average Midfoot Lateral Pressures'); hline(0,'k');grid on;
%             plot(x, Trial_Mean.Pressure.R.LatArch.Avg ,'g-', 'LineWidth',2');
%             plot(x, Trial_Mean.Pressure.R.LatArch.Avg + Trial_Mean.Pressure.R.LatArch.STD ,'g--', 'LineWidth',1');
%             plot(x, Trial_Mean.Pressure.R.LatArch.Avg - Trial_Mean.Pressure.R.LatArch.STD,'g--', 'LineWidth',1');
%             subplot(425); hold on; title('Average Hindfoot Medial Pressures'); hline(0,'k');grid on;
%             plot(x, Trial_Mean.Pressure.R.MedHeel.Avg ,'g-', 'LineWidth',2');
%             plot(x, Trial_Mean.Pressure.R.MedHeel.Avg + Trial_Mean.Pressure.R.MedHeel.STD ,'g--', 'LineWidth',1');
%             plot(x, Trial_Mean.Pressure.R.MedHeel.Avg - Trial_Mean.Pressure.R.MedHeel.STD,'g--', 'LineWidth',1');
%             subplot(426); hold on; title('Average Hindfoot Lateral Pressures'); hline(0,'k');grid on;
%             plot(x, Trial_Mean.Pressure.R.LatHeel.Avg ,'g-', 'LineWidth',2');
%             plot(x, Trial_Mean.Pressure.R.LatHeel.Avg + Trial_Mean.Pressure.R.LatHeel.STD ,'g--', 'LineWidth',1');
%             plot(x, Trial_Mean.Pressure.R.LatHeel.Avg - Trial_Mean.Pressure.R.LatHeel.STD,'g--', 'LineWidth',1');
%             subplot(427); hold on; title('Average Whole Foot Pressures'); hline(0,'k'); grid on;
%             plot(x, Trial_Mean.Pressure.R.WholeFoot.Avg ,'g-', 'LineWidth',2');
%             plot(x, Trial_Mean.Pressure.R.WholeFoot.Avg + Trial_Mean.Pressure.R.WholeFoot.STD ,'g--', 'LineWidth',1');
%             plot(x, Trial_Mean.Pressure.R.WholeFoot.Avg - Trial_Mean.Pressure.R.WholeFoot.STD,'g--', 'LineWidth',1');
%         end
%     end
% end
% if strcmp(PPPlots.AvgForces, 'Yes')
%     figure;  % Forces
%     for j = 1: length(Selection)
%         for i = 1:DynamicPPTrials(Selection(j)).NumFullStepsL
%             subplot(421); hold on; title('Average Forefoot Medial Forces'); hline(0,'k');
%             plot(x, Trial_Mean.Force.L.MedFore.Avg ,'r-', 'LineWidth',2');
%             plot(x, Trial_Mean.Force.L.MedFore.Avg + Trial_Mean.Force.L.MedFore.STD ,'r--', 'LineWidth',1');
%             plot(x, Trial_Mean.Force.L.MedFore.Avg - Trial_Mean.Force.L.MedFore.STD,'r--', 'LineWidth',1');
%             subplot(422); hold on; title('Average Forefoot Lateral Forces'); hline(0,'k');
%             plot(x, Trial_Mean.Force.L.LatFore.Avg ,'r-', 'LineWidth',2');
%             plot(x, Trial_Mean.Force.L.LatFore.Avg + Trial_Mean.Force.L.LatFore.STD ,'r--', 'LineWidth',1');
%             plot(x, Trial_Mean.Force.L.LatFore.Avg - Trial_Mean.Force.L.LatFore.STD,'r--', 'LineWidth',1');
%             subplot(423); hold on; title('Average Midfoot Medial Forces'); hline(0,'k');
%             plot(x, Trial_Mean.Force.L.MedArch.Avg ,'r-', 'LineWidth',2');
%             plot(x, Trial_Mean.Force.L.MedArch.Avg + Trial_Mean.Force.L.MedArch.STD ,'r--', 'LineWidth',1');
%             plot(x, Trial_Mean.Force.L.MedArch.Avg - Trial_Mean.Force.L.MedArch.STD,'r--', 'LineWidth',1');
%             subplot(424); hold on; title('Average Midfoot Lateral Forces'); hline(0,'k');
%             plot(x, Trial_Mean.Force.L.LatArch.Avg ,'r-', 'LineWidth',2');
%             plot(x, Trial_Mean.Force.L.LatArch.Avg + Trial_Mean.Force.L.LatArch.STD ,'r--', 'LineWidth',1');
%             plot(x, Trial_Mean.Force.L.LatArch.Avg - Trial_Mean.Force.L.LatArch.STD,'r--', 'LineWidth',1');
%             subplot(425); hold on; title('Average Hindfoot Medial Forces'); hline(0,'k');
%             plot(x, Trial_Mean.Force.L.MedHeel.Avg ,'r-', 'LineWidth',2');
%             plot(x, Trial_Mean.Force.L.MedHeel.Avg + Trial_Mean.Force.L.MedHeel.STD ,'r--', 'LineWidth',1');
%             plot(x, Trial_Mean.Force.L.MedHeel.Avg - Trial_Mean.Force.L.MedHeel.STD,'r--', 'LineWidth',1');
%             subplot(426); hold on; title('Average Hindfoot Lateral Forces'); hline(0,'k');
%             plot(x, Trial_Mean.Force.L.LatHeel.Avg ,'r-', 'LineWidth',2');
%             plot(x, Trial_Mean.Force.L.LatHeel.Avg + Trial_Mean.Force.L.LatHeel.STD ,'r--', 'LineWidth',1');
%             plot(x, Trial_Mean.Force.L.LatHeel.Avg - Trial_Mean.Force.L.LatHeel.STD,'r--', 'LineWidth',1');
%             subplot(427); hold on; title('Average Whole Foot Forces'); hline(0,'k');
%             plot(x, Trial_Mean.Force.L.WholeFoot.Avg ,'r-', 'LineWidth',2');
%             plot(x, Trial_Mean.Force.L.WholeFoot.Avg + Trial_Mean.Force.L.WholeFoot.STD ,'r--', 'LineWidth',1');
%             plot(x, Trial_Mean.Force.L.WholeFoot.Avg - Trial_Mean.Force.L.WholeFoot.STD,'r--', 'LineWidth',1');
%         end
%         for i = 1:DynamicPPTrials(Selection(j)).NumFullStepsR
%             subplot(421); hold on; title('Average Forefoot Medial Forces'); hline(0,'k'); grid on;
%             plot(x, Trial_Mean.Force.R.MedFore.Avg ,'g-', 'LineWidth',2');
%             plot(x, Trial_Mean.Force.R.MedFore.Avg + Trial_Mean.Force.R.MedFore.STD ,'g--', 'LineWidth',1');
%             plot(x, Trial_Mean.Force.R.MedFore.Avg - Trial_Mean.Force.R.MedFore.STD,'g--', 'LineWidth',1');
%             subplot(422); hold on; title('Average Forefoot Lateral Forces'); hline(0,'k'); grid on;
%             plot(x, Trial_Mean.Force.R.LatFore.Avg ,'g-', 'LineWidth',2');
%             plot(x, Trial_Mean.Force.R.LatFore.Avg + Trial_Mean.Force.R.LatFore.STD ,'g--', 'LineWidth',1');
%             plot(x, Trial_Mean.Force.R.LatFore.Avg - Trial_Mean.Force.R.LatFore.STD,'g--', 'LineWidth',1');
%             subplot(423); hold on; title('Average Midfoot Medial Forces'); hline(0,'k'); grid on;
%             plot(x, Trial_Mean.Force.R.MedArch.Avg ,'g-', 'LineWidth',2');
%             plot(x, Trial_Mean.Force.R.MedArch.Avg + Trial_Mean.Force.R.MedArch.STD ,'g--', 'LineWidth',1');
%             plot(x, Trial_Mean.Force.R.MedArch.Avg - Trial_Mean.Force.R.MedArch.STD,'g--', 'LineWidth',1');
%             subplot(424); hold on; title('Average Midfoot Lateral Forces'); hline(0,'k'); grid on;
%             plot(x, Trial_Mean.Force.R.LatArch.Avg ,'g-', 'LineWidth',2');
%             plot(x, Trial_Mean.Force.R.LatArch.Avg + Trial_Mean.Force.R.LatArch.STD ,'g--', 'LineWidth',1');
%             plot(x, Trial_Mean.Force.R.LatArch.Avg - Trial_Mean.Force.R.LatArch.STD,'g--', 'LineWidth',1');
%             subplot(425); hold on; title('Average Hindfoot Medial Forces'); hline(0,'k'); grid on;
%             plot(x, Trial_Mean.Force.R.MedHeel.Avg ,'g-', 'LineWidth',2');
%             plot(x, Trial_Mean.Force.R.MedHeel.Avg + Trial_Mean.Force.R.MedHeel.STD ,'g--', 'LineWidth',1');
%             plot(x, Trial_Mean.Force.R.MedHeel.Avg - Trial_Mean.Force.R.MedHeel.STD,'g--', 'LineWidth',1');
%             subplot(426); hold on; title('Average Hindfoot Lateral Forces'); hline(0,'k'); grid on;
%             plot(x, Trial_Mean.Force.R.LatHeel.Avg ,'g-', 'LineWidth',2');
%             plot(x, Trial_Mean.Force.R.LatHeel.Avg + Trial_Mean.Force.R.LatHeel.STD ,'g--', 'LineWidth',1');
%             plot(x, Trial_Mean.Force.R.LatHeel.Avg - Trial_Mean.Force.R.LatHeel.STD,'g--', 'LineWidth',1');
%             subplot(427); hold on; title('Average Whole Foot Forces'); hline(0,'k'); grid on;
%             plot(x, Trial_Mean.Force.R.WholeFoot.Avg ,'g-', 'LineWidth',2');
%             plot(x, Trial_Mean.Force.R.WholeFoot.Avg + Trial_Mean.Force.R.WholeFoot.STD ,'g--', 'LineWidth',1');
%             plot(x, Trial_Mean.Force.R.WholeFoot.Avg - Trial_Mean.Force.R.WholeFoot.STD,'g--', 'LineWidth',1');
%         end
%     end
% end
% t-test to determine differences between L and R sides?

% find instance of peak loading pressure, peak propulsion pressure, and overall peak pressure
% Thresh = 0.10; % peak force threshold set to 10% of body weight
% for j = 1:length(Selection)
%     % LEFT Foot
%     for i = 1:DynamicPPTrials(Selection(j)).NumFullStepsL
%         % Peak Force, Pressure, and Area
%         % minimum promenance as 10% of weight sorted in descending order
%         [DynamicPPTrials(Selection(j)).Peak.Left(i).Val, DynamicPPTrials(Selection(j)).Peak.Left(i).Ind] = findpeaks(DynamicPPTrials(Selection(j)).Force.ReSam.L.WholeFoot{i},...
%             'MinPeakProminence',Subject.Newtons * Thresh, 'SortStr','descend');
% %         if strcmp(PPSettings.ToeWalker,'No') == 1
%             %Loading force peak from Heel regions
%             % find instance of loading peak force and pressure
%             [DynamicPPTrials(Selection(j)).LoadPeak.Left(i).Val, DynamicPPTrials(Selection(j)).LoadPeak.Left(i).Ind] = findpeaks(DynamicPPTrials(Selection(j)).Force.ReSam.L.LatHeel{i} +...
%                 DynamicPPTrials(Selection(j)).Force.ReSam.L.MedHeel{i} , 'MinPeakProminence',Subject.Newtons * Thresh, 'SortStr','descend');
%             %             [DynamicPPTrials(Selection(j)).LoadPeak.Pressure.Left.Val{i}, DynamicPPTrials(Selection(j)).LoadPeak.Pressure.Left.Ind{i}] = findpeaks(DynamicPPTrials(Selection(j)).Pressure.ReSam.L.LatHeel{i} +...
%             %                 DynamicPPTrials(Selection(j)).Pressure.ReSam.L.MedHeel{i} , 'MinPeakProminence',Subject.Newtons * Thresh, 'SortStr','descend');
%             %             % find whole foot force and pressure at the point of max heel force
%             DynamicPPTrials(Selection(j)).LoadPeak.Left(i).Val = DynamicPPTrials(Selection(j)).Force.ReSam.L.WholeFoot{i}(DynamicPPTrials(Selection(j)).LoadPeak.Left(i).Ind);
%             % propulsion force peak from forefoot
%             [DynamicPPTrials(Selection(j)).PropPeak.Left(i).Val, DynamicPPTrials(Selection(j)).PropPeak.Left(i).Ind] = findpeaks(DynamicPPTrials(Selection(j)).Force.ReSam.L.LatFore{i} +...
%                 DynamicPPTrials(Selection(j)).Force.ReSam.L.MedFore{i} , 'MinPeakProminence',Subject.Newtons * Thresh, 'SortStr','descend');
%             % find whole foot force at the point of max forefoot force
%             DynamicPPTrials(Selection(j)).PropPeak.Left(i).Val = DynamicPPTrials(Selection(j)).Force.ReSam.L.WholeFoot{i}(DynamicPPTrials(Selection(j)).PropPeak.Left(i).Ind);
% %         end
%     end
%
%     % RIGHT Foot
%     for i = 1:DynamicPPTrials(Selection(j)).NumFullStepsR
%         % Peak Force, Pressure, and Area
%         [DynamicPPTrials(Selection(j)).Peak.Right(i).Val, DynamicPPTrials(Selection(j)).Peak.Right(i).Ind] = findpeaks(DynamicPPTrials(Selection(j)).Force.ReSam.R.WholeFoot{i},...
%             'MinPeakProminence',Subject.Newtons* Thresh, 'SortStr','descend');
% %         if strcmp(PPSettings.ToeWalker,'No') == 1
%             %Loading force peak from Heel regions
%             % find instance of loading peak
%             [DynamicPPTrials(Selection(j)).LoadPeak.Right(i).Val, DynamicPPTrials(Selection(j)).LoadPeak.Right(i).Ind] = findpeaks(DynamicPPTrials(Selection(j)).Force.ReSam.R.LatHeel{i} +...
%                 DynamicPPTrials(Selection(j)).Force.ReSam.R.MedHeel{i} , 'MinPeakProminence',Subject.Newtons * Thresh, 'SortStr','descend');
%             % find whole foot force at the point of max heel force
%             DynamicPPTrials(Selection(j)).LoadPeak.Right(i).Val = DynamicPPTrials(Selection(j)).Force.ReSam.R.WholeFoot{i}(DynamicPPTrials(Selection(j)).LoadPeak.Right(i).Ind);
%             % propulsion force peak from forefoot
%             [DynamicPPTrials(Selection(j)).PropPeak.Right(i).Val, DynamicPPTrials(Selection(j)).PropPeak.Right(i).Ind] = findpeaks(DynamicPPTrials(Selection(j)).Force.ReSam.R.LatFore{i} +...
%                 DynamicPPTrials(Selection(j)).Force.ReSam.R.MedFore{i} , 'MinPeakProminence',Subject.Newtons * Thresh, 'SortStr','descend');
%             % find whole foot force at the point of max forefoot force
%             DynamicPPTrials(Selection(j)).PropPeak.Right(i).Val = DynamicPPTrials(Selection(j)).Force.ReSam.R.WholeFoot{i}(DynamicPPTrials(Selection(j)).PropPeak.Right(i).Ind);
% %         end
%     end
% end

%% Computation of overall pressure and force ratios
% peaks are determined by overall force, not pressure
aL = 0; aR = 0;
for  j = 1:length(Selection)
    % LEFT
    % Pressure ratios
    for i = 1:DynamicPPTrials(Selection(j)).NumLeft
        aL = aL + 1;
        % averaged pressure ratios over the entire stance phase
        RatioPress.MedRatio.Left(aL).StepAvg = (sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.L(i).MedHeel)) + sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.L(i).MedArch)) ...
            + sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.L(i).MedFore))) ./ sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.L(i).Whole));
        RatioPress.LatRatio.Left(aL).StepAvg = (sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.L(i).LatHeel)) + sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.L(i).LatArch)) ...
            + sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.L(i).LatFore))) ./ sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.L(i).Whole));
        RatioPress.ArchIndex.Left(aL).StepAvg = (sum(sum(DynamicPPTrials(Selection(j)).MainMask  .* DynamicPPTrials(Selection(j)).Mask.L(i).MedArch)) + sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.L(i).LatArch)))...
            ./ sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.L(i).Whole));
        RatioPress.ArchMedLat.Left(aL).StepAvg = sum(sum(DynamicPPTrials(Selection(j)).MainMask  .* DynamicPPTrials(Selection(j)).Mask.L(i).MedArch)) ./ sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.L(i).LatArch));

        if j == Trial2Cont(1) && i == Lcont % if chosen step values
            RatioPress.MedRatio.Left(1).Selection = (sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.L(i).MedHeel)) + sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.L(i).MedArch)) ...
                + sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.L(i).MedFore))) ./ sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.L(i).Whole));
            RatioPress.LatRatio.Left(1).Selection = (sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.L(i).LatHeel)) + sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.L(i).LatArch)) ...
                + sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.L(i).LatFore))) ./ sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.L(i).Whole));
            RatioPress.ArchIndex.Left(1).Selection = (sum(sum(DynamicPPTrials(Selection(j)).MainMask  .* DynamicPPTrials(Selection(j)).Mask.L(i).MedArch)) + sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.L(i).LatArch)))...
                ./ sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.L(i).Whole));
            RatioPress.ArchMedLat.Left(1).Selection = sum(sum(DynamicPPTrials(Selection(j)).MainMask  .* DynamicPPTrials(Selection(j)).Mask.L(i).MedArch)) ./ sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.L(i).LatArch));
        end
    end
    
    % RIGHT
    % Pressure Ratios
    for i = 1:DynamicPPTrials(Selection(j)).NumRight
        aR = aR + 1;
        % average pressure ratios over the entire stance phase
        RatioPress.MedRatio.Right(aR).StepAvg = (sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.R(i).MedHeel)) + sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.R(i).MedArch)) ...
            + sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.R(i).MedFore))) ./ sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.R(i).Whole));
        RatioPress.LatRatio.Right(aR).StepAvg = (sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.R(i).LatHeel)) + sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.R(i).LatArch)) ...
            + sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.R(i).LatFore))) ./ sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.R(i).Whole));
        RatioPress.ArchIndex.Right(aR).StepAvg = (sum(sum(DynamicPPTrials(Selection(j)).MainMask  .* DynamicPPTrials(Selection(j)).Mask.R(i).MedArch)) + sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.R(i).LatArch)))...
            ./ sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.R(i).Whole));
        RatioPress.ArchMedLat.Right(aR).StepAvg = sum(sum(DynamicPPTrials(Selection(j)).MainMask  .* DynamicPPTrials(Selection(j)).Mask.R(i).MedArch)) ./ sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.R(i).LatArch));
        

        if j == Trial2Cont(2) && i == Rcont % if chosen step values
          RatioPress.MedRatio.Right(1).Selection = (sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.R(i).MedHeel)) + sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.R(i).MedArch)) ...
                + sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.R(i).MedFore))) ./ sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.R(i).Whole));
            RatioPress.LatRatio.Right(1).Selection = (sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.R(i).LatHeel)) + sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.R(i).LatArch)) ...
                + sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.R(i).LatFore))) ./ sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.R(i).Whole));
            RatioPress.ArchIndex.Right(1).Selection = (sum(sum(DynamicPPTrials(Selection(j)).MainMask  .* DynamicPPTrials(Selection(j)).Mask.R(i).MedArch)) + sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.R(i).LatArch)))...
                ./ sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.R(i).Whole));
            RatioPress.ArchMedLat.Right(1).Selection = sum(sum(DynamicPPTrials(Selection(j)).MainMask  .* DynamicPPTrials(Selection(j)).Mask.R(i).MedArch)) ./ sum(sum(DynamicPPTrials(Selection(j)).MainMask .* DynamicPPTrials(Selection(j)).Mask.R(i).LatArch));
    
        end
    end
end

% calcualte averave and standard deviations for step ratios
for i = 1:aL % LEFT SIDE
    % average pressure ratios over the entire stance phase
    RatioPress.MedRatio.LeftAvg.StepAvg = [nanmean([RatioPress.MedRatio.Left.StepAvg]), nanstd([RatioPress.MedRatio.Left.StepAvg])];
    RatioPress.LatRatio.LeftAvg.StepAvg = [nanmean([RatioPress.LatRatio.Left.StepAvg]), nanstd([RatioPress.LatRatio.Left.StepAvg])];
    RatioPress.ArchIndex.LeftAvg.StepAvg = [nanmean([RatioPress.ArchIndex.Left.StepAvg]), nanstd([RatioPress.ArchIndex.Left.StepAvg])];
    RatioPress.ArchMedLat.LeftAvg.StepAvg = [nanmean([RatioPress.ArchMedLat.Left.StepAvg]), nanstd([RatioPress.ArchMedLat.Left.StepAvg])];
%     % average force ratios over the entire stance phase
%     RatioForce.MedRatio.LeftAvg.StepAvg = [nanmean([RatioForce.MedRatio.Left.StepAvg]), nanstd([RatioForce.MedRatio.Left.StepAvg])];
%     RatioForce.LatRatio.LeftAvg.StepAvg = [nanmean([RatioForce.LatRatio.Left.StepAvg]), nanstd([RatioForce.LatRatio.Left.StepAvg])];
%     RatioForce.ArchIndex.LeftAvg.StepAvg = [nanmean([RatioForce.ArchIndex.Left.StepAvg]), nanstd([RatioForce.ArchIndex.Left.StepAvg])];
%     RatioForce.ArchMedLat.LeftAvg.StepAvg = [nanmean([RatioForce.ArchMedLat.Left.StepAvg]), nanstd([RatioForce.ArchMedLat.Left.StepAvg])];
end

for i = 1:aR % RIGHT SIDE
    % average pressure ratios over the entire stance phase
    RatioPress.MedRatio.RightAvg.StepAvg = [nanmean([RatioPress.MedRatio.Right.StepAvg]), nanstd([RatioPress.MedRatio.Right.StepAvg])];
    RatioPress.LatRatio.RightAvg.StepAvg = [nanmean([RatioPress.LatRatio.Right.StepAvg]), nanstd([RatioPress.LatRatio.Right.StepAvg])];
    RatioPress.ArchIndex.RightAvg.StepAvg = [nanmean([RatioPress.ArchIndex.Right.StepAvg]), nanstd([RatioPress.ArchIndex.Right.StepAvg])];
    RatioPress.ArchMedLat.RightAvg.StepAvg = [nanmean([RatioPress.ArchMedLat.Right.StepAvg]), nanstd([RatioPress.ArchMedLat.Right.StepAvg])];
%     % average force ratios over the entire stance phase
%     RatioForce.MedRatio.RightAvg.StepAvg = [nanmean([RatioForce.MedRatio.Right.StepAvg]), nanstd([RatioForce.MedRatio.Right.StepAvg])];
%     RatioForce.LatRatio.RightAvg.StepAvg = [nanmean([RatioForce.LatRatio.Right.StepAvg]), nanstd([RatioForce.LatRatio.Right.StepAvg])];
%     RatioForce.ArchIndex.RightAvg.StepAvg = [nanmean([RatioForce.ArchIndex.Right.StepAvg]), nanstd([RatioForce.ArchIndex.Right.StepAvg])];
%     RatioForce.ArchMedLat.RightAvg.StepAvg = [nanmean([RatioForce.ArchMedLat.Right.StepAvg]), nanstd([RatioForce.ArchMedLat.Right.StepAvg])];
end
clearvars stride RightAdj LeftAdj PkThresh m n MpS i ind j k LogDyn a ans ax Lefties Righties DisplayDim ReEditFPAngles PressForward Thresh Test FLest IPmethods Block.IP HCmethods Block.HC IM

%%  Export and save data
clearvars aL aR DblCheck ForwardRegion i j ListStr NumLeftStrides NumRightStrides SumSteps RGB x txt Bottom
clc; disp('Saving Workspace to MAT file');
% Export report as excel file
% all figures already loaded. this section will take a while due to the multiple sections that have to be written to the excel file.
if strcmp(PPSettings.ExportReport, 'Yes') == 1
    clc; disp('Exporting Report');
    % Page 1
    % subject qualities and date
    Identifiers = {date; folder; Subject.Age; Subject.Height; Subject.Mass};
    xlswrite('ImpressionsOutput.xlsx', Identifiers ,'Sheet1', 'B2:B6');
    
    % Settings Display
    Settings = {PPSettings.MaskChoice; PPSettings.PPMatType; PPSettings.LoadNew; PPSettings.AutoSelectAll; num2str(Trial2Cont)};
    xlswrite('ImpressionsOutput.xlsx', Settings ,'Sheet1', 'M2:M6');
    xlswrite('ImpressionsOutput.xlsx',{Lcont, Rcont} ,'Sheet1', 'M7');
    
    % foot progression angles
    for i = 1:length(Selection)
        % LEFT
        if DynamicPPTrials(Selection(i)).NumLeft > 1
            for j = 1:2
                LProgAng(j,i) = cell2mat(DynamicPPTrials(Selection(i)).LProg.Ang(j));
            end
        elseif DynamicPPTrials(Selection(i)).NumLeft == 0
            LProgAng(1:2,i) = [NaN:NaN];
        else
            LProgAng(1:2,i) = [cell2mat(DynamicPPTrials(Selection(i)).LProg.Ang(1)); NaN];
        end
        % RIGHT
        if DynamicPPTrials(Selection(i)).NumRight > 1
            for j = 1:2
                RProgAng(j,i) = cell2mat(DynamicPPTrials(Selection(i)).RProg.Ang(j));
            end
        elseif DynamicPPTrials(Selection(i)).NumRight == 0
            RProgAng(1:2,i) = [NaN: NaN];
        else
            RProgAng(1:2,i) = [cell2mat(DynamicPPTrials(Selection(i)).RProg.Ang(1)); NaN];
        end
    end
    
    FP_Data = [LProgAng; RProgAng];
    xlswrite('ImpressionsOutput.xlsx', FP_Data,'Sheet1', 'C10');
    
    
    % Page 2
    % Temporal Spatial Data
    if TempSpat == 1
        TSdata = [TS.GaitSpeed.MpS.Avg, TS.GaitSpeed.MpS.Std, TSN.TimeDist(3)/60;...
            TS.GaitSpeed.MpM.Avg, TS.GaitSpeed.MpM.Std, TSN.TimeDist(3);...
            TS.Cadence.Avg, TS.Cadence.Std, TSN.TimeDist(1);...
            TS.StrideLength.Avg, TS.StrideLength.Std, TSN.TimeDist(4);...
            TS.StepWidth.Avg, TS.StepWidth.Std, TSN.TimeDist(6)];
        xlswrite('ImpressionsOutput.xlsx', TSdata  ,'Sheet1', 'Q3:S7');
        % Step Timing
        StepTiming = [TS.Stance.Avg, TS.Stance.Std, .62;...
            TS.Swing.Avg, TS.Swing.Std, .38;...
            TS.IDS.Avg, TS.IDS.Std, .12;...
            TS.SS.Avg, TS.SS.Std, .38;...
            TS.SDS.Avg, TS.SDS.Std, .12];
        xlswrite('ImpressionsOutput.xlsx', 100*StepTiming  ,'Sheet1', 'Q10:S14');
    end
    % Page 3
    % Pressure Ratios
    PressRatios = [RatioPress.MedRatio.Left.Selection, RatioPress.MedRatio.Right.Selection,...
        RatioPress.MedRatio.LeftAvg.StepAvg, RatioPress.MedRatio.RightAvg.StepAvg;...
        RatioPress.LatRatio.Left.Selection, RatioPress.LatRatio.Right.Selection,...
        RatioPress.LatRatio.LeftAvg.StepAvg, RatioPress.LatRatio.RightAvg.StepAvg;...
        RatioPress.ArchIndex.Left.Selection, RatioPress.ArchIndex.Right.Selection,...
        RatioPress.ArchIndex.LeftAvg.StepAvg, RatioPress.ArchIndex.RightAvg.StepAvg;...
        RatioPress.ArchMedLat.Left.Selection, RatioPress.ArchMedLat.Right.Selection,...
        RatioPress.ArchMedLat.LeftAvg.StepAvg, RatioPress.ArchMedLat.RightAvg.StepAvg];
    xlswrite('ImpressionsOutput.xlsx', PressRatios  ,'Sheet1', 'V5');
    
    %    CoP Metrics
    aL = 0; aR = 0;
    for i = 1:length(Regions)
        for j = 1:length(Regions{i})
            if strcmp(Regions{i}(j).Side, 'Left')
                aL = aL + 1;
                CoPIndex.Left(aL).Start = Regions{i}(j).StepCoP(1,2) / Regions{i}(j).BoundingBox(4);
                CoPIndex.Left(aL).End = Regions{i}(j).StepCoP(end,2) / Regions{i}(j).BoundingBox(4);
            else
                aR = aR + 1;
                CoPIndex.Right(aR).Start = Regions{i}(j).StepCoP(1,2) / Regions{i}(j).BoundingBox(4);
                CoPIndex.Right(aR).End = Regions{i}(j).StepCoP(end,2) / Regions{i}(j).BoundingBox(4);
            end
        end
    end
    
    CoPData = [Regions{Trial2Cont(1)}(LcontReg).StepCoP(1,2) / Regions{Trial2Cont(1)}(LcontReg).BoundingBox(4), Regions{Trial2Cont(2)}(RcontReg).StepCoP(1,2) / Regions{Trial2Cont(2)}(RcontReg).BoundingBox(4), ...
        nanmean([CoPIndex.Left.Start]), nanstd([CoPIndex.Left.Start]), nanmean([CoPIndex.Right.Start]), nanstd([CoPIndex.Right.Start]);...
        Regions{Trial2Cont(1)}(LcontReg).StepCoP(end,2) / Regions{Trial2Cont(1)}(LcontReg).BoundingBox(4), Regions{Trial2Cont(2)}(RcontReg).StepCoP(end,2) / Regions{Trial2Cont(2)}(RcontReg).BoundingBox(4), ...
        nanmean([CoPIndex.Left.End]), nanstd([CoPIndex.Left.End]), nanmean([CoPIndex.Right.End]), nanstd([CoPIndex.Right.End])];
    
    % cop index metrics?
    xlswrite('ImpressionsOutput.xlsx', CoPData ,'Sheet1', 'V10:AA11');
    
    
    
    % Page 4
    
    % data to display here?
    
    
    % export as PDF
    XL_path = GetFullPath('ImpressionsOutput.xlsx');
    NewName = strsplit(XL_path, '\');
    NewName(end) = [];
    FolderPath = strcat(strjoin(NewName, '\'),'\', folder);
    PDF_path = strcat(FolderPath,'\', 'ImpressionsOutput.pdf');
    xls2pdf(XL_path, 'Sheet1',PDF_path, 'A1: AQ46');
%     delete ImpressionsOutput.xlsx
end
clearvars NewName XL_path

%% Validation Plots
% if strcmp(PPPlots.Validation, 'Yes')
%     if strcmp(PPSettings.PPMatType, 'Novel')
%         for i = 1:NumTrials
%             [~,~,p] = size(DynamicPPTrials(i).TM);
%             for j = 1:p
%                 Overall(i).PeakPressure(j) = max(max(DynamicPPTrials(Selection(i)).TM(:,:,j))); % peak pressure in kPa
%                 Overall(i).Pressure(j) = sum(sum(DynamicPPTrials(Selection(i)).TM(:,:,j))); % peak pressure in kPa
%                 A = DynamicPPTrials(Selection(i)).TM(:,:,j);
%                 A(A==0) = NaN;
%                 Overall(i).MeanPressure(j) = nanmean(nanmean(A)); % peak pressure in kPa
%                 Overall(i).Area(j) = sum(sum(DynamicPPTrials(Selection(i)).TMLog(:,:,j))); % area in cells
%                 Overall(i).Force(j) =  Overall(i).Pressure(j) .* Overall(i).Area(j) / 10000; % force in N
%             end
%             Overall(i).Area = Overall(i).Area  .* Adj.Area; % area in cm^2
%             % filter force data
%             windowSize = 3;
%             b = (1/windowSize)*ones(1,windowSize);
%             a = 1;
%             Overall(i).ForceFilt = filter(b,a,Overall(i).Force);
%             figure;
%             subplot(311);
%             plot(Overall(i).PeakPressure); grid on; hold on;
%             plot(Overall(i).MeanPressure, 'r');
%             legend('Peak Pressure','Mean Pressure');
%             title(AnalyzePPImages(i).file_name);
%             ylabel('Pressure (kPa)');
%             xlabel('Frame');
%             subplot(312);
%             plot(Overall(i).Force, 'r'); grid on; hold on;
%             plot(Overall(i).ForceFilt);
%             legend('Raw Force','Filtered Force');
%             ylabel('Force (N)');
%             xlabel('Frame');
%             subplot(313);
%             plot(Overall(i).Area);  grid on;
%             ylabel('Area (cm^2)');
%             xlabel('Frame');
%             
%             f = char(AnalyzePPImages(i).file_name);
%             savefig(strcat(f(1:end-4), '_Validation.fig'));
%             
%         end
%     end
% end
% 
% 
% %% create a movie of each trial
% % close all;
% if strcmp(PPSettings.MovieOutput,'Yes') == 1
%     clc; disp('Saving GIFs');
%     %     addpath(FolderPath);
%     addpath(folder)
%     if strcmp(PPSettings.PPMatType, 'Novel')
%         SampFreq = 100;
%     elseif strcmp(PPSettings.PPMatType, 'RSScan')
%         SampFreq = 126;
%     end
%     ReplaySpeed = 1/5;
%     if Trial2Cont(1) == Trial2Cont(2)
%         z = Trial2Cont(1);
%         if strcmp(PPSettings.MovieTypes, 'Full Trial') || strcmp(PPSettings.MovieTypes, 'Both')
%             figure('Position',[100 100 200 600]);
%             axis equal;
%             contour(PPTrials(Selection(z)).TM(:,:,1), 25);
%             text(3,5,strcat(num2str(round(ReplaySpeed*100)),'% speed'), 'FontSize', 6);
%             hold on;
%             plot(PPTrials(Selection(z)).CoP(1,1), PPTrials(Selection(z)).CoP(1,2), '.k', 'MarkerSize',14);
%             gif(strcat(FolderPath,'\','Trial.gif'),'DelayTime',1/(SampFreq*ReplaySpeed));
%             hold off;
%             ax = gca;
%             ax.XTick = [];
%             ax.YTick = [];
%             for i = 2:size(PPTrials(Selection(z)).TM,3)
%                 contour(PPTrials(Selection(z)).TM(:,:,i), 25);
%                 hold on;
%                 text(3,5,strcat(num2str(round(ReplaySpeed*100)),'% speed'), 'FontSize', 6);
%                 ax.XTick = [];
%                 ax.YTick = [];
%                 plot(PPTrials(Selection(z)).CoP(i-1,1), PPTrials(Selection(z)).CoP(i-1,2), '.r', 'MarkerSize',10);
%                 if i > 3
%                     plot(PPTrials(Selection(z)).CoP(i-2,1), PPTrials(Selection(z)).CoP(i-2,2), '.r', 'MarkerSize',10);
%                     plot(PPTrials(Selection(z)).CoP(i-3,1), PPTrials(Selection(z)).CoP(i-3,2), '.r', 'MarkerSize',10);
%                 end
%                 plot(PPTrials(Selection(z)).CoP(i,1), PPTrials(Selection(z)).CoP(i,2), '.k', 'MarkerSize',14);
%                 hold off;
%                 gif;
%             end
%             clf;
%             close;
%         end
%     else
%         uiwait(msgbox('Two separate trials were used for close ups, therefore the full trial cannot be exported'));
%     end
%     
%     % plot foot contour and CoP
%     if strcmp(PPSettings.MovieTypes, 'Single Steps') || strcmp(PPSettings.MovieTypes, 'Both')
%         % combine L and R steps into one matrix
%         [Lm, Ln, Lp] = size(Regions{Selection(z)}(LcontReg).Step(:,:,:));
%         [Rm, Rn, Rp] = size(Regions{Selection(z)}(RcontReg).Step(:,:,:));
%         
%         L_Step = Regions{Selection(z)}(LcontReg).Step(:,:,:);
%         R_Step = Regions{Selection(z)}(RcontReg).Step(:,:,:);
%         % align height
%         if Lm > Rm
%             d = Lm - Rm;
%             M = Lm;
% %             clearvars R_Step
% %             for i = 1:Rp
% %                 R_Step(:,:,i) = vertcat(Regions{Selection(z)}(RcontReg).Step(:,:,i), zeros(d,Rn));
% %             end
%             R_Step(1:Rm+d,:,:) = vertcat(R_Step(:,:,:), zeros(d,Rn,Rp));
%             Rm = Rm + d; 
%         elseif Rm > Lm
%             d = Rm - Lm;
%             M = Rm;
% %             clearvars L_Step
% %             for i = 1:Lp
% %                 L_Step(:,:,i) = vertcat(Regions{Selection(z)}(LcontReg).Step(:,:,i), zeros(d,Ln));
% %             end
%             L_Step(1:Lm+d,:,:) = vertcat(L_Step(:,:,:), zeros(d,Ln,Lp));
%             Lm = Lm + d;
%         elseif Rm == Lm
%             M = Rm;
%         end
%         % align width
%         if Ln > Rn
%             d = Ln - Rn;
%             N = Ln;
%             R_Step(:,1:Rn+d,:) = horzcat(R_Step(:,:,:), zeros(Rm,d,Rp));
%             Rn = Rn + d; 
%         elseif Rn > Ln
%             N = Rn;
%             d = Rn - Ln;
%             L_Step(:,1:Ln+d,:) = horzcat(L_Step(:,:,:), zeros(Lm,d,Lp));
%             Ln = Ln + d; 
%         elseif Rn == Ln
%             N = Rn;
%         end
%         
%         % align length (time)
%         if Lp > Rp
% %             d = Lp - Rp;
%             P = Lp; 
%             for i = Rp:Lp
%                 R_Step(:,:,i) = zeros(M,N);
%             end
%         elseif Rp > Lp
%             P = Rp; 
% %             d = Rp - Lp;
%             for i = Lp:Rp
%                 L_Step(:,:,i) = zeros(M,N);
%             end
%         end
%         
%         STEP = horzcat(L_Step, R_Step);
%         
%         figure('Position',[100 100 800 800]);
%         contour(STEP(:,:,1), 50); axis equal; hold on;
%         vline(Ln-0.5, 'k'); 
%         text(1.5,1.5,strcat(num2str(round(ReplaySpeed*100)),'% speed'), 'FontSize', 6);
%         plot(Regions{Selection(z)}(LcontReg).StepCoP(1,1), Regions{Selection(z)}(LcontReg).StepCoP(1,2), '.k', 'MarkerSize',14); % left CoP
%         plot(N + Regions{Selection(z)}(RcontReg).StepCoP(1,1), Regions{Selection(z)}(RcontReg).StepCoP(1,2), '.k', 'MarkerSize',14); % right CoP
%         axis equal;
%         ax = gca;
%         ax.XTick = [];
%         ax.YTick = [];
%         gif(strcat(FolderPath,'\', 'L-RSteps.gif'),'DelayTime',1/(SampFreq*ReplaySpeed));
%         for i = 2:P
%             contour(STEP(:,:,i), 50);
%             hold on;
%             vline(Ln-0.5, 'k');
%             text(1.5,1.5,strcat(num2str(round(ReplaySpeed*100)),'% speed'), 'FontSize', 6);
%             ax.XTick = [];
%             ax.YTick = [];
%             if i <= Lp
%                 plot(Regions{Selection(z)}(LcontReg).StepCoP(i-1,1), Regions{Selection(z)}(LcontReg).StepCoP(i-1,2), '.r', 'MarkerSize',10); % left CoP
%                 if i > 3
%                     plot(Regions{Selection(z)}(LcontReg).StepCoP(i-2,1), Regions{Selection(z)}(LcontReg).StepCoP(i-2,2), '.r', 'MarkerSize',10); % left CoP
%                     plot(Regions{Selection(z)}(LcontReg).StepCoP(i-3,1), Regions{Selection(z)}(LcontReg).StepCoP(i-3,2), '.r', 'MarkerSize',10); % left CoP
%                 end
%                 plot(Regions{Selection(z)}(LcontReg).StepCoP(i,1), Regions{Selection(z)}(LcontReg).StepCoP(i,2), '.k', 'MarkerSize',14); % left CoP
%             end
%             if i <= Rp
%                 plot(N + Regions{Selection(z)}(RcontReg).StepCoP(i-1,1), Regions{Selection(z)}(RcontReg).StepCoP(i-1,2), '.r', 'MarkerSize',10); % right CoP
%                 if i > 3
%                     plot(N + Regions{Selection(z)}(RcontReg).StepCoP(i-2,1), Regions{Selection(z)}(RcontReg).StepCoP(i-2,2), '.r', 'MarkerSize',10); % right CoP
%                     plot(N + Regions{Selection(z)}(RcontReg).StepCoP(i-3,1), Regions{Selection(z)}(RcontReg).StepCoP(i-3,2), '.r', 'MarkerSize',10); % right CoP
%                 end
%                 plot(N + Regions{Selection(z)}(RcontReg).StepCoP(i,1), Regions{Selection(z)}(RcontReg).StepCoP(i,2), '.k', 'MarkerSize',14); % right CoP
%             end
%             hold off; axis equal;
%             gif;
%         end
%         clf;
%         close;
%         
%         
%     end
% 
% end

%% Export PP Data as .mat file and finish
clearvars aL aR ans i j L R z filename Identifiers NumSteps indX indY IND Empty AddX AddY TStrials W y
save(strcat(char(folder),'\', char(folder),'_Data'));
close all;
clc; disp('Impressions has finished analyzing all plantar pressures. See Excel output, matlab data, figures, and GIFs for results.');




