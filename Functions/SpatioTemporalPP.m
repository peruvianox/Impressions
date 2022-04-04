function [Regions, ToeDrag] = SpatioTemporalPP(TM, PixelThresh, PPSettings, PPDim, CoPLine)

% Using a time series of plantar pressure matricies, this function will
% determine the temporal spatial parameters for a subject who walks across
% the plantar pressure mat. It will also give the bounding boxes

TrialType = PPSettings.PPMatType; 

%% Image processing of trial
SumTM = sum(TM, 3); % create summed image
CC = bwconncomp(SumTM,8);
% delete groups that are too small
for i = 1:CC.NumObjects
    [M,~] = size(CC.PixelIdxList{i});
    if M < 10
        CC.PixelIdxList{i} = [];
    end
end

%% create logical array of connected areas
BW = labelmatrix(CC);
% figure;
% contour(BW,100);
% axis equal;
Reg = regionprops(BW,'Area','Centroid','BoundingBox');

%% find distances between the clusters and merge if they are within the threshold
L = length(Reg);
% PixelThresh set to estimated foot length
D = zeros(L,1);
Pair = zeros(L,2);
ToMerge = zeros(L,2);
n =  1;
for i = 1:L-1
    for j = i+1:L
        Pair(n,:)  = [i,j];
        D(n) = sqrt((abs(Reg(i).Centroid(1) -  Reg(j).Centroid(1)))^2 + (abs(Reg(i).Centroid(2) -  Reg(j).Centroid(2)))^2);
        if D(n) < PixelThresh
            ToMerge(n,:) = Pair(n,:);
        else
            ToMerge(n,:) = [0,0];
        end
        n = n+1;
    end
end

% A = find(ToMerge(:,1));
Merge = ToMerge(ToMerge(:,1) > 0,:);

%% Merge regions if needed
BWnew = BW; % copy to a new variable
% find any intersection points, where there are > 2 regions combined together
[IntersectReg, ia, ib] = intersect(Merge(:,1),Merge(:,2));
% if there are intersected regions, combine them
if isempty(IntersectReg) == 0
    for i = 1:length(IntersectReg)
        % if there are > 3 regions that need to be joined, then
        % Join first 2 regions
        BWlog = logical(BW == Merge(ia(i),2)); % find number(s) to merge
        % sum(sum(BWlog)) % uncomment to check area of altered region
        BWnew(BWlog == 1) = Merge(ia(i),1);
        % Join second 2 regions
        BWlog = logical(BW == Merge(ib(i),1));
        % sum(sum(BWlog)) % uncomment to check area of altered region
        BWnew(BWlog == 1) = Merge(ib(i),2);
    end
end
% identify rows to delete to clean out intersected regions
[S,~] = size(Merge);
MergeLog = zeros(S, length(IntersectReg));
for j = 1:S
    for i = 1:length(IntersectReg)
        if Merge(j,1) == Merge(ib(i),1) || Merge(j,2) == Merge(ib(i),1)
            MergeLog(j, i) = 1;
        elseif Merge(j,1) == Merge(ia(i),2) || Merge(j,2) == Merge(ia(i),2)
            MergeLog(j, i) = 1;
        elseif Merge(j,1) == IntersectReg(i) || Merge(j,2) == IntersectReg(i)
            MergeLog(j, i) = 1;
        else
            MergeLog(j, i) = 0;
        end
    end
end
% delete intersected rows of Merge
MergeDel = sum(MergeLog, 2);
Merge(MergeDel > 0,:) = [];

% Delete non-intersected regions, when only 2 are being combined
[S,~] = size(Merge);
if S >= 1
    for i = 1:S
        % create logical of column 2
        BWlog = logical(BW == Merge(i,2));
        sum(sum(BWlog)); % uncomment to check area of altered region
        BWnew(BWlog == 1) = Merge(i,1);
    end
end

clearvars Merge MergeDel MergeLog S 

RegNew = regionprops(BWnew,'Area','Centroid','BoundingBox','MajorAxisLength','MinorAxisLength','Orientation','Eccentricity','ConvexHull');
Reg = RegNew;

%% search for toe drag and delete if present
% if strcmp(PPSettings.ToeWalker,'No')
A = 0;
for i = 1:length(Reg)
    if  Reg(i).BoundingBox(4) > PixelThresh * 1.33 % toe drag defined as an anterior-posterior bounding box of 1.33 times greater than the estimated foot length
        A = A+1;
        uiwait(msgbox(['Impressions has detected foot drag within a trial. Select drag area(s) to ignore.']));
        RemoveToeDrag = 'Yes';
        j = 1;
        if A == 1
            figure( 'Position', [100, 100, 300, 500]);
        end
        while strcmp(RemoveToeDrag, 'Yes') == 1
            contour(BWnew, 100);
            h = imrect;
            Delete(j).Coord = round(getPosition(h));
            if Delete(j).Coord(1) < 1
                Delete(j).Coord(1) = 1;
            end
            if Delete(j).Coord(2) < 1
                Delete(j).Coord(2) = 1;
            end
            BWnew(Delete(j).Coord(2):Delete(j).Coord(2)+Delete(j).Coord(4), Delete(j).Coord(1):Delete(j).Coord(1)+Delete(j).Coord(3)) = 0;
            Question = 'Would you like to remove more toe drags?';
            RemoveToeDrag = questdlg(Question,'Remove Toe Drag','Yes','No','No');
            ToeDrag(i).Delete(j).Coord = Delete(j).Coord;
            clearvars h
            clf; j = j + 1;
        end
    else
        ToeDrag(i).Delete = 0;
    end
end
if A > 0
    RegNew = regionprops(BWnew,'Area','Centroid','BoundingBox','MajorAxisLength','MinorAxisLength','Orientation','Eccentricity','ConvexHull');
    Reg = RegNew;
end
close;

%% Double check about merged regions before continuing
BW = BWnew;
L = length(Reg);
D = zeros(L,1);
Pair = zeros(L,2);
ToMerge = zeros(L,2);
n =  1;
for i = 1:L-1
    for j = i+1:L
        Pair(n,:)  = [i,j];
        D(n) = sqrt((abs(Reg(i).Centroid(1) -  Reg(j).Centroid(1)))^2 + (abs(Reg(i).Centroid(2) -  Reg(j).Centroid(2)))^2);
        if D(n) < PixelThresh
            ToMerge(n,:) = Pair(n,:);
        else
            ToMerge(n,:) = [0,0];
        end
        n = n+1;
    end
end

Merge = ToMerge(ToMerge(:,1) > 0,:);

%% Merge regions if needed
% find any intersection points, where there are > 2 regions combined together
[IntersectReg, ia, ib] = intersect(Merge(:,1),Merge(:,2));
% if there are intersected regions, combine them
if isempty(IntersectReg) == 0
    for i = 1:length(IntersectReg)
        % if there are > 3 regions that need to be joined, then
        % Join first 2 regions
        BWlog = logical(BW == Merge(ia(i),2)); % find number(s) to merge
        % sum(sum(BWlog)) % uncomment to check area of altered region
        BWnew(BWlog == 1) = Merge(ia(i),1);
        % Join second 2 regions
        BWlog = logical(BW== Merge(ib(i),1));
        % sum(sum(BWlog)) % uncomment to check area of altered region
        BWnew(BWlog == 1) = Merge(ib(i),2);
    end
    % identify rows to delete to clean out intersected regions
    [S,~] = size(Merge);
    MergeLog = zeros(S, length(IntersectReg));
    for j = 1:S
        for i = 1:length(IntersectReg)
            if Merge(j,1) == Merge(ib(i),1) || Merge(j,2) == Merge(ib(i),1)
                MergeLog(j, i) = 1;
            elseif Merge(j,1) == Merge(ia(i),2) || Merge(j,2) == Merge(ia(i),2)
                MergeLog(j, i) = 1;
            elseif Merge(j,1) == IntersectReg(i) || Merge(j,2) == IntersectReg(i)
                MergeLog(j, i) = 1;
            else
                MergeLog(j, i) = 0;
            end
        end
    end
    % delete intersected rows of Merge
    MergeDel = sum(MergeLog, 2);
    Merge(MergeDel > 0,:) = [];
end

%% Delete non-intersected regions, when only 2 are being combined
[S,~] = size(Merge);
if S >= 1
    for i = 1:S
        % create logical of column 2
        BWlog = logical(BWnew == Merge(i,2));
        sum(sum(BWlog)); % uncomment to check area of altered region
        BWnew(BWlog == 1) = Merge(i,1);
    end
end

clearvars Merge MergeDel MergeLog S 
RegNew = regionprops(BWnew,'Area','Centroid','BoundingBox','MajorAxisLength','MinorAxisLength','Orientation','Eccentricity','ConvexHull');

%% delete merged areas
for i = 1:length(RegNew)
    if RegNew(i).Area > 0
        ToDel(i) = 0;
    else
        ToDel(i) = 1; 
    end
end
Reg = RegNew(~ToDel); 
clearvars IntersectReg ia ib Pair D ToMerge Merge

%% Check PPs if no AutoSelect
if strcmp(PPSettings.AutoSelectAll, 'No')
    % plot for visual feedback
        figure( 'Position', [100, 100, 300, 500]);
    contour(BWnew, 100);
    hold on;
    for i = 1:length(RegNew)
        rectangle('Position', RegNew(i).BoundingBox, 'EdgeColor','k','LineWidth',2);
    end
    axis equal;
    RedoRegions = questdlg('Does the pressure grouping need to be redone?','Redo Regions','Yes','No','No');
    if strcmp(RedoRegions,'Yes') == 1
        clearvars IntersectReg ia ib Pair D ToMerge Merge
        % reload regions
        %     Reg = regionprops(BWnew,'Area','Centroid','BoundingBox');
        % find distances between the clusters and merge if they are within the threshold
        L = length(Reg);
        PixelThresh = PixelThresh * 1.2;
        % Loop for measuring distances between pressures
        D = zeros(L,1);
        Pair = zeros(L,2);
        ToMerge = zeros(6,2);
        n =  1;
        for i = 1:L-1
            for j = i+1:L
                Pair(n,:)  = [i,j];
                D(n) = sqrt((abs(Reg(i).Centroid(1) -  Reg(j).Centroid(1)))^2 + (abs(Reg(i).Centroid(2) -  Reg(j).Centroid(2)))^2);
                if D(n) < PixelThresh
                    ToMerge(n,:) = Pair(n,:);
                else
                    ToMerge(n,:) = [0,0];
                end
                n = n+1;
            end
        end
        
        Merge = ToMerge(ToMerge(:,1) > 0,:);
        % Merge regions if needed
        BWnew2 = BWnew;
        % find any intersection points, where there are > 2 regions combined together
        [IntersectReg, ia, ib] = intersect(Merge(:,1),Merge(:,2));
        % if there are intersected regions, combine them
        if isempty(IntersectReg) == 0
            for i = 1:length(IntersectReg)
                % if there are > 3 regions that need to be joined, then
                % Join first 2 regions
                BWlog = logical(BWnew == Merge(ia(i),2));
                % sum(sum(BWlog)) % uncomment to check area of altered region
                BWnew2(BWlog == 1) = Merge(ia(i),1);
                % Join second 2 regions
                BWlog = logical(BWnew == Merge(ib(i),1));
                % sum(sum(BWlog)) % uncomment to check area of altered region
                BWnew2(BWlog == 1) = Merge(ib(i),2);
            end
            % delete intersected rows of Merge, to clean out intersected regions
            [S,~] = size(Merge);
            MergeLog = zeros(S, length(IntersectReg));
            for j = 1:S
                for i = 1:length(IntersectReg)
                    if Merge(j,1) == Merge(ib(i),1) || Merge(j,2) == Merge(ib(i),1)
                        MergeLog(j, i) = 1;
                    elseif Merge(j,1) == Merge(ia(i),2) || Merge(j,2) == Merge(ia(i),2)
                        MergeLog(j, i) = 1;
                    elseif Merge(j,1) == IntersectReg(i) || Merge(j,2) == IntersectReg(i)
                        MergeLog(j, i) = 1;
                    else
                        MergeLog(j, i) = 0;
                    end
                end
            end
            MergeDel = sum(MergeLog, 2);
            Merge(MergeDel > 0,:) = [];
        end
        
        [S,~] = size(Merge);
        if S >= 1
            for i = 1:S
                % create logical of column 2
                BWlog = logical(BW == Merge(i,2));
                sum(sum(BWlog)); 
                BWnew(BWlog == 1) = Merge(i,1);
            end
        end
        
        RegNew2 = regionprops(BWnew2,'Area','Centroid','BoundingBox','MajorAxisLength','MinorAxisLength','Orientation','Eccentricity','ConvexHull');
        Reg = RegNew2;
        
        % verify merging
        clf;
        contour(BWnew2, 100);
        hold on;
        for i = 1:length(RegNew2)
            rectangle('Position', RegNew2(i).BoundingBox, 'EdgeColor','k','LineWidth',2);
        end
        axis equal;
        
        %% if pressures still need correcting, manually select foot pressures
        RedoRegionsAg = questdlg('Does the pressure grouping need to be redone, again?','Redo Regions Again');
        % if is still doesn't work, treat each pressure as an individual foot
        if strcmp(RedoRegionsAg, 'Yes') == 1
            clf;
            %     figure( 'Position', [100, 100, 300, 500]);
            %     contour(BWnew, 100);
            %     hold on;
            CC = bwconncomp(SumTM,8);
            % delete groups that are too small
            for i = 1:CC.NumObjects
                [M,~] = size(CC.PixelIdxList{i});
                if M < 10
                    CC.PixelIdxList{i} = [];
                end
            end
            % create logical array of connected areas
            BW = labelmatrix(CC);
            %     figure( 'Position', [100, 100, 300, 500]);
            hold on;
            contour(BW,100);
            axis equal;
            
            ListStr = {'1','2','3','4','5','6'};
            NumPress = listdlg( 'PromptString','How many pressires in trial?','ListString',ListStr); % determine which trials to re-edit
            uiwait(msgbox('Using your cursor, draw a box around the foot to identify.', 'Instructions'));
            
            BWlog = logical(BW > 0); % create logical matrix
            
            for i = 1:NumPress % loop through number of foot pressures
                [m,n] = size(BWlog);
                BWclass = zeros(m,n); % create empty matrix
                
                h = imrect; % identify foot pressures
                Comb(i).Coord = round(getPosition(h));
                
                bot = Comb(i).Coord(2); % identify bounding box
                top = Comb(i).Coord(2) + Comb(i).Coord(4);
                L = Comb(i).Coord(1);
                R = Comb(i).Coord(1) + Comb(i).Coord(3);
                if bot < 1
                    bot = 1;
                end
                if top > PPDim.Length
                    top = PPDim.Length;
                end
                if L < 1
                    L = 1;
                end
                if R > PPDim.Width
                    R = PPDim.Width;
                end
                
                BWclass(bot:top, L:R) = 1;
                BWadd = BWlog + BWclass; % add matricies todether
                BWnew = BWadd == 2; % identify region for foot pressure
                BW_new(i).mat = i*BWnew; % multiply by loop number to classify
                
                clearvars bot top L R h BWclass
            end
            
            if length(BW_new) == 1 % add together all pressures
                BWNEW = BW_new(1).mat;
            elseif  length(BW_new) == 2
                BWNEW = BW_new(1).mat + BW_new(2).mat;
            elseif  length(BW_new) == 3
                BWNEW = BW_new(1).mat + BW_new(2).mat + BW_new(3).mat;
            elseif  length(BW_new) == 4
                BWNEW = BW_new(1).mat + BW_new(2).mat + BW_new(3).mat + BW_new(4).mat;
            elseif  length(BW_new) == 5
                BWNEW = BW_new(1).mat + BW_new(2).mat + BW_new(3).mat + BW_new(4).mat + BW_new(5).mat;
            end
            
%                  figure( 'Position', [100, 100, 300, 500]);
%                  hold on;
%                  contour(BWNEW,100);
            
            Reg = regionprops(BWNEW,'Area','Centroid','BoundingBox','MajorAxisLength','MinorAxisLength','Orientation','Eccentricity','ConvexHull');
            axis equal;
            % loop to plot bounding boxes identifying feet
            for i = 1:length(Reg)
                rectangle('Position', Reg(i).BoundingBox, 'EdgeColor','k','LineWidth',2);
            end
        end
    end
    close;
end

%% once regions are all set, delete the rows with too small areas, or that are too close to ends of mat
% figure( 'Position', [100, 100, 300, 500]);
% hold on;
% contour(BW,100);
Z = zeros(1,length(Reg)); 
for i = 1:length(Reg)
    if Reg(i).Area < 30 % assigns to delete regions that are too small to be feet (area < 30)
        Z(i) = 1;
    end
    if Reg(i).BoundingBox(1) < 1
        Z(i) = 1; 
    elseif Reg(i).BoundingBox(2) < 1
        Z(i) = 1;
    elseif Reg(i).BoundingBox(1) + Reg(i).BoundingBox(3) > PPDim.Width - 1
        Z(i) = 1; 
    elseif Reg(i).BoundingBox(2) + Reg(i).BoundingBox(4) > PPDim.Length - 1
        Z(i) = 1;
    end
end
% Delete empty rows of the regions
Reg(Z==1) = [];
clearvars D Z ToDel

%% Find foot strikes and foot offs for each step
[L,W,P] = size(TM);
Zones = length(Reg);
k = 1; 
for j = 1:Zones % define # of regions to be analyzed
    % define search zone based on region
    Search = round([Reg(j).BoundingBox(1), Reg(j).BoundingBox(2), Reg(j).BoundingBox(1) + Reg(j).BoundingBox(3), Reg(j).BoundingBox(2) + Reg(j).BoundingBox(4)]);
    % Correct any of the search regions if they are out of the matrix region
    if Search(1) < 0
        Search(1) = 0;
    end
    if Search(2) < 0
        Search(2) = 0;
    end
    if Search(3) > W
        Search(3) = W;
    end
    %     OffTime = 0;
    if Search(4) > L
        Search(4) = L;
        %         OffTime = 1;
    end
    
    SUM = zeros(1,P); %pre-allocate sum region
    for i = 1:P % find the sum of the region to determine if a step is present
        SUM(i) = sum(sum(TM(Search(2):Search(4), Search(1):Search(3),i)));
    end
    % find the foot strike and foot off of the step
    Reg(j).Strike  = find(SUM,1,'first');
    Reg(j).Off = find(SUM,1,'last');
    if strcmp(PPSettings.PPMatType, 'RSScan')
        if Reg(j).Off == 248 &&  Reg(j).Strike > 200 % if foot off is defines as last time-matrix -> delete to ensure data accuracy
            ToDel(k) = j; 
            k = k +1; 
        end
    end
end
if exist('ToDel', 'var')
    Reg(ToDel) = [];
end

%% Re-define step characteristics in seconds
if strcmp(TrialType, 'Novel') == 1
    TimeAdj = 100; % Novel measures at 100 Hz
elseif  strcmp(TrialType, 'RSScan') == 1
    TimeAdj = 126; % RSscan measures at 126 Hz
end

for i = 1:length(Reg)
    Reg(i).TimeStrike = Reg(i).Strike / TimeAdj;
    Reg(i).TimeOff = Reg(i).Off / TimeAdj;
end

%% Re-order output based on event times
% RegSorted = sort([Reg(:).Strike], 'ascend');
RegNames = fieldnames(Reg);
RegCell  = struct2cell(Reg);

[A,B,p3] = size(RegCell);
for i = 1:A
    for j = 1:B
        if isempty(RegCell{i,j})
            RegCell{i,j} = 0;
        end
    end
end
if p3 == 1
    RegSorted = sortrows(RegCell', 7);
else
    RegSorted = sortrows(permute(RegCell,[3 1 2]), 7);
end

RegStruct = cell2struct(RegSorted', RegNames, 1);
Reg = RegStruct; 

%% Label Right and left feet
% find centroids for all feet
for i = 1:length(Reg)
    Locs(i) = Reg(i).Centroid(1);
end
 MeanCentLoc = mean(Locs);
 
if length(Reg) > 2 || length(Reg) == 1 % if  >2 or only 1 full step(s) on mat 
    % create linear trendline down long axis of mat
%       CoPLine = DynamicPPTrials(3).CoP;
%       BW = DynamicPPTrials(3).SumTM; 
    pCoP = polyfit(CoPLine(:,2), CoPLine(:,1),1); % create polynomial fit to the 1st degree = trendline
    x = linspace(1,PPDim.Length,PPDim.Length);
    y = polyval(pCoP, x);
%         figure( 'Position', [100, 100, 300, 500]);
%         hold on;
%         contour(BW,100);
%         plot(y,x);
    for i = 1:length(Reg)
        % find average position of all centroids and use that location to determine L and R feet
        if Reg(i).Centroid(1) < y(round(Reg(i).Centroid(2))) % if horiz position of centroid is less than trendline -> LEFT
            Reg(i).Side = 'Left';
            Reg(i).SideErr = abs(Reg(i).Centroid(1) - y(round(Reg(i).Centroid(2)))); 
        else
            Reg(i).Side = 'Right'; % if horiz position of centroid is greater than trendline -> RIGHT
             Reg(i).SideErr = abs(Reg(i).Centroid(1) - y(round(Reg(i).Centroid(2)))); 
        end
    end
elseif length(Reg) == 2 % if only 2 pressures, use centroid locations to label
    for i = 1:length(Reg)
        % find average position of all centroids and use that location to determine L and R feet
        if Reg(i).Centroid(1) < MeanCentLoc % if horiz position of centroid is less than trendline -> LEFT
            Reg(i).Side = 'Left';
             Reg(i).SideErr = abs(Reg(i).Centroid(1) - MeanCentLoc); 
        else
            Reg(i).Side = 'Right'; % if horiz position of centroid is greater than trendline -> RIGHT
             Reg(i).SideErr = abs(Reg(i).Centroid(1) - MeanCentLoc); 
        end
    end
end

%% reoder steps in order of occurrence
for i =1:length(Reg)
    Strikes(i) = Reg(i).Strike;
end
if exist('Strikes', 'var')
    [~,Order] = sort(Strikes);
    Reg = Reg(Order);
end

%% double check classification
for i = 1:length(Reg)-1
    if strcmp(Reg(i).Side, 'Left')
        if strcmp(Reg(i+1).Side, 'Right')
            AltErr(i) = 0; 
        else
            AltErr(i) = 1;
        end
    else
          if strcmp(Reg(i+1).Side, 'Right')
            AltErr(i) = 1; 
        else
            AltErr(i) = 0;
          end
    end
end

% if not alternating, use centroid positions to correct
if exist('Alt','var') == 1
    if sum(AltErr) > 0 && length(Reg) > 2
        for i = 1:length(Reg)
            % find average position of all centroids and use that location to determine L and R feet
            if Reg(i).Centroid(1) < MeanCentLoc % if horiz position of centroid is less than trendline -> LEFT
                Reg(i).Side = 'Left';
%                 Reg(i).SideErr = abs(Reg(i).Centroid(1) - MeanCentLoc);
            else
                Reg(i).Side = 'Right'; % if horiz position of centroid is greater than trendline -> RIGHT
%                 Reg(i).SideErr = abs(Reg(i).Centroid(1) - MeanCentLoc);
            end
        end
    end
end
clearvars AltErr

% if there are still mis classifications, use error from deviation line as indicator
for i = 1:length(Reg)-1
    if strcmp(Reg(i).Side, 'Left')
        if strcmp(Reg(i+1).Side, 'Right')
            AltErr(i) = 0; 
        else
            AltErr(i) = 1;
        end
    else
          if strcmp(Reg(i+1).Side, 'Right')
            AltErr(i) = 1; 
        else
            AltErr(i) = 0;
          end
    end
end

if exist('Alt','var') == 1
    a = find(AltErr);
    [~,ind] =  min([Reg(a:a+1).SideErr]);
    if strcmp(Reg(ind).Side, 'Left')
        Reg(ind).Side = 'Right';
    else
        Reg(ind).Side = 'Left';
    end
end

%% delete footstrikes that are too close to the ends of the walkway
clearvars ToDel
for j = 1:length(Reg)
    % automatically generate foot progression angle using orientation line through the centroid.
    ydist = Reg(j).MajorAxisLength .* sind(Reg(j).Orientation);
    xdist = Reg(j).MajorAxisLength .* cosd(Reg(j).Orientation);
    Pt1 = [Reg(j).Centroid(1) + (xdist/2), Reg(j).Centroid(2) - (ydist/2)];
    Pt2 = [Reg(j).Centroid(1) - (xdist/2), Reg(j).Centroid(2) + (ydist/2)];
    if Pt1(1) < 0 ||  Pt1(2) < 0 % if heel point is off the mat
        ToDel(j) = 1;
%     elseif Reg(j).MajorAxisLength <  PixelThresh .* 0.4 % if foot is too short
%         ToDel(j) = 1;
    elseif isnan(Pt1(1)) || isnan(Pt2(1)) % if heel/toe points dont exist
        ToDel(j) = 1;
    elseif Pt1(1) > PPDim.Length ||  Pt1(2) > PPDim.Length % if toe point is off the side of the mat
        ToDel(j) = 1;
    elseif  Reg(j).BoundingBox(2) < 1 % if heel point is off the end of the mat
        ToDel(j) = 1;
    elseif  Reg(j).BoundingBox(2) + Reg(j).BoundingBox(4) > PPDim.Length - 1 % if toe point is off the end of the mat
        ToDel(j) = 1;
    else
        ToDel(j) = 0;
    end
end

% delete steps that are not fully on mat
if exist('ToDel','var')
    Reg(ToDel==1) = [];
end
clearvars ToDel


%% Verify Selections
% figure; 
% contour(TM, 25);
% hold on;
% for i = 1:length(Reg)
%     if strcmp(Reg(i).Side, 'Left')
%         rectangle('Position',Reg(i).BoundingBox,'Color', 'r');
%     else
%         rectangle('Position',Reg(i).BoundingBox,'Color', 'g');
%     end
% end
% axis equal;
% 
% ReClassify = questdlg('Has Impressions correctly identified and labelled pressures?');
% if strcmp(ReClassify, 'No')
%     SideClass = questdlg('Are there mislabeled left/right steps?','Side Classification Error'); 
%     IdentificationErr = questdlg('Is there an improper identification?','Identification Error'); 
% end

%% Save data
Regions = Reg;

end



