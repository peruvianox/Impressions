function[TimeNorm, CoP, NormSum] = NormalizeStep(Step, FPAngle, Side)
% this function will input a 3D matrix of foot pressure data, a general
% foot progression angle, and which side it is on, and normalized the step
% into a 100x100x100 matrix of the foot in which each dimension is
% normalized to 100% of its length, width, and duration.

%% initializing
[~, ~, p] = size(Step);
TargetSize = [100 100];

%% alter inputs and settigns for left vs right foot
if strcmp(Side, 'Left')
    FPAngle = -FPAngle;
end

%% Rotate and Resize PP image dimensions
for i = 1:p % loop through number of frames
    RotStep(:,:,i) = imrotate(Step(:,:,i), FPAngle, 'bilinear');
end

% prepare to crop image by checking for empty rows and columns and find
% regions to crop
SumImage = sum(RotStep, 3);
ColSum = sum(SumImage, 1);
RowSum = sum(SumImage, 2);
RowBot = find(RowSum,1);
RowTop = length(RowSum) - find(flipud(RowSum),1);
ColLeft = find(ColSum, 1);
ColRight = length(ColSum) - find(fliplr(ColSum), 1);

% figure; subplot(121)
% contour(SumImage)
% Crpd = imcrop(SumImage, [ColLeft,RowBot,  ColRight, RowTop]);
% subplot(122);
% contour(Crpd);

for i = 1:p
    CropStep(:,:,i) = imcrop(RotStep(:,:,i), [ColLeft,RowBot,  ColRight, RowTop]);
    NormSteps(:,:,i) = imresize(CropStep(:,:,i), TargetSize);
    NormSteps(NormSteps<0) = 0;
end

%% uncomment plot rotated and resized step
% figure;
% for i = 1:p % loop through number of frames
%      subplot(131);
%     contour(Step(:,:,i));
%     title('Old');
%     subplot(132);
%     contour(RotStep(:,:,i));
%     title('Rotated');
%     subplot(133);
%     contour(NormSteps(:,:,i));
%     title('ReSized');
%     pause(0.05);
% end

%% resample to 100 frames in length
TimeNorm = zeros(100,100,100);
for i = 1:100
    Row = permute(NormSteps(i,:,:), [3 2 1]);
    NormRow= resample(Row, 100, p);
    TimeNorm(i,:,:) = permute(NormRow, [3 2 1]);
    TimeNorm(TimeNorm<0.1) = 0;
end

%% re-calcualte CoP and save
CoP = zeros(100,2);
for k = 1:100
    % manual CoP calculation
    Cx = sum(sum(TimeNorm(:,:,k)).*(1:100))/sum(sum(TimeNorm(:,:,k))); % x location centroid
    Cy = sum(sum(TimeNorm(:,:,k),2)'.*(1:100))/sum(sum(TimeNorm(:,:,k),2)'); % y location centroid
    CoP(k,:) = [Cx,Cy];
    clearvars Cx Cy
end

NormSum = sum(TimeNorm,3);

%% Quality check to make sure step was rotated the correct way
% Orientation for new normalized step
CC = bwconncomp(NormSum,8);
BW = labelmatrix(CC);
BW(BW>0) = 1; 
Reg = regionprops(BW, 'Orientation');

% Orientation for original step
OrigSum = sum(Step, 3); 
CC = bwconncomp(OrigSum,8);
BW = labelmatrix(CC);
BW(BW>0) = 1; 
OrigReg = regionprops(BW, 'Orientation');

% contour(OrigSum)

% % get original CoP
% for k = 1:p
%     % manual CoP calculation
%     OrigCx = sum(sum(Step(1:m,1:n,k)).*(1:n))/sum(sum(Step(1:m,1:n,k))); % x location centroid
%     OrigCy = sum(sum(Step(1:m,1:n,k),2)'.*(1:m))/sum(sum(Step(1:m,1:n,k),2)'); % y location centroid
%     OrigCoP(k,:) = [OrigCx,OrigCy];
%     clearvars Cx Cy
% end
% % linear trendline of original CoP
% OrigP = polyfit(OrigCoP(:,1), OrigCoP(:,2), 1);
% OrigX = 1:100;
% OrigY = polyval(OrigP,OrigX);
% 
% NewP = polyfit(CoP(:,1), CoP(:,2), 1); % new CoP
% LEFT
if strcmp(Side,'Left')
    if Reg.Orientation > OrigReg.Orientation % if new foot prog angle is greater than the other side, flip the other way.
        ReDo = 'Yes';
        FPAngle = -FPAngle; % flip angle to rotate
    else
        ReDo = 'No';
    end
elseif strcmp(Side, 'Right')
    if Reg.Orientation < OrigReg.Orientation % if new foot prog angle is greater than the other side, flip the other way.
        ReDo = 'Yes';
        FPAngle = -FPAngle; % flip angle to rotate
    else
        ReDo = 'No';
    end
end

%% Redo rotation the other direction if necessary
if strcmp(ReDo, 'Yes')
    %% Rotate and Resize PP image dimensions
    clearvars RotStep ColSum RowSum RowBot RowTop ColLeft ColRight CropStep NormSteps TimeNorm
    [~, ~, p] = size(Step);
    for i = 1:p % loop through number of frames
        RotStep(:,:,i) = imrotate(Step(:,:,i), FPAngle, 'bilinear');
    end
    
    % prepare to crop image by checking for empty rows and columns and find regions to crop
    SumImage = sum(RotStep, 3);
    ColSum = sum(SumImage, 1);
    RowSum = sum(SumImage, 2);
    RowBot = find(RowSum,1);
    RowTop = length(RowSum) - find(flipud(RowSum),1);
    ColLeft = find(ColSum, 1);
    ColRight = length(ColSum) - find(fliplr(ColSum), 1);
    
    % figure; subplot(121)
    % contour(SumImage)
    % Crpd = imcrop(SumImage, [ColLeft,RowBot,  ColRight, RowTop]);
    % subplot(122);
    % contour(Crpd);
    
    for i = 1:p
        CropStep(:,:,i) = imcrop(RotStep(:,:,i), [ColLeft,RowBot,  ColRight, RowTop]);
        NormSteps(:,:,i) = imresize(CropStep(:,:,i), TargetSize);
        NormSteps(NormSteps<0) = 0;
    end
    
    %% resample to 100 frames in length
    TimeNorm = zeros(100,100,100);
    for i = 1:100
        Row = permute(NormSteps(i,:,:), [3 2 1]);
        NormRow= resample(Row, 100, p);
        TimeNorm(i,:,:) = permute(NormRow, [3 2 1]);
        TimeNorm(TimeNorm<0.1) = 0;
    end
    
    %% re-calcualte CoP and save
    CoP = zeros(100,2);
    for k = 1:100
        % manual CoP calculation
        Cx = sum(sum(TimeNorm(:,:,k)).*(1:100))/sum(sum(TimeNorm(:,:,k))); % x location centroid
        Cy = sum(sum(TimeNorm(:,:,k),2)'.*(1:100))/sum(sum(TimeNorm(:,:,k),2)'); % y location centroid
        CoP(k,:) = [Cx,Cy];
        clearvars Cx Cy
    end
end

% figure;
% for i = 1:p % loop through number of frames
%      subplot(131);
%     contour(Step(:,:,i));
%     title('Old');
%     subplot(132);
%     contour(RotStep(:,:,i));
%     title('Rotated');
%     subplot(133);
%     contour(NormSteps(:,:,i));
%     title('ReSized');
%     pause(0.05);
% end

%% uncomment plot time-adjusted step with CoPs
% figure;
% for i = 1:100
%     hold on;
%     contour(TimeNorm(:,:,i));
%     plot(CoP(i,1), CoP(i,2), '.k');
%     if i > 3
%         plot(CoP(i-1,1), CoP(i-1,2), '.r');
%         plot(CoP(i-2,1), CoP(i-2,2), '.r');
%         plot(CoP(i-3,1), CoP(i-3,2), '.r');
%     end
%     pause(0.05);
%     clf;
% end

end

