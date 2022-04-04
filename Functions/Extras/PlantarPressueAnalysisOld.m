%% Plantar Pressure Analysis
%% To Do 

% filename = 'NR';

% convert cells to cm! 
% fix COPI issue with missing data
% static trial
% intoeing issues
% normalize pressure to body wt
% normalize walking surface area to standing area
% quicker inputs
% COP data
% correct arch segmentation = 1/3rds not with toes

% STEP WIDTH!

%% Create .mat file for image processing and save

% AnalyzePPImages = plantar_Excel2struct_PRS;


%% Get data from .mat file and organize for analysis

[filename, ~] = uigetfile('*.mat','Select .mat file for PP Analysis.');
% filename = 'XX.mat';
load(filename);

TimeMatrix = AnalyzePPImages.time_matrix(:,:,:);
TM = TimeMatrix;
SumTM = sum(TimeMatrix,3); 
MAT = AnalyzePPImages;

%% Create contour figure
Cont = figure;
set(Cont, 'Position', [100, 100, 400, 600]);
contour(SumTM,100); hold on; 
axis equal; 
% if image is upside down - correct it
button = questdlg('Does the image need to be flipped?');
if strcmp(button,'Yes') ==1
    Flip = 1;
    close
    SumTM = imrotate(SumTM,180);
    TM = imrotate(TM,180);
    Cont = figure;
    set(Cont, 'Position', [100, 100, 300, 600]);
    contour(SumTM,25); hold on;
    axis equal;
else
    Flip = 0;
end

%% Select feet for analysis
% draw foot progression angle lines
% LEFT side
Question = 'How many LEFT feet would you like to analyze?';
NumLeftButton = questdlg(Question, 'NumLeftFeet','1','2','3','1');
if NumLeftButton == '1'
NumLeft = 1;
end
if NumLeftButton == '2'
NumLeft = 2;
end
if NumLeftButton == '3'
NumLeft = 3;
end

for i = 1:NumLeft
uiwait(msgbox({'Using the cursor & displayed figure, select points just posterior to the heel and just anterior to the 2nd toe of the LEFT foot'}));
Lht{i} = ginput(2); 
line(Lht{i}(:,1),Lht{i}(:,2),'Color','k','LineWidth',1.5);
end

% RIGHT side
Question = 'How many RIGHT feet would you like to analyze?';
NumRightButton = questdlg(Question, 'NumLeftFeet','1','2','3','1');
if NumRightButton == '1'
NumRight = 1;
end
if NumRightButton == '2'
NumRight = 2;
end
if NumRightButton == '3'
NumRight = 3;
end

for i = 1:NumRight
uiwait(msgbox({'Using the cursor & displayed figure, select points just posterior to the heel and just anterior to the 2nd toe of the RIGHT foot'}));
Rht{i} = ginput(2); 
line(Rht{i}(:,1),Rht{i}(:,2),'Color','k','LineWidth',1.5);
end

%% calculate prog angles
for i = 1:NumLeft
LProgAng{i} = round(abs(180/pi*(atan((Lht{i}(2,1)-Lht{i}(1,1))/(Lht{i}(2,2)-Lht{i}(1,2)))))); 
% display prog angles on plot 
Ltxt = ['L Prog = ',num2str(LProgAng{i})]; 
text(Lht{i}(2,1), Lht{i}(2,2)+10, Ltxt,'VerticalAlignment','top', 'FontSize',7); 
% calculate inverted slope of prog angles
LhtV = [(Lht{i}(2,1)-Lht{i}(1,1)) (Lht{i}(2,2)-Lht{i}(1,2))]  / norm([(Lht{i}(2,1)-Lht{i}(1,1)) (Lht{i}(2,2)-Lht{i}(1,2))]);
LProg3Slope{i} = atan(-1/((Lht{i}(2,1)-Lht{i}(1,1))/(Lht{i}(2,2)-Lht{i}(1,2)))) .* 180/pi;
LhtInV{i} = [-LhtV(2) LhtV(1)];
end

for i = 1:NumRight
RProgAng{i} = round(abs(180/pi*(atan((Rht{i}(2,1)-Rht{i}(1,1))/(Rht{i}(2,2)-Rht{i}(1,2)))))); 
% display prog angles on plot
Rtxt = ['R Prog = ',num2str(RProgAng{i})]; 
text(Rht{i}(2,1), Rht{i}(2,2)+10, Rtxt,'VerticalAlignment','top', 'FontSize',7); 
% calculate inverted slope of prog angles
RProg3Slope{i} = atan(-1/((Rht{i}(2,1)-Rht{i}(1,1))/(Rht{i}(2,2)-Rht{i}(1,2)))) .* 180/pi;
RhtV = [(Rht{i}(2,1)-Rht{i}(1,1)) (Rht{i}(2,2)-Rht{i}(1,2))]  / norm([(Rht{i}(2,1)-Rht{i}(1,1)) (Rht{i}(2,2)-Rht{i}(1,2))]);
RhtInV{i} = [-RhtV(2) RhtV(1)];
end

%% search perpendicularly to each prog line 
MainMask = SumTM; 
MaskLog = logical(MainMask); 
[m,n] = size(MainMask);
% Search in one direction for the lateral border of the foot
% LEFT
for j = 1:NumLeft
    for i = 1:1:20
        a = [Lht{j}(:,1)+(i*LhtInV{j}(:,1))];
        b = [Lht{j}(:,2)+(i*LhtInV{j}(:,2))];
        SrchAreaX = [a(1)+2; a(1)-2; a(2)+2; a(2)-2];
        SrchAreaY = [b(1)+2; b(1)-2; b(2)+2; b(2)-2];
        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
        Region = MainMask .* Zone;
        if sum(sum(Region)) == 0
            line(a,b,'Color','k');
            break
        end
    end
    if a < Lht{j}(1,1)
        LLat{j} = [a b];
    else
        LMed{j} = [a b];
    end
    % search in the other direction
    for i = 1:1:20
        a = [Lht{j}(:,1)-(i*LhtInV{j}(:,1))];
        b = [Lht{j}(:,2)-(i*LhtInV{j}(:,2))];
        SrchAreaX = [a(1)+2; a(1)-2; a(2)+2; a(2)-2];
        SrchAreaY = [b(1)+2; b(1)-2; b(2)+2; b(2)-2];
        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
        Region = MainMask .* Zone;
        if sum(sum(Region)) == 0
            line(a,b,'Color','k');
            break
        end
    end
    if a < Lht{j}(1,1)
        LLat{j} = [a b];
    else
        LMed{j} = [a b];
    end
end

%% RIGHT
for j  = 1:NumRight
    for i = 1:1:20
        a = [Rht{j}(:,1)+(i*RhtInV{j}(:,1))];
        b = [Rht{j}(:,2)+(i*RhtInV{j}(:,2))];
        SrchAreaX = [a(1)+2; a(1)-2; a(2)+2; a(2)-2];
        SrchAreaY = [b(1)+2; b(1)-2; b(2)+2; b(2)-2];
        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
        Region = MainMask .* Zone;
        if sum(sum(Region)) == 0
            line(a,b,'Color','k');
            break
        end
    end
    if a > Rht{j}(1,1)
        RLat{j} = [a b];
    else
        RMed{j} = [a b];
    end
    % search in the other direction
    for i = 1:1:20
        a = [Rht{j}(:,1)-(i*RhtInV{j}(:,1))];
        b = [Rht{j}(:,2)-(i*RhtInV{j}(:,2))];
        SrchAreaX = [a(1)+2; a(1)-2; a(2)+2; a(2)-2];
        SrchAreaY = [b(1)+2; b(1)-2; b(2)+2; b(2)-2];
        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
        Region = MainMask .* Zone;
        if sum(sum(Region)) == 0
            line(a,b,'Color','k');
            break
        end
    end
    if a > Rht{j}(1,1)
        RLat{j} = [a b];
    else
        RMed{j} = [a b];
    end
end

clearvars Zone Region SrchAreaX SrchAreaY i button a b NumLeftButton NumRightButton

%% 1/3 sectioning
% LEFT
for i = 1:NumLeft
    % find the 1/3 of foot prog line points
    LProg13{i} = [Lht{i}(1,1) + ((Lht{i}(2,1) - Lht{i}(1,1))/3),  Lht{i}(1,2) + ((Lht{i}(2,2) - Lht{i}(1,2))/3)];
    LProg23{i} = [Lht{i}(1,1) + 2*((Lht{i}(2,1) - Lht{i}(1,1))/3),  Lht{i}(1,2) + 2*((Lht{i}(2,2) - Lht{i}(1,2))/3)];
    % define the inverse slope for the 1/3 sectioning lines
    LProg13line{i} = [LProg13{i} + LhtInV{i}; LProg13{i} - LhtInV{i}];
    LProg23line{i} = [LProg23{i} + LhtInV{i}; LProg23{i} - LhtInV{i}];
    % find intersections between third sectioning lines and lat/med borders
    LProg13Lat{i} = linlinintersect([LProg13line{i};LLat{i}]);
    LProg13Med{i} = linlinintersect([LProg13line{i};LMed{i}]);
    LProg23Lat{i} = linlinintersect([LProg23line{i};LLat{i}]);
    LProg23Med{i} = linlinintersect([LProg23line{i};LMed{i}]);
    % show the third-sectioning lines
    line([LProg13Lat{i}(1) LProg13Med{i}(1) ],[LProg13Lat{i}(2) LProg13Med{i}(2)],'Color','k');
    line([LProg23Lat{i}(1) LProg23Med{i}(1)],[LProg23Lat{i}(2) LProg23Med{i}(2)],'Color','k');
end

% RIGHT
for i = 1:NumRight
    % find the 1/3 of foot prog line points
    RProg13{i} = [Rht{i}(1,1) + ((Rht{i}(2,1) - Rht{i}(1,1))/3),  Rht{i}(1,2) + ((Rht{i}(2,2) - Rht{i}(1,2))/3)];
    RProg23{i} = [Rht{i}(1,1) + 2*((Rht{i}(2,1) - Rht{i}(1,1))/3),  Rht{i}(1,2) + 2*((Rht{i}(2,2) - Rht{i}(1,2))/3)];
    % define the inverse slope for the 1/3 sectioning lines
    RProg13line{i} = [RProg13{i} + RhtInV{i}; RProg13{i} - RhtInV{i}];
    RProg23line{i} = [RProg23{i} + RhtInV{i}; RProg23{i} - RhtInV{i}];
    % find intersections between third sectioning lines and lat/med borders
    RProg13Lat{i} = linlinintersect([RProg13line{i};RLat{i}]);
    RProg13Med{i} = linlinintersect([RProg13line{i};RMed{i}]);
    RProg23Lat{i} = linlinintersect([RProg23line{i};RLat{i}]);
    RProg23Med{i} = linlinintersect([RProg23line{i};RMed{i}]);
    % show the third-sectioning lines
    line([RProg13Lat{i}(1) RProg13Med{i}(1) ],[RProg13Lat{i}(2) RProg13Med{i}(2)],'Color','k');
    line([RProg23Lat{i}(1) RProg23Med{i}(1)],[RProg23Lat{i}(2) RProg23Med{i}(2)],'Color','k');
end

clearvars L13 L23 R13 R23 LProg3Slope RProg3Slope Ltxt Rtxt LhtV RhtV

%% find intersection points for zone classification
% Split up the LEFT Foot into the 6 zones
for i = 1:NumLeft
LFootPts{i} = round([LLat{i};LMed{i}]);
% Lateral and Medial Heel
Zones.L(i).LatHeel = round([LLat{i}(1,:); Lht{i}(1,:); LProg13{i}(1,:); LProg13Lat{i}(1,:)]);
Zones.L(i).MedHeel = round([LMed{i}(1,:); Lht{i}(1,:); LProg13{i}(1,:); LProg13Med{i}(1,:)]);
% Lateral and Medial Arch
Zones.L(i).LatArch = round([LProg13Lat{i}(1,:); LProg13{i}(1,:); LProg23{i}(1,:); LProg23Lat{i}(1,:)]);
Zones.L(i).MedArch = round([LProg13Med{i}(1,:); LProg23Med{i}(1,:); LProg23{i}(1,:); LProg13{i}(1,:)]);
% Lateral and Medial Forefoot
Zones.L(i).LatFore = round([LProg23{i}(1,:); Lht{i}(2,:); LLat{i}(2,:); LProg23Lat{i}(1,:)]);
Zones.L(i).MedFore = round([LProg23Med{i}(1,:); LMed{i}(2,:); Lht{i}(2,:); LProg23{i}(1,:)]);
end

% Split up the RIGHT Foot into the 6 zones
for i = 1:NumRight
RFootPts{i} = round([RLat{i};RMed{i}]);
% Lateral and Medial Heel
Zones.R(i).LatHeel = round([RLat{i}(1,:); RProg13Lat{i}(1,:); RProg13{i}(1,:); Rht{i}(1,:)]);
Zones.R(i).MedHeel = round([RMed{i}(1,:); Rht{i}(1,:); RProg13{i}(1,:); RProg13Med{i}(1,:)]);
% Lateral and Medial Arch
Zones.R(i).LatArch = round([RProg13Lat{i}(1,:); RProg23Lat{i}(1,:); RProg23{i}(1,:); RProg13{i}(1,:)]);
Zones.R(i).MedArch = round([RProg13{i}(1,:); RProg23{i}(1,:); RProg23Med{i}(1,:); RProg13Med{i}(1,:)]);
% Lateral and Medial Forefoot
Zones.R(i).LatFore = round([RProg23Lat{i}(1,:); RLat{i}(2,:); Rht{i}(2,:); RProg23{i}(1,:)]);
Zones.R(i).MedFore = round([RProg23Med{i}(1,:);  RProg23{i}(1,:); Rht{i}(2,:); RMed{i}(2,:)]);
end

% *uncomment to check to make sure the polygon points are correct 
% plot(RLatHeel(:,1),RLatHeel(:,2), 'om'); 
% plot(RMedFore(:,1),RMedFore(:,2), '*m'); 

%% Identify foot strikes for Spatio Temporal Measurements
Reg = SpatioTemporalPP(AnalyzePPImages,Flip); 

%% calculate step width and step length
% StepWidth = 

%% Generate a L & R plots with COP data
% determine if image was flipped. If yes, flip the COP measures
if Flip == 1
    COP = (MAT.center_of_pressure);
    COP(:,1) = 64 - COP(:,1);
    COP(:,2) = 257 - COP(:,2);
else
    COP = MAT.center_of_pressure;
end% define x and y COP
COPx = COP(:,1); 
COPy = COP(:,2); 

% LEFT side
for i = 1:NumLeft
% Prepare to merge COP and pressure plots
LFootBox{i} = [min(LFootPts{i}(:,1)) , min(LFootPts{i}(:,2)) ; max(LFootPts{i}(:,1)) , max(LFootPts{i}(:,2))];
% Left Foot - define times when COP is in foot area
Del(:,1) = COP(:,1) < min(LFootBox{i}(:,1));
Del(:,2) = COP(:,1) > max(LFootBox{i}(:,1));
Del(:,3) = COP(:,2) < min(LFootBox{i}(:,2));
Del(:,4) = COP(:,2) > max(LFootBox{i}(:,2));
SumDel = sum(Del,2); 
% Create logical of times when the cop is in the foot area
LogDel = logical(~SumDel);
% Apply logical to COP, and subtract out lower left corner of plotting zone
COP_L{i} = [COPx(LogDel) - LFootBox{i}(1,1), COPy(LogDel) - LFootBox{i}(1,2)];
end

% RIGHT side
for i = 1:NumRight
% Prepare to merge COP and pressure plots
RFootBox{i} = [min(RFootPts{i}(:,1)) , min(RFootPts{i}(:,2)) ; max(RFootPts{i}(:,1)) , max(RFootPts{i}(:,2))];
% Right Foot  - define times when COP is in foot area
Del(:,1) = COP(:,1) < min(RFootBox{i}(:,1));
Del(:,2) = COP(:,1) > max(RFootBox{i}(:,1));
Del(:,3) = COP(:,2) < min(RFootBox{i}(:,2));
Del(:,4) = COP(:,2) > max(RFootBox{i}(:,2));
SumDel = sum(Del,2); 
% Create logical of times when the cop is in the foot area
LogDel = logical(~SumDel);
% Apply logical to COP, and subtract out lower left corner of plotting zone
COP_R{i} = [COPx(LogDel) - RFootBox{i}(1,1), COPy(LogDel) - RFootBox{i}(1,2)];
end 

% Adjust foot boxes if necessary
for j = 1:NumLeft
    for i = 1:2 % if LEFT side box is defined outside of PP matrix
        if LFootBox{j}(i,1) < 0  
            LFootBox{j}(i,1) = 0;
        end
        if LFootBox{j}(i,1) > 63
            LFootBox{j}(i,1) = 63;
        end
        if LFootBox{j}(i,2) < 0
            LFootBox{j}(i,2) = 0;
        end
        if LFootBox{j}(i,2) > 256
            LFootBox{j}(i,2) = 256;
        end
    end
end
    
for j = 1:NumRight
    for i = 1:2 % if RIGHT side box is defined outside of PP matrix
        if RFootBox{j}(i,1) < 0  
            RFootBox{j}(i,1) = 0;
        end
        if RFootBox{j}(i,1) > 63
            RFootBox{j}(i,1) = 63;
        end
        if RFootBox{j}(i,2) < 0
            RFootBox{j}(i,2) = 0;
        end
        if RFootBox{j}(i,2) > 256
            RFootBox{j}(i,2) = 256;
        end
    end
end

% define whole plot foot contours 
% LEFT side
if NumLeft == 1
    Lcont = 1;
else
    Question = 'How many LEFT foot contour and COP would you like to plot?';
    NumLeftCont = questdlg(Question, 'LeftFootPlot','1','2','3','1');
    if NumLeftCont == '1'
        Lcont = 1;
    end
    if NumLeftCont == '2'
        Lcont = 2;
    end
    if NumLeftCont == '3'
        Lcont = 3;
    end
end
LFootContour  = SumTM(LFootBox{Lcont}(1,2):LFootBox{Lcont}(2,2),LFootBox{Lcont}(1,1):LFootBox{Lcont}(2,1));

% RIGHT side
if NumRight == 1
    Rcont = 1;
else
    Question = 'How many RIGHT foot contour and COP would you like to plot?';
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
RFootContour  = SumTM(RFootBox{Rcont}(1,2):RFootBox{Rcont}(2,2),RFootBox{Rcont}(1,1):RFootBox{Rcont}(2,1));

% Create COP figure
COPplots = figure; 

% Left Foot Plot
subplot(1,2,1); 
contour(LFootContour, 100); 
axis equal; hold on;
plot(COP_L{Lcont}(:,1), COP_L{Lcont}(:,2), 'k.', 'MarkerSize',10);
% Right Foot Plot
subplot(1,2,2); 
contour(RFootContour, 100); 
axis equal; hold on;
plot(COP_R{Rcont}(:,1), COP_R{Rcont}(:,2), 'k.', 'MarkerSize',10);

% Create COP analyses
% COP trendline vs prog angle
% polynomial fit on the 1st degree
% LEFT
pL = polyfit(COP_L{Lcont}(:,1), COP_L{Lcont}(:,2),1);
ypL = polyval(pL, COP_L{Lcont}(:,1)); 

xpL = linspace(Lht{Lcont}(1,1) - LFootBox{1}(1,1), Lht{Lcont}(2,1) - LFootBox{1}(1,1));
ypL = polyval(pL, xpL); 
subplot (121); 
%plot (COP_L{Lcont}(:,1), ypL,'r', 'LineWidth',3);

% re-do box reference from global orientation
plot (xpL, ypL,'r', 'LineWidth',3);

% RIGHT
pR = polyfit(COP_R{Rcont}(:,1), COP_R{Rcont}(:,2),1);
ypR = polyval(pR, COP_R{Rcont}(:,1)); 
subplot (122); 
plot (COP_R{Rcont}(:,1), ypR,'r','LineWidth',3);

clearvars Del SumDel LogDel

figure;
plot(COP_R{1}(:,1), COP_R{1}(:,2),'.k'); 
hold on;
plot(COP_R{1}(:,1), ypR,'.b'); 

%% Calculate Pressures and Areas 
% LEFT Foot
for i = 1:NumLeft
% Create Masks
% Heel
Mask.L(i).LatHeel = poly2mask(Zones.L(i).LatHeel(:,1),Zones.L(i).LatHeel(:,2),m,n); 
LLatHeelPress{i} = sum(sum(Mask.L(i).LatHeel .* MainMask));
LLatHeelArea{i} = sum(sum(Mask.L(i).LatHeel .* MaskLog));
Mask.L(i).MedHeel = poly2mask(Zones.L(i).MedHeel(:,1),Zones.L(i).MedHeel(:,2),m,n); 
LMedHeelPress{i} = sum(sum(Mask.L(i).MedHeel .* MainMask));
LMedHeelArea{i} = sum(sum(Mask.L(i).MedHeel .* MaskLog));
% Arch
Mask.L(i).LatArch = poly2mask(Zones.L(i).LatArch(:,1),Zones.L(i).LatArch(:,2),m,n); 
LLatArchPress{i} = sum(sum(Mask.L(i).LatArch .* MainMask));
LLatArchArea{i} = sum(sum(Mask.L(i).LatArch .* MaskLog));
Mask.L(i).MedArch = poly2mask(Zones.L(i).MedArch(:,1),Zones.L(i).MedArch(:,2),m,n); 
LMedArchPress{i} = sum(sum(Mask.L(i).MedArch .* MainMask));
LMedArchArea{i} = sum(sum(Mask.L(i).MedArch .* MaskLog));
% Forefoot
Mask.L(i).LatFore = poly2mask(Zones.L(i).LatFore(:,1),Zones.L(i).LatFore(:,2),m,n); 
LLatForePress{i} = sum(sum(Mask.L(i).LatFore .* MainMask));
LLatForeArea{i} = sum(sum(Mask.L(i).LatFore .* MaskLog));
Mask.L(i).MedFore = poly2mask(Zones.L(i).MedFore(:,1),Zones.L(i).MedFore(:,2),m,n); 
LMedForePress{i} = sum(sum(Mask.L(i).MedFore .* MainMask));
LMedForeArea{i} = sum(sum(Mask.L(i).MedFore .* MaskLog));
% Pressures
LLatPress{i} = LLatForePress{i} + LLatArchPress{i} + LLatHeelPress{i};
LMedPress{i} = LMedForePress{i} + LMedArchPress{i} + LMedHeelPress{i};
LFootTotPress{i} = LLatPress{i} + LMedPress{i};
% Areas
LLatArea{i} = LLatForeArea{i} + LLatArchArea{i} + LLatHeelArea{i};
LMedArea{i} = LMedForeArea{i} + LMedArchArea{i} + LMedHeelArea{i};
LFootTotArea{i} = LLatArea{i} + LMedArea{i};
% Compute Foot Ratios
LMedRatio{i} = LMedArea{i}/LFootTotArea{i}; % Medial to Whole area
LLatRatio{i} = LLatArea{i}/LFootTotArea{i}; % % Lateral to Whole area
LRatioArea{i} =  LMedArea{i}/LLatArea{i}; % Medial to Lateral Area
LArchIndex{i} = (LMedArchArea{i}+LLatArchArea{i})/LFootTotArea{i};  % Arch to Whole Area
LCOPI{i} = LLatArchArea{i}/LMedArchArea{i}; % Lateral Arch to Medial Arch Area
end

%% RIGHT Foot
for i = 1:NumRight
% Heel
Mask.R(i).LatHeel = poly2mask(Zones.R(i).LatHeel(:,1),Zones.R(i).LatHeel(:,2),m,n); 
RLatHeelPress{i} = sum(sum(Mask.R(i).LatHeel .* MainMask));
RLatHeelArea{i} = sum(sum(Mask.R(i).LatHeel .* MaskLog));
Mask.R(i).MedHeel = poly2mask(Zones.R(i).MedHeel(:,1),Zones.R(i).MedHeel(:,2),m,n); 
RMedHeelPress{i} = sum(sum(Mask.R(i).MedHeel .* MainMask));
RMedHeelArea{i} = sum(sum(Mask.R(i).MedHeel .* MaskLog));
% Arch
Mask.R(i).LatArch = poly2mask(Zones.R(i).LatArch(:,1),Zones.R(i).LatArch(:,2),m,n); 
RLatArchPress{i} = sum(sum(Mask.R(i).LatArch .* MainMask));
RLatArchArea{i} = sum(sum(Mask.R(i).LatArch .* MaskLog));
Mask.R(i).MedArch = poly2mask(Zones.R(i).MedArch(:,1),Zones.R(i).MedArch(:,2),m,n); 
RMedArchPress{i} = sum(sum(Mask.R(i).MedArch .* MainMask));
RMedArchArea{i} = sum(sum(Mask.R(i).MedArch .* MaskLog));
% Forefoot
Mask.R(i).LatFore = poly2mask(Zones.R(i).LatFore(:,1),Zones.R(i).LatFore(:,2),m,n); 
RLatForePress{i} = sum(sum(Mask.R(i).LatFore .* MainMask));
RLatForeArea{i} = sum(sum(Mask.R(i).LatFore .* MaskLog));
Mask.R(i).MedFore = poly2mask(Zones.R(i).MedFore(:,1),Zones.R(i).MedFore(:,2),m,n); 
RMedForePress{i} = sum(sum(Mask.R(i).MedFore .* MainMask));
RMedForeArea{i} = sum(sum(Mask.R(i).MedFore .* MaskLog));  
% Pressures
RLatPress{i} = RLatForePress{i} + RLatArchPress{i} + RLatHeelPress{i};
RMedPress{i} = RMedForePress{i} + RMedArchPress{i} + RMedHeelPress{i};
RFootTotPress{i} = RLatPress{i} + RMedPress{i};
% Areas
RLatArea{i} = RLatForeArea{i} + RLatArchArea{i} + RLatHeelArea{i};
RMedArea{i} = RMedForeArea{i} + RMedArchArea{i} + RMedHeelArea{i};
RFootTotArea{i} = RLatArea{i} + RMedArea{i};
% Compute Foot Ratios
RMedRatio{i} = RMedArea{i}/RFootTotArea{i}; % Medial to Whole area
RLatRatio{i} = RLatArea{i}/RFootTotArea{i}; % % Lateral to Whole area
RRatioArea{i} =  RMedArea{i}/RLatArea{i}; % Medial to Lateral Area
RArchIndex{i} = (RMedArchArea{i}+RLatArchArea{i})/RFootTotArea{i};  % Arch to Whole Area
RCOPI{i} = RLatArchArea{i}/RMedArchArea{i}; % Lateral Arch to Medial Arch Area
end
%% Organization for Export
PressureRows = {'Press';'Fore';'Arch';'Heel';'Side';'Foot'};
AreaRows = {'Areas';'Fore';'Arch';'Heel';'Side';'Foot'};
RatioRows = {'Ratios';'Med Ratio';'Lat Ratio';'Ratio Area';'Arch Index';'COPI'};

% LEFT side
for i = 1:NumLeft
LFootPress{i} = {'Lat L','Med L';LLatForePress{i} LMedForePress{i}; LLatArchPress{i} LMedArchPress{i}; LLatHeelPress{i} LMedHeelPress{i}; LLatPress{i} LMedPress{i}; LFootTotPress{i} LFootTotPress{i}};
LFootAreas{i} = {'Lat L','Med L';LLatForeArea{i} LMedForeArea{i}; LLatArchArea{i} LMedArchArea{i}; LLatHeelArea{i} LMedHeelArea{i}; LLatArea{i} LMedArea{i}; LFootTotArea{i} LFootTotArea{i}};
LRatios{i} = {'Left';LMedRatio{i};LLatRatio{i};LRatioArea{i};LArchIndex{i};LCOPI{i}};
end

% RIGHT side
for i = 1:NumRight
RFootPress{i} = {'Lat R','Med R';RLatForePress{i} RMedForePress{i}; RLatArchPress{i} RMedArchPress{i}; RLatHeelPress{i} RMedHeelPress{i}; RLatPress{i} RMedPress{i}; RFootTotPress{i} RFootTotPress{i}};
RFootAreas{i} = {'Lat R','Med R';RLatForeArea{i} RMedForeArea{i}; RLatArchArea{i} RMedArchArea{i}; RLatHeelArea{i} RMedHeelArea{i}; RLatArea{i} RMedArea{i}; RFootTotArea{i} RFootTotArea{i}};
RRatios{i} = {'Right';RMedRatio{i};RLatRatio{i};RRatioArea{i};RArchIndex{i};RCOPI{i}};
end

%% Convert measurements into units
% RSscan plate has 16384 sensors with dimensions of 5.08 mm x 7.62 mm.
% Each element of time_matrix also has a dimension of 5.0794 mm x 7.6172
% mm. Therefore, the area of each element containing the plantar pressure
% measurement is 38.6908 mm^2 = 0.386908 cm^2 = 3.86908 x 10^-5 m^2
%
% Since PRESSURE = FORCE/AREA, then FORCE = PRESSURE * (3.86908 x 10^-5 m^2)
Area = 0.386908; % cm^2
Pressure = 3.86908 * 10^-5; % m^2 - what are the force units? N? lbs?

for i = 1:NumLeft
LFootAreasN{i} = cell2mat(LFootAreas{i}(2:6,:)) * Area;
LFootPressN{i} = cell2mat(LFootPress{i}(2:6,:)) * Pressure; 
end

for i = 1:NumRight
RFootAreasN{i} = cell2mat(RFootAreas{i}(2:6,:)) * Area;
RFootPressN{i} = cell2mat(RFootPress{i}(2:6,:)) * Pressure; 
end

%% Determine what step(s) should be outputed
Pressures{i} = [PressureRows,LFootPress{i},RFootPress{i}];
Areas = [AreaRows,LFootAreas,RFootAreas];
Ratios = [RatioRows,LRatios,RRatios];

clearvars LLatForePress LMedForePress  LLatArchPress LMedArchPress LLatHeelPress LMedHeelPress LLatPress LMedPress LFootTotPress LFootTotPress ...
    RLatForePress RMedForePress RLatArchPress RMedArchPress RLatHeelPress RMedHeelPress RLatPress RMedPress RFootTotPress RFootTotPress ...
    LLatForeArea LMedForeArea LLatArchArea LMedArchArea LLatHeelArea LMedHeelArea LLatArea LMedArea LFootTotArea LFootTotArea ...
    RLatForeArea RMedForeArea RLatArchArea RMedArchArea RLatHeelArea RMedHeelArea RLatArea RMedArea RFootTotArea RFootTotArea ...
    LMedRatio LLatRatio LRatioArea LArchIndex LCOPI ...
    RMedRatio RLatRatio RRatioArea RArchIndex RCOPI


%% Trial and Subject Qualities
% is it a static trial?
button = questdlg('Is this a static trial?');
% if yes, then compute standing foot area for each foot for normalization
if strcmp(button,'Yes') == 1
    Type = 'Static';
else
    Type = 'Dynamic';
end
% Age of subject
SubjectAge = input('What is the age of the subject?');

% Mass of subject
SubjectMass = input('What is the mass of the subject?'); 

Qualities1 = {'Trial Type';'Subject Age';'L Foot Prog';'R Foot Prog';'Subject Mass';0};
Qualities2 = {Type; SubjectAge; LProgAng; RProgAng;SubjectMass;0};

Qualities = [Qualities1, Qualities2];

clearvars Qualities1 Qualities2 SubjectAge SubjectMass Type

%%  Export data
% Output = [Qualities, Pressures, Areas, Ratios]; 
% % save as excel spreadsheet
% %PPOutput = regexprep(filename,'.mat','.xlsx');
% PPOutput = strcat(filename,'.xlsx');
% xlswrite(PPOutput,Output);
% % save as figure
% %PPOutput = regexprep(filename,'.mat','.jpg');
% PPOutput = strcat(filename,'.jpg');
% saveas(Cont,PPOutput);


