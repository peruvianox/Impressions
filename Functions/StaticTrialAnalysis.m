function [FootLength] = StaticTrialAnalysis (StaticPPTrials, DisplayDim)

% if there are multiple static trials, choose one to be analyzed
% if length(StaticPPTrials) > 1
%     figure( 'Position', [100, 100, 600, 600]);
%     for i = 1:length(StaticPPTrials)
%         subplot(1, length(StaticPPTrials), i);
%         contour(StaticPPTrials(i).SumTM,250); 
%         axis equal; set(gca, 'XTick', []); set(gca, 'YTick', []);
%         title(['Trial ' num2str(i)]); 
%     end
%     
%     if length(StaticPPTrials) == 1
%         ListStr = {'Trial 1'};
%     elseif length(StaticPPTrials) == 2
%         ListStr = {'Trial 1', 'Trial 2'};
%     elseif length(StaticPPTrials) == 3
%         ListStr = {'Trial 1', 'Trial 2', 'Trial 3'};
%     end
% 
%     ChoiceStaticTrial = listdlg( 'PromptString','Choose a static trial to analyze','ListString',ListStr); 
%     StaticPPTrials = StaticPPTrials(ChoiceStaticTrial); 
%     close; 
% end

%% Select left and Right feet
figure('Position',  [DisplayDim(1) DisplayDim(2) 250 DisplayDim(4)]);
contour(StaticPPTrials.SumTM,100);
title('Static Trial');
set(gca, 'XTick', []);
set(gca, 'YTick', []);

%% manual foot length
% LEFT
uiwait(msgbox(['Please select the heel and toe points for the LEFT foot'], 'Instructions'));
Locs = ginput(2); 
StaticPPTrials.LeftHeel = Locs(1,:);  
StaticPPTrials.LeftToe = Locs(2,:);  

uiwait(msgbox(['Please select the heel and toe points for the RIGHT foot'], 'Instructions'));
Locs = ginput(2); 
StaticPPTrials.RightHeel = Locs(1,:);  
StaticPPTrials.RightToe = Locs(2,:);  

%% Calculate foot length
% left 
LeftFootLength = norm(StaticPPTrials.LeftToe - StaticPPTrials.LeftHeel) / 2; 
% right
RightFootLength = norm(StaticPPTrials.RightToe - StaticPPTrials.RightHeel) / 2; 

% average lengths
FootLength = mean([LeftFootLength, RightFootLength]); 


%% Auto foot length
% not curently functional

% % LEFT
% uiwait(msgbox(['Please select the left foot'], 'Instructions'));
% h = imrect;
% StaticPPTrials.Lrect = round(getPosition(h)); % draw box around step
% 
% % RIGHT
% uiwait(msgbox(['Please select the right foot'], 'Instructions'));
% h = imrect;
% StaticPPTrials.Rrect = round(getPosition(h)); % draw box around step
% 
% %% extract left and right feet from area
% % LEFT
% Left = StaticPPTrials.Lrect(1); 
% Right =  StaticPPTrials.Lrect(1) +  StaticPPTrials.Lrect(3);
% Lo = StaticPPTrials.Lrect(2); 
% Hi = StaticPPTrials.Lrect(2) +  StaticPPTrials.Lrect(4);
%  
% LogMat = logical(StaticPPTrials.SumTM(Lo:Hi, Left:Right) > 0);
% StaticPPTrials.LLog = zeros(size(StaticPPTrials.SumTM));
% StaticPPTrials.LLog(Lo:Hi, Left:Right) = LogMat;
% 
% clearvars Left Right Lo Hi LogMat
% 
% % RIGHT
% Left = StaticPPTrials.Rrect(1); 
% Right =  StaticPPTrials.Rrect(1) +  StaticPPTrials.Rrect(3);
% Lo = StaticPPTrials.Rrect(2); 
% Hi = StaticPPTrials.Rrect(2) +  StaticPPTrials.Rrect(4);
%  
% LogMat = logical(StaticPPTrials.SumTM(Lo:Hi, Left:Right) > 0);
% StaticPPTrials.RLog = zeros(size(StaticPPTrials.SumTM));
% StaticPPTrials.RLog(Lo:Hi, Left:Right) = LogMat;
% 
% clearvars Left Right Lo Hi LogMat
%  
% % uncomment to plot
% % figure; subplot(121);
% % contour(StaticPPTrials.LLog); 
% % title('Left Static'); 
% % subplot(122);
% % contour(StaticPPTrials.RLog)
% % title('Right Static');
% 
% 
% %% create logical array of connected areas
% % LEFT
% CC = bwconncomp(StaticPPTrials.LLog);
% BW = labelmatrix(CC);
% BW(BW>0) = 1; 
% StaticPPTrials.LReg = regionprops(BW,'BoundingBox', 'Orientation', 'Centroid');
% clearvars CC BW
% 
% % RIGHT
% 
% 
% %% make initial foot progression line
% % draw it shorter than necessary
% TinyFt = SubjDemo(5).Choice * 0.6;  % make foot 60% of original length and will search out to find edges
% TinyFt = TinyFt *2;  % change foot from cm to half-cm - same units as plate
% % get unit vector from orientation
% StaticPPTrials.LUV(2) = sind(StaticPPTrials.LReg.Orientation); 
% StaticPPTrials.LUV(1) = cosd(StaticPPTrials.LReg.Orientation); 
% if StaticPPTrials.LReg.Orientation < 90
%     StaticPPTrials.LUV(1) = -StaticPPTrials.LUV(1); 
% end
% 
% %% Determine length of foot by searching until there are no 
% figure; hold on; 
% contour(StaticPPTrials.SumTM, 100); 
% axis equal;
% 
% % LEFT
% ydist = TinyFt .* sind(StaticPPTrials.LReg.Orientation);
% xdist = TinyFt .* cosd(StaticPPTrials.LReg.Orientation);
% Pt1 = [StaticPPTrials.LReg.Centroid(1) + (xdist/2), StaticPPTrials.LReg.Centroid(2) - (ydist/2)];
% Pt2 = [StaticPPTrials.LReg.Centroid(1) - (xdist/2), StaticPPTrials.LReg.Centroid(2) + (ydist/2)];
% StaticPPTrials.LHT = [Pt1(1) Pt2(1); Pt1(2) Pt2(2)];
% StaticPPTrials.L_FPLine = line([Pt1(1) Pt2(1)],[Pt1(2) Pt2(2)] , 'Color','k');
% 
% % estimate medial and lateral border fo the foot 
% % TOE
% StaticPPTrials.LMedEstT(1,2) = StaticPPTrials.LHT(1,2) + sind(StaticPPTrials.LReg.Orientation)*(TinyFt*0.4); % X loc
% StaticPPTrials.LMedEstT(2,2) = StaticPPTrials.LHT(2,2) + cosd(StaticPPTrials.LReg.Orientation)*(TinyFt*0.4); % Y loc
% 
% StaticPPTrials.LLatEstT(1,2) = StaticPPTrials.LHT(1,2) - sind(StaticPPTrials.LReg.Orientation)*(TinyFt*0.4); % X loc
% StaticPPTrials.LLatEstT(2,2) = StaticPPTrials.LHT(2,2) - cosd(StaticPPTrials.LReg.Orientation)*(TinyFt*0.4); % Y loc
% 
% plot(StaticPPTrials.LMedEstT(1,2), StaticPPTrials.LMedEstT(2,2), '.m');
% plot(StaticPPTrials.LLatEstT(1,2), StaticPPTrials.LLatEstT(2,2), '.m');
% 
% %% do zone searching 
% SrchArea = 2; 
% [m,n] = size(StaticPPTrials.SumTM); 
% % LEFT
% % TOE
% for i = 1:0.5:25
%     a = [StaticPPTrials.LMedEstT(:,2)+(i*StaticPPTrials.LUV(1))];
%     b = [StaticPPTrials.LLatEstT(:,2)+(i*StaticPPTrials.LUV(2))];
%     SrchAreaX = [a(1)+SrchArea; a(1)-SrchArea; a(2)+SrchArea; a(2)-SrchArea];
%     SrchAreaY = [b(1)+SrchArea; b(1)-SrchArea; b(2)+SrchArea; b(2)-SrchArea];
%     Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
%     Region = StaticPPTrials.LLog .* Zone;
%     hold on; plot(a, b, '.m'); 
%     if sum(sum(Region)) == 0
%         line(a,b,'Color','m');
%         break
%         if a < DynamicPPTrials(Selection(q)).Lht{j}(1,1)
%             DynamicPPTrials(Selection(q)).LLat{j} = [a b];
%         else
%             DynamicPPTrials(Selection(q)).LMed{j} = [a b];
%         end
%     end
% end
% 
% 
% %%
% 
% % search in the other direction
% for i = 1:1:20
%     a = [DynamicPPTrials(Selection(q)).Lht{j}(:,1)-(i*DynamicPPTrials(Selection(q)).LhtInV{j}(:,1))];
%     b = [DynamicPPTrials(Selection(q)).Lht{j}(:,2)-(i*DynamicPPTrials(Selection(q)).LhtInV{j}(:,2))];
%     SrchAreaX = [a(1)+2; a(1)-2; a(2)+2; a(2)-2];
%     SrchAreaY = [b(1)+2; b(1)-2; b(2)+2; b(2)-2];
%     Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
%     Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
%     if sum(sum(Region)) == 0
%         subplot(1,NumDynTrials, Selection(q));
%         line(a,b,'Color','k');
%         break
%     end
% end
% if a < DynamicPPTrials(Selection(q)).Lht{j}(1,1)
%     DynamicPPTrials(Selection(q)).LLat{j} = [a b];
% else
%     DynamicPPTrials(Selection(q)).LMed{j} = [a b];
% end
% 
% 


end