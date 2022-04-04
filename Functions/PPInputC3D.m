function[C3Ddata,TrimSpat, DynamicPPTrials] = PPInputC3D(PPOutput, DynamicPPTrials, Subject, Regions)

%% load c3d data
C3Ddata = ProcessC3D;

%% Analyze static trials if present
for i = 1:length(C3Ddata)
    if strcmp(C3Ddata(i).Direction,'Static') == 1 || isfield(C3Ddata(i), 'Static') 
        StaticTrial(i) = 1;
    else
        StaticTrial(i) = 0;
    end
end
% StaticTrial = logical(StaticTrial);

if sum(StaticTrial) > 0
    StaticData = C3Ddata(StaticTrial==1); % save static data
    Aggregate = C3Ddata(1).Aggregate; % save aggregate data!
    C3Ddata(StaticTrial==1) = []; % delete static trial from dynamic structure
    C3Ddata(1).Static = StaticData; % save static trial as new structure
    C3Ddata(1).Aggregate = Aggregate; 
     % Static = find(StaticTrial);
    % end
    % if isempty(StaticData) == 0
    %         Static = find(StaticTrial);
    for i = 1:length(StaticData.Trajectories.sortedPoints)
        if strcmp(StaticData.Trajectories.sortedPoints(i).name, 'LHEE') == 1
            LHEEind = i;
        end
        if strcmp(StaticData.Trajectories.sortedPoints(i).name, 'RHEE') == 1
            RHEEind = i;
        end
        if strcmp(StaticData.Trajectories.sortedPoints(i).name, 'LTOE') == 1
            LTOEind = i;
        end
        if strcmp(StaticData.Trajectories.sortedPoints(i).name, 'RTOE') == 1
            RTOEind = i;
        end
    end
    
    LHEE = StaticData.Trajectories.sortedPoints(LHEEind).data(3,1,:);
    StaticData.LHEE = reshape(LHEE,[size(LHEE,3) 1]);
    LTOE = StaticData.Trajectories.sortedPoints(LTOEind).data(3,1,:);
    StaticData.LTOE = reshape(LTOE,[size(LTOE,3) 1]);
    RHEE = StaticData.Trajectories.sortedPoints(RHEEind).data(3,1,:);
    StaticData.RHEE = reshape(RHEE,[size(RHEE,3) 1]);
    RTOE = StaticData.Trajectories.sortedPoints(RTOEind).data(3,1,:);
    StaticData.RTOE = reshape(RTOE,[size(RTOE,3) 1]);
    
    % Plot to visualize marker heights over time
    % figure; hold on; grid on;
    % plot(StaticData.LHEE);
    % plot(StaticData.LTOE);
    % plot(StaticData.RHEE);
    % plot(StaticData.RTOE);
    % legend({'LHEE','LTOE','RHEE','RTOE'})
    
    % save average marker heights during cal trial and add 2 mm
    zThresh.LHEE = round((mean(StaticData.LHEE)) ./ 5)  + 3;
    zThresh.LTOE = round((mean(StaticData.LTOE)) ./ 5) + 3;
    zThresh.RHEE = round((mean(StaticData.RHEE)) ./ 5) + 3;
    zThresh.RTOE = round((mean(StaticData.RTOE)) ./ 5) + 3;
    
else % there are no cal trials
    zThresh.LHEE = 14; % if no cal trial, set height threshold to 14 mm for when in contact with PP mat
    zThresh.RHEE = 14;
    zThresh.LTOE = 14;
    zThresh.RTOE = 14;
end
xThresh = [-211 229]; % threshold settings determined from PP Global Calibration
yThresh = [2103 3548];

clearvars LHEE RHEE LTOE RTOE LHEEind RHEEind LTOEind RTOEind StaticTrial

%% Extract Marker Trajectories
NumTrials = length(C3Ddata);
for z = 1:NumTrials
    % find trajectories
    % LEFT
    LHEELog = strcmp(C3Ddata(z).Trajectories.PointLabels,'LHEE'); % heel marker
    LTOELog = strcmp(C3Ddata(z).Trajectories.PointLabels,'LTOE'); % toe
    LD1MLog = strcmp(C3Ddata(z).Trajectories.PointLabels,'LD1M'); % distal 1st met
    LD5MLog = strcmp(C3Ddata(z).Trajectories.PointLabels,'LD5M'); % distal 5th met
    LP1MLog = strcmp(C3Ddata(z).Trajectories.PointLabels,'LP1M'); % proximal 1st met
    LP5MLog = strcmp(C3Ddata(z).Trajectories.PointLabels,'LP5M'); % proximal 5th met
    LLCALog = strcmp(C3Ddata(z).Trajectories.PointLabels,'LLCA'); % lateral calcaneous
    LMCALog = strcmp(C3Ddata(z).Trajectories.PointLabels,'LMCA'); % medial calcaneous
    
    RHEELog = strcmp(C3Ddata(z).Trajectories.PointLabels,'RHEE'); % heel marker
    RTOELog = strcmp(C3Ddata(z).Trajectories.PointLabels,'RTOE'); % toe
    RD1MLog = strcmp(C3Ddata(z).Trajectories.PointLabels,'RD1M'); % distal 1st met
    RD5MLog = strcmp(C3Ddata(z).Trajectories.PointLabels,'RD5M'); % distal 5th met
    RP1MLog = strcmp(C3Ddata(z).Trajectories.PointLabels,'RP1M'); % proximal 1st met
    RP5MLog = strcmp(C3Ddata(z).Trajectories.PointLabels,'RP5M'); % proximal 5th met
    RLCALog = strcmp(C3Ddata(z).Trajectories.PointLabels,'RLCA'); % lateral calcaneous
    RMCALog = strcmp(C3Ddata(z).Trajectories.PointLabels,'RMCA'); % medial calcaneous
    
    [~,~,p] = size(C3Ddata(z).Trajectories.Point_xyz);
    for j = 1:p
        % LEFT
        C3Ddata(z).LHEE.raw(j,1:3) = C3Ddata(z).Trajectories.Point_xyz(1:3,LHEELog,j);
        C3Ddata(z).LTOE.raw(j,1:3) = C3Ddata(z).Trajectories.Point_xyz(1:3,LTOELog,j);
        C3Ddata(z).LD1M.raw(j,1:3) = C3Ddata(z).Trajectories.Point_xyz(1:3,LD1MLog,j);
        C3Ddata(z).LD5M.raw(j,1:3) = C3Ddata(z).Trajectories.Point_xyz(1:3,LD5MLog,j);
        C3Ddata(z).LP1M.raw(j,1:3) = C3Ddata(z).Trajectories.Point_xyz(1:3,LP1MLog,j);
        C3Ddata(z).LP5M.raw(j,1:3) = C3Ddata(z).Trajectories.Point_xyz(1:3,LP5MLog,j);
        C3Ddata(z).LLCA.raw(j,1:3) = C3Ddata(z).Trajectories.Point_xyz(1:3,LLCALog,j);
        C3Ddata(z).LMCA.raw(j,1:3) = C3Ddata(z).Trajectories.Point_xyz(1:3,LMCALog,j);
        % RIGHT
        C3Ddata(z).RHEE.raw(j,1:3) = C3Ddata(z).Trajectories.Point_xyz(1:3,RHEELog,j);
        C3Ddata(z).RTOE.raw(j,1:3) = C3Ddata(z).Trajectories.Point_xyz(1:3,RTOELog,j);
        C3Ddata(z).RD1M.raw(j,1:3) = C3Ddata(z).Trajectories.Point_xyz(1:3,RD1MLog,j);
        C3Ddata(z).RD5M.raw(j,1:3) = C3Ddata(z).Trajectories.Point_xyz(1:3,RD5MLog,j);
        C3Ddata(z).RP1M.raw(j,1:3) = C3Ddata(z).Trajectories.Point_xyz(1:3,RP1MLog,j);
        C3Ddata(z).RP5M.raw(j,1:3) = C3Ddata(z).Trajectories.Point_xyz(1:3,RP5MLog,j);
        C3Ddata(z).RLCA.raw(j,1:3) = C3Ddata(z).Trajectories.Point_xyz(1:3,RLCALog,j);
        C3Ddata(z).RMCA.raw(j,1:3) = C3Ddata(z).Trajectories.Point_xyz(1:3,RMCALog,j);
    end
    clearvars LHEELog LTOELog LD1MLog LD5MLog LP1MLog LP5MLog LLCALog LMCALog ...
        RHEELog RTOELog RD1MLog RD5MLog RP1MLog RP5MLog RLCALog RMCALog OnPlate
end


%% Align Coordinate Systems between PP and MoCap
Xoffset = 0;
Yoffset = PPOutput.Width;
for z = 1:NumTrials
    % LEFT
    % Extract marker times when they are within the plate range (X and Y values)
    % plate ends located at 2103 and 3548 mm from mocap origin, set 10 mm buffer for each side
    Xadj = 2103;
    Yadj = 220;
    OnPlate = find(C3Ddata(z).LHEE.raw(:,1) > 2093 & C3Ddata(z).LTOE.raw(:,1) < 3558); % Find LEFT indicies of foot within plate range
    C3Ddata(z).LHEE.Adj(:,1) = ((C3Ddata(z).LHEE.raw(OnPlate(1):OnPlate(end),1) - Xadj) ./ 5) + Xoffset;  % adjust X values by subtracting out edge of PP mat (2103) and then dividing by 5 b/c PP mat units are in 5mm sensors
    C3Ddata(z).LHEE.Adj(:,2) = ((C3Ddata(z).LHEE.raw(OnPlate(1):OnPlate(end),2) - Yadj) ./ 5 ) + Yoffset; % adjust Y values by subtracting out edge of PP mat (2103) and then dividing by 5 b/c PP mat units are in 5mm sensors, then add 88 cells back in to align correctly
    C3Ddata(z).LHEE.Adj(:,3) = C3Ddata(z).LHEE.raw(OnPlate(1):OnPlate(end),3) ./ 5; % divide z dimension by 5 to get all coordinates in the same units.
    C3Ddata(z).LTOE.Adj(:,1) = ((C3Ddata(z).LTOE.raw(OnPlate(1):OnPlate(end),1) - Xadj) ./ 5) + Xoffset; %
    C3Ddata(z).LTOE.Adj(:,2) = ((C3Ddata(z).LTOE.raw(OnPlate(1):OnPlate(end),2) - Yadj) ./ 5) + Yoffset;
    C3Ddata(z).LTOE.Adj(:,3) = C3Ddata(z).LTOE.raw(OnPlate(1):OnPlate(end),3) ./ 5;
    C3Ddata(z).LD1M.Adj(:,1) = ((C3Ddata(z).LD1M.raw(OnPlate(1):OnPlate(end),1) - Xadj) ./ 5) + Xoffset;%
    C3Ddata(z).LD1M.Adj(:,2) = ((C3Ddata(z).LD1M.raw(OnPlate(1):OnPlate(end),2) - Yadj) ./ 5) + Yoffset;
    C3Ddata(z).LD1M.Adj(:,3) = C3Ddata(z).LD1M.raw(OnPlate(1):OnPlate(end),3) ./ 5;
    C3Ddata(z).LD5M.Adj(:,1) = ((C3Ddata(z).LD5M.raw(OnPlate(1):OnPlate(end),1) - Xadj) ./ 5) + Xoffset;%
    C3Ddata(z).LD5M.Adj(:,2) = ((C3Ddata(z).LD5M.raw(OnPlate(1):OnPlate(end),2) - Yadj) ./ 5) + Yoffset;
    C3Ddata(z).LD5M.Adj(:,3) = C3Ddata(z).LD5M.raw(OnPlate(1):OnPlate(end),3) ./ 5;
    C3Ddata(z).LP1M.Adj(:,1) = ((C3Ddata(z).LP1M.raw(OnPlate(1):OnPlate(end),1) - Xadj) ./ 5) + Xoffset;%
    C3Ddata(z).LP1M.Adj(:,2) = ((C3Ddata(z).LP1M.raw(OnPlate(1):OnPlate(end),2) - Yadj) ./ 5) + Yoffset;
    C3Ddata(z).LP1M.Adj(:,3) = C3Ddata(z).LP1M.raw(OnPlate(1):OnPlate(end),3) ./ 5;
    C3Ddata(z).LP5M.Adj(:,1) = ((C3Ddata(z).LP5M.raw(OnPlate(1):OnPlate(end),1) - Xadj) ./ 5) + Xoffset;%
    C3Ddata(z).LP5M.Adj(:,2) = ((C3Ddata(z).LP5M.raw(OnPlate(1):OnPlate(end),2) - Yadj) ./ 5) + Yoffset;
    C3Ddata(z).LP5M.Adj(:,3) = C3Ddata(z).LP5M.raw(OnPlate(1):OnPlate(end),3) ./ 5;
    C3Ddata(z).LLCA.Adj(:,1) = ((C3Ddata(z).LLCA.raw(OnPlate(1):OnPlate(end),1) - Xadj) ./ 5) + Xoffset;%
    C3Ddata(z).LLCA.Adj(:,2) = ((C3Ddata(z).LLCA.raw(OnPlate(1):OnPlate(end),2) - Yadj) ./ 5) + Yoffset;
    C3Ddata(z).LLCA.Adj(:,3) = C3Ddata(z).LLCA.raw(OnPlate(1):OnPlate(end),3) ./ 5;
    C3Ddata(z).LMCA.Adj(:,1) = ((C3Ddata(z).LMCA.raw(OnPlate(1):OnPlate(end),1) - Xadj) ./ 5) + Xoffset;%
    C3Ddata(z).LMCA.Adj(:,2) = ((C3Ddata(z).LMCA.raw(OnPlate(1):OnPlate(end),2) - Yadj) ./ 5) + Yoffset;
    C3Ddata(z).LMCA.Adj(:,3) = C3Ddata(z).LMCA.raw(OnPlate(1):OnPlate(end),3) ./ 5;
    
    % flip Y values to match PP coordinates
    C3Ddata(z).LHEE.Adj(:,2) = PPOutput.Width - C3Ddata(z).LHEE.Adj(:,2);
    C3Ddata(z).LTOE.Adj(:,2) = PPOutput.Width - C3Ddata(z).LTOE.Adj(:,2);
    C3Ddata(z).LD1M.Adj(:,2) = PPOutput.Width - C3Ddata(z).LD1M.Adj(:,2);
    C3Ddata(z).LD5M.Adj(:,2) = PPOutput.Width - C3Ddata(z).LD5M.Adj(:,2);
    C3Ddata(z).LP1M.Adj(:,2) = PPOutput.Width - C3Ddata(z).LP1M.Adj(:,2);
    C3Ddata(z).LP5M.Adj(:,2) = PPOutput.Width - C3Ddata(z).LP5M.Adj(:,2);
    C3Ddata(z).LLCA.Adj(:,2) = PPOutput.Width - C3Ddata(z).LLCA.Adj(:,2);
    C3Ddata(z).LMCA.Adj(:,2) = PPOutput.Width - C3Ddata(z).LMCA.Adj(:,2);
    
    Trim(z).Left = OnPlate;
    clearvars OnPlate Press
    % compute trajectory velocity for heel and toe markers
    C3Ddata(z).LHEE.vel(1) = NaN;
    C3Ddata(z).LTOE.vel(1) = NaN;
    for j = 2:length(C3Ddata(z).LHEE.Adj)
        C3Ddata(z).LHEE.vel(j) = sqrt(sum(C3Ddata(z).LHEE.Adj(j,:) - C3Ddata(z).LHEE.Adj(j-1,:)).^2);
        C3Ddata(z).LTOE.vel(j) = sqrt(sum(C3Ddata(z).LTOE.Adj(j,:) - C3Ddata(z).LTOE.Adj(j-1,:)).^2);
    end
    % filter velocities with a 3-point moving mean
    v = datenum(version('-date'));
    if v < datenum('June 1, 2018')
        p = length(C3Ddata(z).LHEE.vel);
        C3Ddata(z).LHEE.velF(1) = C3Ddata(z).LHEE.vel(1);
        C3Ddata(z).LHEE.velF(p) = C3Ddata(z).LHEE.vel(p);
        for i = 2:p - 1
            C3Ddata(z).LHEE.velF(i) = mean(C3Ddata(z).LHEE.vel(i-1:i+1));
        end
        p = length(C3Ddata(z).LTOE.vel);
        C3Ddata(z).LTOE.velF(1) = C3Ddata(z).LTOE.vel(1);
        C3Ddata(z).LTOE.velF(p) = C3Ddata(z).LTOE.vel(p);
        for i = 2:p - 1
            C3Ddata(z).LTOE.velF(i) = mean(C3Ddata(z).LTOE.vel(i-1:i+1));
        end
    else
        C3Ddata(z).LHEE.velF = movmean(C3Ddata(z).LHEE.vel,3);
        C3Ddata(z).LTOE.velF = movmean(C3Ddata(z).LTOE.vel,3);
    end
    %     figure;
    %     plot(C3Ddata(z).LHEE.vel); hold on;
    %     plot(C3Ddata(z).LTOE.vel);
    %     plot(C3Ddata(z).LHEE.velF);
    %     plot(C3Ddata(z).LTOE.velF);
    
    % RIGHT
    % plate ends located at 2103 and 3548 mm from mocap origin, set 10 mm buffer for each side
    OnPlate = find(C3Ddata(z).RHEE.raw(:,1) > 2090 & C3Ddata(z).RTOE.raw(:,1) < 3560); % Find RIGHT indicies of foot within plate range
    C3Ddata(z).RHEE.Adj(:,1) = ((C3Ddata(z).RHEE.raw(OnPlate(1):OnPlate(end),1) - Xadj) ./ 5) + Xoffset;
    C3Ddata(z).RHEE.Adj(:,2) = ((C3Ddata(z).RHEE.raw(OnPlate(1):OnPlate(end),2) - Yadj) ./ 5) + Yoffset;
    C3Ddata(z).RHEE.Adj(:,3) = C3Ddata(z).RHEE.raw(OnPlate(1):OnPlate(end),3) ./ 5;
    C3Ddata(z).RTOE.Adj(:,1) = ((C3Ddata(z).RTOE.raw(OnPlate(1):OnPlate(end),1) - Xadj) ./ 5) + Xoffset;
    C3Ddata(z).RTOE.Adj(:,2) = ((C3Ddata(z).RTOE.raw(OnPlate(1):OnPlate(end),2) - Yadj) ./ 5) + Yoffset;
    C3Ddata(z).RTOE.Adj(:,3) = C3Ddata(z).RTOE.raw(OnPlate(1):OnPlate(end),3) ./ 5;
    C3Ddata(z).RD1M.Adj(:,1) = ((C3Ddata(z).RD1M.raw(OnPlate(1):OnPlate(end),1) - Xadj) ./ 5) + Xoffset;%
    C3Ddata(z).RD1M.Adj(:,2) = ((C3Ddata(z).RD1M.raw(OnPlate(1):OnPlate(end),2) - Yadj) ./ 5) + Yoffset;
    C3Ddata(z).RD1M.Adj(:,3) = C3Ddata(z).RD1M.raw(OnPlate(1):OnPlate(end),3) ./ 5;
    C3Ddata(z).RD5M.Adj(:,1) = ((C3Ddata(z).RD5M.raw(OnPlate(1):OnPlate(end),1) - Xadj) ./ 5) + Xoffset;%
    C3Ddata(z).RD5M.Adj(:,2) = ((C3Ddata(z).RD5M.raw(OnPlate(1):OnPlate(end),2) - Yadj) ./ 5) + Yoffset;
    C3Ddata(z).RD5M.Adj(:,3) = C3Ddata(z).RD5M.raw(OnPlate(1):OnPlate(end),3) ./ 5;
    C3Ddata(z).RP1M.Adj(:,1) = ((C3Ddata(z).RP1M.raw(OnPlate(1):OnPlate(end),1) - Xadj) ./ 5) + Xoffset;%
    C3Ddata(z).RP1M.Adj(:,2) = ((C3Ddata(z).RP1M.raw(OnPlate(1):OnPlate(end),2) - Yadj) ./ 5) + Yoffset;
    C3Ddata(z).RP1M.Adj(:,3) = C3Ddata(z).RP1M.raw(OnPlate(1):OnPlate(end),3) ./ 5;
    C3Ddata(z).RP5M.Adj(:,1) = ((C3Ddata(z).RP5M.raw(OnPlate(1):OnPlate(end),1) - Xadj) ./ 5) + Xoffset;%
    C3Ddata(z).RP5M.Adj(:,2) = ((C3Ddata(z).RP5M.raw(OnPlate(1):OnPlate(end),2) - Yadj) ./ 5) + Yoffset;
    C3Ddata(z).RP5M.Adj(:,3) = C3Ddata(z).RP5M.raw(OnPlate(1):OnPlate(end),3) ./ 5;
    C3Ddata(z).RLCA.Adj(:,1) = ((C3Ddata(z).RLCA.raw(OnPlate(1):OnPlate(end),1) - Xadj) ./ 5) + Xoffset;%
    C3Ddata(z).RLCA.Adj(:,2) = ((C3Ddata(z).RLCA.raw(OnPlate(1):OnPlate(end),2) - Yadj) ./ 5) + Yoffset;
    C3Ddata(z).RLCA.Adj(:,3) = C3Ddata(z).RLCA.raw(OnPlate(1):OnPlate(end),3) ./ 5;
    C3Ddata(z).RMCA.Adj(:,1) = ((C3Ddata(z).RMCA.raw(OnPlate(1):OnPlate(end),1) - Xadj) ./ 5) + Xoffset;%
    C3Ddata(z).RMCA.Adj(:,2) = ((C3Ddata(z).RMCA.raw(OnPlate(1):OnPlate(end),2) - Yadj) ./ 5) + Yoffset;
    C3Ddata(z).RMCA.Adj(:,3) = C3Ddata(z).RMCA.raw(OnPlate(1):OnPlate(end),3) ./ 5;
    
    % flip Y values to match PP coordinates
    C3Ddata(z).RHEE.Adj(:,2) = PPOutput.Width - C3Ddata(z).RHEE.Adj(:,2);
    C3Ddata(z).RTOE.Adj(:,2) = PPOutput.Width - C3Ddata(z).RTOE.Adj(:,2);
    C3Ddata(z).RD1M.Adj(:,2) = PPOutput.Width - C3Ddata(z).RD1M.Adj(:,2);
    C3Ddata(z).RD5M.Adj(:,2) = PPOutput.Width - C3Ddata(z).RD5M.Adj(:,2);
    C3Ddata(z).RP1M.Adj(:,2) = PPOutput.Width - C3Ddata(z).RP1M.Adj(:,2);
    C3Ddata(z).RP5M.Adj(:,2) = PPOutput.Width - C3Ddata(z).RP5M.Adj(:,2);
    C3Ddata(z).RLCA.Adj(:,2) = PPOutput.Width - C3Ddata(z).RLCA.Adj(:,2);
    C3Ddata(z).RMCA.Adj(:,2) = PPOutput.Width - C3Ddata(z).RMCA.Adj(:,2);
    
    Trim(z).Right = OnPlate;
    clearvars OnPlate
    % compute trajectory velocity for heel and toe markers
    C3Ddata(z).RHEE.vel(1) = NaN;
    C3Ddata(z).RTOE.vel(1) = NaN;
    for j = 2:length(C3Ddata(z).RHEE.Adj)
        C3Ddata(z).RHEE.vel(j) = sqrt(sum(C3Ddata(z).RHEE.Adj(j,:) - C3Ddata(z).RHEE.Adj(j-1,:)).^2);
        C3Ddata(z).RTOE.vel(j) = sqrt(sum(C3Ddata(z).RTOE.Adj(j,:) - C3Ddata(z).RTOE.Adj(j-1,:)).^2);
    end
    % filter velocities with a 3-point moving mean
    % filter velocities with a 3-point moving mean
    v = datenum(version('-date'));
    if v < datenum('June 1, 2018')
        p = length(C3Ddata(z).RHEE.vel);
        C3Ddata(z).RHEE.velF(1) = C3Ddata(z).RHEE.vel(1);
        C3Ddata(z).RHEE.velF(p) = C3Ddata(z).RHEE.vel(p);
        for i = 2:p - 1
            C3Ddata(z).RHEE.velF(i) = mean(C3Ddata(z).RHEE.vel(i-1:i+1));
        end
        p = length(C3Ddata(z).RTOE.vel);
        C3Ddata(z).RTOE.velF(1) = C3Ddata(z).RTOE.vel(1);
        C3Ddata(z).RTOE.velF(p) = C3Ddata(z).RTOE.vel(p);
        for i = 2:p - 1
            C3Ddata(z).RTOE.velF(i) = mean(C3Ddata(z).RTOE.vel(i-1:i+1));
        end
    else
        C3Ddata(z).RHEE.velF = movmean(C3Ddata(z).RHEE.vel,3);
        C3Ddata(z).RTOE.velF = movmean(C3Ddata(z).RTOE.vel,3);
    end
    
    if strcmp(C3Ddata(z).Direction, 'South') == 1 % if they are walking south, flip marker data;
        L = PPOutput.Length;
        W = PPOutput.Width;
        % LEFT
        C3Ddata(z).LHEE.Adj(:,1) = L - C3Ddata(z).LHEE.Adj(:,1) ;  % flip X values by subtracting from the length of the mat
        C3Ddata(z).LHEE.Adj(:,2) = W - C3Ddata(z).LHEE.Adj(:,2); % flip Y values by subtracting from the width of the mat
        C3Ddata(z).LTOE.Adj(:,1) = L - C3Ddata(z).LTOE.Adj(:,1);
        C3Ddata(z).LTOE.Adj(:,2) = W - C3Ddata(z).LTOE.Adj(:,2);
        C3Ddata(z).LD1M.Adj(:,1) = L - C3Ddata(z).LD1M.Adj(:,1);
        C3Ddata(z).LD1M.Adj(:,2) = W - C3Ddata(z).LD1M.Adj(:,2);
        C3Ddata(z).LD5M.Adj(:,1) = L - C3Ddata(z).LD5M.Adj(:,1);
        C3Ddata(z).LD5M.Adj(:,2) = W - C3Ddata(z).LD5M.Adj(:,2);
        C3Ddata(z).LP1M.Adj(:,1) = L - C3Ddata(z).LP1M.Adj(:,1);
        C3Ddata(z).LP1M.Adj(:,2) = W - C3Ddata(z).LP1M.Adj(:,2);
        C3Ddata(z).LP5M.Adj(:,1) = L - C3Ddata(z).LP5M.Adj(:,1);
        C3Ddata(z).LP5M.Adj(:,2) = W - C3Ddata(z).LP5M.Adj(:,2);
        C3Ddata(z).LLCA.Adj(:,1) = L - C3Ddata(z).LLCA.Adj(:,1);
        C3Ddata(z).LLCA.Adj(:,2) = W - C3Ddata(z).LLCA.Adj(:,2);
        C3Ddata(z).LMCA.Adj(:,1) = L - C3Ddata(z).LMCA.Adj(:,1);
        C3Ddata(z).LMCA.Adj(:,2) = W - C3Ddata(z).LMCA.Adj(:,2);
        % RIGHT
        C3Ddata(z).RHEE.Adj(:,1) = L - C3Ddata(z).RHEE.Adj(:,1) ;  % flip X values by subtracting from the length of the mat
        C3Ddata(z).RHEE.Adj(:,2) = W - C3Ddata(z).RHEE.Adj(:,2); % flip Y values by subtracting from the width of the mat
        C3Ddata(z).RTOE.Adj(:,1) = L - C3Ddata(z).RTOE.Adj(:,1);
        C3Ddata(z).RTOE.Adj(:,2) = W - C3Ddata(z).RTOE.Adj(:,2);
        C3Ddata(z).RD1M.Adj(:,1) = L - C3Ddata(z).RD1M.Adj(:,1);
        C3Ddata(z).RD1M.Adj(:,2) = W - C3Ddata(z).RD1M.Adj(:,2);
        C3Ddata(z).RD5M.Adj(:,1) = L - C3Ddata(z).RD5M.Adj(:,1);
        C3Ddata(z).RD5M.Adj(:,2) = W - C3Ddata(z).RD5M.Adj(:,2);
        C3Ddata(z).RP1M.Adj(:,1) = L - C3Ddata(z).RP1M.Adj(:,1);
        C3Ddata(z).RP1M.Adj(:,2) = W - C3Ddata(z).RP1M.Adj(:,2);
        C3Ddata(z).RP5M.Adj(:,1) = L - C3Ddata(z).RP5M.Adj(:,1);
        C3Ddata(z).RP5M.Adj(:,2) = W - C3Ddata(z).RP5M.Adj(:,2);
        C3Ddata(z).RLCA.Adj(:,1) = L - C3Ddata(z).RLCA.Adj(:,1);
        C3Ddata(z).RLCA.Adj(:,2) = W - C3Ddata(z).RLCA.Adj(:,2);
        C3Ddata(z).RMCA.Adj(:,1) = L - C3Ddata(z).RMCA.Adj(:,1);
        C3Ddata(z).RMCA.Adj(:,2) = W - C3Ddata(z).RMCA.Adj(:,2);
    else
        L = 0;
        W = 4; % offset between mat and lab
         % LEFT
        C3Ddata(z).LHEE.Adj(:,1) = L + C3Ddata(z).LHEE.Adj(:,1) ;  % flip X values by subtracting from the length of the mat
        C3Ddata(z).LHEE.Adj(:,2) = W + C3Ddata(z).LHEE.Adj(:,2); % flip Y values by subtracting from the width of the mat
        C3Ddata(z).LTOE.Adj(:,1) = L + C3Ddata(z).LTOE.Adj(:,1);
        C3Ddata(z).LTOE.Adj(:,2) = W + C3Ddata(z).LTOE.Adj(:,2);
        C3Ddata(z).LD1M.Adj(:,1) = L + C3Ddata(z).LD1M.Adj(:,1);
        C3Ddata(z).LD1M.Adj(:,2) = W + C3Ddata(z).LD1M.Adj(:,2);
        C3Ddata(z).LD5M.Adj(:,1) = L + C3Ddata(z).LD5M.Adj(:,1);
        C3Ddata(z).LD5M.Adj(:,2) = W + C3Ddata(z).LD5M.Adj(:,2);
        C3Ddata(z).LP1M.Adj(:,1) = L + C3Ddata(z).LP1M.Adj(:,1);
        C3Ddata(z).LP1M.Adj(:,2) = W + C3Ddata(z).LP1M.Adj(:,2);
        C3Ddata(z).LP5M.Adj(:,1) = L + C3Ddata(z).LP5M.Adj(:,1);
        C3Ddata(z).LP5M.Adj(:,2) = W + C3Ddata(z).LP5M.Adj(:,2);
        C3Ddata(z).LLCA.Adj(:,1) = L + C3Ddata(z).LLCA.Adj(:,1);
        C3Ddata(z).LLCA.Adj(:,2) = W + C3Ddata(z).LLCA.Adj(:,2);
        C3Ddata(z).LMCA.Adj(:,1) = L + C3Ddata(z).LMCA.Adj(:,1);
        C3Ddata(z).LMCA.Adj(:,2) = W + C3Ddata(z).LMCA.Adj(:,2);
        % RIGHT
        C3Ddata(z).RHEE.Adj(:,1) = L + C3Ddata(z).RHEE.Adj(:,1) ;  % flip X values by subtracting from the length of the mat
        C3Ddata(z).RHEE.Adj(:,2) = W + C3Ddata(z).RHEE.Adj(:,2); % flip Y values by subtracting from the width of the mat
        C3Ddata(z).RTOE.Adj(:,1) = L + C3Ddata(z).RTOE.Adj(:,1);
        C3Ddata(z).RTOE.Adj(:,2) = W + C3Ddata(z).RTOE.Adj(:,2);
        C3Ddata(z).RD1M.Adj(:,1) = L + C3Ddata(z).RD1M.Adj(:,1);
        C3Ddata(z).RD1M.Adj(:,2) = W + C3Ddata(z).RD1M.Adj(:,2);
        C3Ddata(z).RD5M.Adj(:,1) = L + C3Ddata(z).RD5M.Adj(:,1);
        C3Ddata(z).RD5M.Adj(:,2) = W + C3Ddata(z).RD5M.Adj(:,2);
        C3Ddata(z).RP1M.Adj(:,1) = L + C3Ddata(z).RP1M.Adj(:,1);
        C3Ddata(z).RP1M.Adj(:,2) = W + C3Ddata(z).RP1M.Adj(:,2);
        C3Ddata(z).RP5M.Adj(:,1) = L + C3Ddata(z).RP5M.Adj(:,1);
        C3Ddata(z).RP5M.Adj(:,2) = W + C3Ddata(z).RP5M.Adj(:,2);
        C3Ddata(z).RLCA.Adj(:,1) = L + C3Ddata(z).RLCA.Adj(:,1);
        C3Ddata(z).RLCA.Adj(:,2) = W + C3Ddata(z).RLCA.Adj(:,2);
        C3Ddata(z).RMCA.Adj(:,1) = L + C3Ddata(z).RMCA.Adj(:,1);
        C3Ddata(z).RMCA.Adj(:,2) = W + C3Ddata(z).RMCA.Adj(:,2);
    end
end

%% determine number of steps on mat
% set velocity threshold for foot contace
velThresh = 0.3;
stanceThresh = 8;
for z = 1:NumTrials
    % LEFT
    % classify stance periods
    for j = 1:length(C3Ddata(z).LHEE.Adj(:,3))
        if C3Ddata(z).LHEE.Adj(j,3) < round(zThresh.LHEE) == 1 || C3Ddata(z).LTOE.Adj(j,3) < round(zThresh.LTOE) == 1
            if C3Ddata(z).LTOE.velF(j) < velThresh == 1 || C3Ddata(z).LHEE.velF(j) < velThresh == 1
                C3Ddata(z).LeftContact(1).All(j) = 1;
            else
                C3Ddata(z).LeftContact(1).All(j) = 0;
            end
        else
            C3Ddata(z).LeftContact(1).All(j) = 0;
        end
    end
    % quality assurance on stance periods
    for j  = 1:length(C3Ddata(z).LHEE.Adj(:,3))
        if C3Ddata(z).LeftContact(1).All(j) == 1 && C3Ddata(z).LeftContact(1).All(j-1) == 0 % if theres a random 1, surrounded by 0s, and no more 1s within 5 frames, delete
            p = find(C3Ddata(z).LeftContact(1).All(j+1:end)==0,1);
            if p < stanceThresh
                C3Ddata(z).LeftContact(1).All(j:j+p-1) = 0;
            end
        end
    end
    % determine number of steps
    [~,locs] =   findpeaks(double(C3Ddata(z).LeftContact(1).All), 'MinPeakHeight',0.5);
    for i = 1:length(locs)
        NextEnd = find(C3Ddata(z).LeftContact(1).All(locs(i):end) == 0, 1);
        % if NextEnd doesnt exist, set to end of C3Ddata(z).LeftContact.All
        C3Ddata(z).LeftContact(1).All(locs(i):NextEnd+locs(i)-2) = i;
        % save stance periods into structure
        C3Ddata(z).LeftContact(i).On = locs(i);
        C3Ddata(z).LeftContact(i).Off = NextEnd + locs(i) - 2;
        C3Ddata(z).LeftContact(i).Mid = round(mean([C3Ddata(z).LeftContact(i).On, C3Ddata(z).LeftContact(i).Off]));
    end
    if isempty(locs)
        C3Ddata(z).LeftContact(1).Mid = NaN;
    end
    % RIGHT
    % classify stance periods
    for j = 1:length(C3Ddata(z).RHEE.Adj(:,3))
        if C3Ddata(z).RHEE.Adj(j,3) < round(zThresh.RHEE) == 1 || C3Ddata(z).RTOE.Adj(j,3) < round(zThresh.RTOE) == 1
            if C3Ddata(z).RTOE.velF(j) < velThresh == 1 || C3Ddata(z).RHEE.velF(j) < velThresh == 1
                C3Ddata(z).RightContact(1).All(j) = 1;
            else
                C3Ddata(z).RightContact(1).All(j) = 0;
            end
        else
            C3Ddata(z).RightContact(1).All(j) = 0;
        end
    end
    % quality assurance on stance periods
    for j  = 1:length(C3Ddata(z).RHEE.Adj(:,3))
        if C3Ddata(z).RightContact(1).All(j) == 1 && C3Ddata(z).RightContact(1).All(j-1) == 0 % if theres a random 1, surrounded by 0s, and no more 1s within 5 frames, delete
            p = find(C3Ddata(z).RightContact(1).All(j+1:end)==0,1);
            if p < stanceThresh
                C3Ddata(z).RightContact(1).All(j:j+p-1) = 0;
            end
        end
    end
    % determine number of steps
    [~,locs] = findpeaks(double(C3Ddata(z).RightContact(1).All), 'MinPeakHeight',0.5);
    for i = 1:length(locs)
        NextEnd = find(C3Ddata(z).RightContact(1).All(locs(i):end) == 0, 1);
        % if NextEnd doesnt exist, set to end of C3Ddata(z).LeftContact.All
        C3Ddata(z).RightContact(1).All(locs(i):NextEnd+locs(i)-2) = i;
        % save stance periods into structure
        C3Ddata(z).RightContact(i).On = locs(i);
        C3Ddata(z).RightContact(i).Off = NextEnd + locs(i) - 2;
        C3Ddata(z).RightContact(i).Mid = round(mean([C3Ddata(z).RightContact(i).On, C3Ddata(z).RightContact(i).Off]));
    end
    if isempty(locs)
        C3Ddata(z).RightContact(1).Mid = NaN;
    end
end
clearvars stanceThresh velThresh i j q z p L W Yoffset Xoffset

%% define foot progression angles and unit vectors
for z = 1:NumTrials
    % LEFT
    for i = 1:length(C3Ddata(z).LeftContact)
        if isnan(C3Ddata(z).LeftContact(i).Mid) == 0
            p = C3Ddata(z).LeftContact(i).Mid;
        else
            break
        end
        L =  sqrt(abs(C3Ddata(z).LTOE.Adj(p,1) - C3Ddata(z).LHEE.Adj(p,1)).^2 + abs(C3Ddata(z).LTOE.Adj(p,2) - C3Ddata(z).LHEE.Adj(p,2)).^2);
        C3Ddata(z).LeftContact(i).UnitVector = [(C3Ddata(z).LTOE.Adj(p,1) - C3Ddata(z).LHEE.Adj(p,1)) / L, (C3Ddata(z).LTOE.Adj(p,2) - C3Ddata(z).LHEE.Adj(p,2)) / L];
        clearvars L p
    end
    % RIGHT
    for i = 1:length(C3Ddata(z).RightContact)
        if isnan(C3Ddata(z).RightContact(i).Mid) == 0
            p = C3Ddata(z).RightContact(i).Mid;
        else
            break
        end
        L =  sqrt(abs(C3Ddata(z).RTOE.Adj(p,1) - C3Ddata(z).RHEE.Adj(p,1)).^2 + abs(C3Ddata(z).RTOE.Adj(p,2) - C3Ddata(z).RHEE.Adj(p,2)).^2);
        C3Ddata(z).RightContact(i).UnitVector = [(C3Ddata(z).RTOE.Adj(p,1) - C3Ddata(z).RHEE.Adj(p,1)) / L, (C3Ddata(z).RTOE.Adj(p,2) - C3Ddata(z).RHEE.Adj(p,2)) / L];
        clearvars L p
    end
end

%% compute actual heel point using marker offset, and foot length
MarkerBase = 0.8 /2 ; % 0.8 cm offset between marker base and centroid - converted to half-cm units to match pp mat
% FLScale = 0.95; % Scaling factor to adjust for fact that whole foot does not contact the ground
% Subject.FootLength = Subject.FootLength .* FLScale;
% LEFT
for z = 1:NumTrials
    % LEFT
    for i = 1: length(C3Ddata(z).LeftContact)
        if isnan(C3Ddata(z).LeftContact(i).Mid)
            M = find(C3Ddata(z).LeftContact(i).All,1);
            N =  length(C3Ddata(z).LeftContact.All) + 1 - find(fliplr(C3Ddata(z).LeftContact(i).All),1);
            p = round(mean([M N]));
            C3Ddata(z).LeftContact(i).Mid = p;
        else
            p = C3Ddata(z).LeftContact(i).Mid;
        end
        
        C3Ddata(z).LHEE.base(i,1) = C3Ddata(z).LHEE.Adj(p,1) +  (MarkerBase .* C3Ddata(z).LeftContact(i).UnitVector(1)); % creating new heel marker at base of heel marker
        C3Ddata(z).LHEE.base(i,2) = C3Ddata(z).LHEE.Adj(p,2) +  (MarkerBase .* C3Ddata(z).LeftContact(i).UnitVector(2));
        C3Ddata(z).LTOE.top(i,1) = C3Ddata(z).LHEE.Adj(p,1) +  ((Subject.FootLength .* 2) .* C3Ddata(z).LeftContact(i).UnitVector(1)); % put in length of foot for total foot length
        C3Ddata(z).LTOE.top(i,2) = C3Ddata(z).LHEE.Adj(p,2) +  ((Subject.FootLength .* 2) .* C3Ddata(z).LeftContact(i).UnitVector(2));
        
        C3Ddata(z).LTOE.step(i,1) = C3Ddata(z).LTOE.Adj(p,1);
        C3Ddata(z).LTOE.step(i,2) = C3Ddata(z).LTOE.Adj(p,2);
        C3Ddata(z).LD1M.step(i,1) = C3Ddata(z).LD1M.Adj(p,1);
        C3Ddata(z).LD1M.step(i,2) = C3Ddata(z).LD1M.Adj(p,2);
        C3Ddata(z).LD5M.step(i,1) = C3Ddata(z).LD5M.Adj(p,1);
        C3Ddata(z).LD5M.step(i,2) = C3Ddata(z).LD5M.Adj(p,2);
        C3Ddata(z).LP1M.step(i,1) = C3Ddata(z).LP1M.Adj(p,1);
        C3Ddata(z).LP1M.step(i,2) = C3Ddata(z).LP1M.Adj(p,2);
        C3Ddata(z).LP5M.step(i,1) = C3Ddata(z).LP5M.Adj(p,1);
        C3Ddata(z).LP5M.step(i,2) = C3Ddata(z).LP5M.Adj(p,2);
        C3Ddata(z).LLCA.step(i,1) = C3Ddata(z).LLCA.Adj(p,1);
        C3Ddata(z).LLCA.step(i,2) = C3Ddata(z).LLCA.Adj(p,2);
        C3Ddata(z).LMCA.step(i,1) = C3Ddata(z).LMCA.Adj(p,1);
        C3Ddata(z).LMCA.step(i,2) = C3Ddata(z).LMCA.Adj(p,2);
        
        C3Ddata(z).LHA(i,1) = C3Ddata(z).LHEE.Adj(p,1) +  (Subject.FootLength .* 2 .* 0.3 .* C3Ddata(z).LeftContact(i).UnitVector(1)); % 1/3rd line - HA = Heel/Arch
        C3Ddata(z).LHA(i,2) = C3Ddata(z).LHEE.Adj(p,2) +  (Subject.FootLength .* 2 .* 0.3 .* C3Ddata(z).LeftContact(i).UnitVector(2));
        C3Ddata(z).LFA(i,1) = C3Ddata(z).LHEE.Adj(p,1) +  (Subject.FootLength .* 2 .* 0.6 .* C3Ddata(z).LeftContact(i).UnitVector(1));  % 2/3rd line - FA = Fore/Arch
        C3Ddata(z).LFA(i,2) = C3Ddata(z).LHEE.Adj(p,2) +  (Subject.FootLength .* 2 .* 0.6 .* C3Ddata(z).LeftContact(i).UnitVector(2));
        
        C3Ddata(z).LHAb(i,1) = C3Ddata(z).LHEE.base(i,1) +  (Subject.FootLength .* 2 .* 0.3 .* C3Ddata(z).LeftContact(i).UnitVector(1)); % 1/3rd line - HA = Heel/Arch
        C3Ddata(z).LHAb(i,2) = C3Ddata(z).LHEE.base(i,2) +  (Subject.FootLength .* 2 .* 0.3 .* C3Ddata(z).LeftContact(i).UnitVector(2));
        C3Ddata(z).LFAb(i,1) = C3Ddata(z).LHEE.base(i,1) +  (Subject.FootLength .* 2 .* 0.6 .* C3Ddata(z).LeftContact(i).UnitVector(1));  % 2/3rd line - FA = Fore/Arch
        C3Ddata(z).LFAb(i,2) = C3Ddata(z).LHEE.base(i,2) +  (Subject.FootLength .* 2 .* 0.6 .* C3Ddata(z).LeftContact(i).UnitVector(2));
        
        C3Ddata(z).LCent(i,1) = C3Ddata(z).LHEE.Adj(p,1) +  (Subject.FootLength .* 2 .* 0.5.* C3Ddata(z).LeftContact(i).UnitVector(1)); % centroid from adj points
        C3Ddata(z).LCent(i,2) = C3Ddata(z).LHEE.Adj(p,2) +  (Subject.FootLength .* 2 .* 0.5 .* C3Ddata(z).LeftContact(i).UnitVector(2));
        C3Ddata(z).LCentB(i,1) = C3Ddata(z).LHEE.base(i,1) +  (Subject.FootLength .* 2 .* 0.5 .* C3Ddata(z).LeftContact(i).UnitVector(1)); % centroid from base points
        C3Ddata(z).LCentB(i,2) = C3Ddata(z).LHEE.base(i,2) +  (Subject.FootLength .* 2 .* 0.5 .* C3Ddata(z).LeftContact(i).UnitVector(2));
        
    end
    % RIGHT
    for i = 1: length(C3Ddata(z).RightContact)
        if isnan(C3Ddata(z).RightContact(i).Mid)
            break
        else
            p = C3Ddata(z).RightContact(i).Mid;
        end
        C3Ddata(z).RHEE.base(i,1) = C3Ddata(z).RHEE.Adj(p,1) +  (MarkerBase .* C3Ddata(z).RightContact(i).UnitVector(1)); % creating new heel marker at base of heel marker
        C3Ddata(z).RHEE.base(i,2) = C3Ddata(z).RHEE.Adj(p,2) +  (MarkerBase .* C3Ddata(z).RightContact(i).UnitVector(2));
        C3Ddata(z).RTOE.top(i,1) = C3Ddata(z).RHEE.base(i,1) +  ((Subject.FootLength .* 2) .* C3Ddata(z).RightContact(i).UnitVector(1)); % put in length of foot for total foot length
        C3Ddata(z).RTOE.top(i,2) = C3Ddata(z).RHEE.base(i,2) +  ((Subject.FootLength .* 2) .* C3Ddata(z).RightContact(i).UnitVector(2));
        
        C3Ddata(z).RTOE.step(i,1) = C3Ddata(z).RTOE.Adj(p,1);
        C3Ddata(z).RTOE.step(i,2) = C3Ddata(z).RTOE.Adj(p,2);
        C3Ddata(z).RD1M.step(i,1) = C3Ddata(z).RD1M.Adj(p,1);
        C3Ddata(z).RD1M.step(i,2) = C3Ddata(z).RD1M.Adj(p,2);
        C3Ddata(z).RD5M.step(i,1) = C3Ddata(z).RD5M.Adj(p,1);
        C3Ddata(z).RD5M.step(i,2) = C3Ddata(z).RD5M.Adj(p,2);
        C3Ddata(z).RP1M.step(i,1) = C3Ddata(z).RP1M.Adj(p,1);
        C3Ddata(z).RP1M.step(i,2) = C3Ddata(z).RP1M.Adj(p,2);
        C3Ddata(z).RP5M.step(i,1) = C3Ddata(z).RP5M.Adj(p,1);
        C3Ddata(z).RP5M.step(i,2) = C3Ddata(z).RP5M.Adj(p,2);
        C3Ddata(z).RLCA.step(i,1) = C3Ddata(z).RLCA.Adj(p,1);
        C3Ddata(z).RLCA.step(i,2) = C3Ddata(z).RLCA.Adj(p,2);
        C3Ddata(z).RMCA.step(i,1) = C3Ddata(z).RMCA.Adj(p,1);
        C3Ddata(z).RMCA.step(i,2) = C3Ddata(z).RMCA.Adj(p,2);
        
        C3Ddata(z).RHA(i,1) = C3Ddata(z).RHEE.Adj(p,1) +  (Subject.FootLength .* 2 .* 0.3 .* C3Ddata(z).RightContact(i).UnitVector(1)); % 1/3rd line - HA = Heel/Arch
        C3Ddata(z).RHA(i,2) = C3Ddata(z).RHEE.Adj(p,2) +  (Subject.FootLength .* 2 .* 0.3 .* C3Ddata(z).RightContact(i).UnitVector(2));
        C3Ddata(z).RFA(i,1) = C3Ddata(z).RHEE.Adj(p,1) +  (Subject.FootLength .* 2 .* 0.6 .* C3Ddata(z).RightContact(i).UnitVector(1));  % 2/3rd line - FA = Fore/Arch
        C3Ddata(z).RFA(i,2) = C3Ddata(z).RHEE.Adj(p,2) +  (Subject.FootLength .* 2 .* 0.6 .* C3Ddata(z).RightContact(i).UnitVector(2));
        
        C3Ddata(z).RHAb(i,1) = C3Ddata(z).RHEE.base(i,1) +  (Subject.FootLength .* 2 .* 0.3 .* C3Ddata(z).RightContact(i).UnitVector(1)); % 1/3rd line - HA = Heel/Arch
        C3Ddata(z).RHAb(i,2) = C3Ddata(z).RHEE.base(i,2) +  (Subject.FootLength .* 2 .* 0.3 .* C3Ddata(z).RightContact(i).UnitVector(2));
        C3Ddata(z).RFAb(i,1) = C3Ddata(z).RHEE.base(i,1) +  (Subject.FootLength .* 2 .* 0.6 .* C3Ddata(z).RightContact(i).UnitVector(1));  % 2/3rd line - FA = Fore/Arch
        C3Ddata(z).RFAb(i,2) = C3Ddata(z).RHEE.base(i,2) +  (Subject.FootLength .* 2 .* 0.6 .* C3Ddata(z).RightContact(i).UnitVector(2));
        
        C3Ddata(z).RCent(i,1) = C3Ddata(z).RHEE.Adj(p,1) +  (Subject.FootLength .* 2 .* 0.5 .* C3Ddata(z).RightContact(i).UnitVector(1)); % centroid from adj points
        C3Ddata(z).RCent(i,2) = C3Ddata(z).RHEE.Adj(p,2) +  (Subject.FootLength .* 2 .* 0.5  .* C3Ddata(z).RightContact(i).UnitVector(2));
        C3Ddata(z).RCentB(i,1) = C3Ddata(z).RHEE.base(i,1) +  (Subject.FootLength .* 2 .* 0.5  .* C3Ddata(z).RightContact(i).UnitVector(1)); % centroid from base points
        C3Ddata(z).RCentB(i,2) = C3Ddata(z).RHEE.base(i,2) +  (Subject.FootLength .* 2 .* 0.5  .* C3Ddata(z).RightContact(i).UnitVector(2));
    end
end

%% Quality check coordinate system alignment ML
% clearvars Dist
% % only checks first step on each side
% % FlipThresh = 10; % threshold for flipping images = 2.5 cm (5 half-cm)
% AdjThresh = 1; % threshold for adjusting C3D locations = 0.5 cm (1 half-cm)
% 
% for z = 1:length(C3Ddata)
%     for j = 1:length(Regions{1,z})
%         if strcmp(Regions{1,z}(j).Side, 'Left') % identify feet
%             LFeet(j) = 1;
%             RFeet(j) = 0;
%         else
%             RFeet(j) = 1;
%             LFeet(j) = 0;
%         end
%     end
%     LFeets = find(LFeet);
%     RFeets = find(RFeet);
%     if isempty(LFeets) == 0
%         Dist.LeftX(z) =  C3Ddata(z).LCentB(1,2) - mean(Regions{1,z}(LFeets(1)).HT.General(1,:));
%     else
%         Dist.LeftX(z) = NaN;
%     end
%     if isempty(RFeets) == 0
%         Dist.RightX(z) = C3Ddata(z).RCentB(1,2) - mean(Regions{1,z}(RFeets(1)).HT.General(1,:));
%     else
%         Dist.RightX(z) = NaN;
%     end
%     Dist.Avg(z) = nanmean([abs(Dist.LeftX(z)), abs(Dist.RightX(z))]);
%     clearvars LFeets RFeets LFeet RFeet
% end
% 
% % flip to ensure correct alignment for all 
% % if Dist.Avg(z) > FlipThresh
% %     for i = 1:length(Regions)
% %         DynamicPPTrials(i).TM = fliplr(DynamicPPTrials(i).TM);
% %         DynamicPPTrials(i).SumTM = fliplr(DynamicPPTrials(i).SumTM);
% %         DynamicPPTrials(i).CoP(:,1) = PPOutput.Width - DynamicPPTrials(i).CoP(:,1) +1;
% %         
% %         % Re-do regions
% %         for j = 1:length(Regions{i})
% %             Regions{i}(j).Centroid(1) = PPOutput.Width - Regions{i}(j).Centroid(1);
% %             Regions{i}(j).BoundingBox(1) = PPOutput.Width - Regions{i}(j).BoundingBox(1);
% %         end
% %     end
% % end
% for z = 1:length(C3Ddata)
%     % LEFT adjustment if necessary
% %     if isnan(Dist.LeftX(z)) == 0
%         if Dist.Avg(z)  > AdjThresh
%             if strcmp(C3Ddata(z).Direction, 'North')
%                 % LEFT
%                 C3Ddata(z).LHEE.Adj(:,2) = C3Ddata(z).LHEE.Adj(:,2) + Dist.LeftX(z); % heel points
%                 C3Ddata(z).LHEE.base(:,2) = C3Ddata(z).LHEE.base(:,2) + Dist.LeftX(z);
%                 C3Ddata(z).LTOE.Adj(:,2) = C3Ddata(z).LTOE.Adj(:,2) + Dist.LeftX(z); % toe points
%                 C3Ddata(z).LTOE.top(:,2) = C3Ddata(z).LTOE.top(:,2) + Dist.LeftX(z);
%                 C3Ddata(z).LTOE.step(:,2) = C3Ddata(z).LTOE.step(:,2) + Dist.LeftX(z);
%                 C3Ddata(z).LD1M.Adj(:,2) = C3Ddata(z).LD1M.Adj(:,2) + Dist.LeftX(z); % D1M marker
%                 C3Ddata(z).LD1M.step(:,2) = C3Ddata(z).LD1M.step(:,2) + Dist.LeftX(z);
%                 C3Ddata(z).LD5M.Adj(:,2) = C3Ddata(z).LD5M.Adj(:,2) + Dist.LeftX(z); % D5M marker
%                 C3Ddata(z).LD5M.step(:,2) = C3Ddata(z).LD5M.step(:,2) + Dist.LeftX(z);
%                 C3Ddata(z).LP1M.Adj(:,2) = C3Ddata(z).LP1M.Adj(:,2) + Dist.LeftX(z); % P1M marker
%                 C3Ddata(z).LP1M.step(:,2) = C3Ddata(z).LP1M.step(:,2) + Dist.LeftX(z);
%                 C3Ddata(z).LP5M.Adj(:,2) = C3Ddata(z).LP5M.Adj(:,2) + Dist.LeftX(z); % P5M marker
%                 C3Ddata(z).LP5M.step(:,2) = C3Ddata(z).LP5M.step(:,2) + Dist.LeftX(z);
%                 C3Ddata(z).LLCA.Adj(:,2) = C3Ddata(z).LLCA.Adj(:,2) + Dist.LeftX(z); % LCA marker
%                 C3Ddata(z).LLCA.step(:,2) = C3Ddata(z).LLCA.step(:,2) + Dist.LeftX(z);
%                 C3Ddata(z).LMCA.Adj(:,2) = C3Ddata(z).LMCA.Adj(:,2) + Dist.LeftX(z); % MCA marker
%                 C3Ddata(z).LMCA.step(:,2) = C3Ddata(z).LMCA.step(:,2) + Dist.LeftX(z);
%                 % RIGHT
%                 C3Ddata(z).RHEE.Adj(:,2) = C3Ddata(z).RHEE.Adj(:,2) + Dist.RightX(z); % heel points
%                 C3Ddata(z).RHEE.base(:,2) = C3Ddata(z).RHEE.base(:,2) + Dist.RightX(z);
%                 C3Ddata(z).RTOE.Adj(:,2) = C3Ddata(z).RTOE.Adj(:,2) + Dist.RightX(z); % toe points
%                 C3Ddata(z).RTOE.top(:,2) = C3Ddata(z).RTOE.top(:,2) + Dist.RightX(z);
%                 C3Ddata(z).RTOE.step(:,2) = C3Ddata(z).RTOE.step(:,2) + Dist.RightX(z);
%                 C3Ddata(z).RD1M.Adj(:,2) = C3Ddata(z).RD1M.Adj(:,2) + Dist.RightX(z); % D1M marker
%                 C3Ddata(z).RD1M.step(:,2) = C3Ddata(z).RD1M.step(:,2) + Dist.RightX(z);
%                 C3Ddata(z).RD5M.Adj(:,2) = C3Ddata(z).RD5M.Adj(:,2) + Dist.RightX(z); % D5M marker
%                 C3Ddata(z).RD5M.step(:,2) = C3Ddata(z).RD5M.step(:,2) + Dist.RightX(z);
%                 C3Ddata(z).RP1M.Adj(:,2) = C3Ddata(z).RP1M.Adj(:,2) + Dist.RightX(z); % P1M marker
%                 C3Ddata(z).RP1M.step(:,2) = C3Ddata(z).RP1M.step(:,2) + Dist.RightX(z);
%                 C3Ddata(z).RP5M.Adj(:,2) = C3Ddata(z).RP5M.Adj(:,2) + Dist.RightX(z); % P5M marker
%                 C3Ddata(z).RP5M.step(:,2) = C3Ddata(z).RP5M.step(:,2) + Dist.RightX(z);
%                 C3Ddata(z).RLCA.Adj(:,2) = C3Ddata(z).RLCA.Adj(:,2) + Dist.RightX(z); % LCA marker
%                 C3Ddata(z).RLCA.step(:,2) = C3Ddata(z).RLCA.step(:,2) + Dist.RightX(z);
%                 C3Ddata(z).RMCA.Adj(:,2) = C3Ddata(z).RMCA.Adj(:,2) + Dist.RightX(z); % MCA marker
%                 C3Ddata(z).RMCA.step(:,2) = C3Ddata(z).RMCA.step(:,2) + Dist.RightX(z);
%             else
%                 % LEFT
%                 C3Ddata(z).LHEE.Adj(:,2) = C3Ddata(z).LHEE.Adj(:,2) - Dist.LeftX(z); % heel points
%                 C3Ddata(z).LHEE.base(:,2) = C3Ddata(z).LHEE.base(:,2) - Dist.LeftX(z);
%                 C3Ddata(z).LTOE.Adj(:,2) = C3Ddata(z).LTOE.Adj(:,2) - Dist.LeftX(z); % toe points
%                 C3Ddata(z).LTOE.top(:,2) = C3Ddata(z).LTOE.top(:,2) - Dist.LeftX(z);
%                 C3Ddata(z).LTOE.step(:,2) = C3Ddata(z).LTOE.step(:,2) - Dist.LeftX(z);
%                 C3Ddata(z).LD1M.Adj(:,2) = C3Ddata(z).LD1M.Adj(:,2) - Dist.LeftX(z); % D1M marker
%                 C3Ddata(z).LD1M.step(:,2) = C3Ddata(z).LD1M.step(:,2) - Dist.LeftX(z);
%                 C3Ddata(z).LD5M.Adj(:,2) = C3Ddata(z).LD5M.Adj(:,2) - Dist.LeftX(z); % D5M marker
%                 C3Ddata(z).LD5M.step(:,2) = C3Ddata(z).LD5M.step(:,2) - Dist.LeftX(z);
%                 C3Ddata(z).LP1M.Adj(:,2) = C3Ddata(z).LP1M.Adj(:,2) - Dist.LeftX(z); % P1M marker
%                 C3Ddata(z).LP1M.step(:,2) = C3Ddata(z).LP1M.step(:,2) - Dist.LeftX(z);
%                 C3Ddata(z).LP5M.Adj(:,2) = C3Ddata(z).LP5M.Adj(:,2) - Dist.LeftX(z); % P5M marker
%                 C3Ddata(z).LP5M.step(:,2) = C3Ddata(z).LP5M.step(:,2) - Dist.LeftX(z);
%                 C3Ddata(z).LLCA.Adj(:,2) = C3Ddata(z).LLCA.Adj(:,2) - Dist.LeftX(z); % LCA marker
%                 C3Ddata(z).LLCA.step(:,2) = C3Ddata(z).LLCA.step(:,2) - Dist.LeftX(z);
%                 C3Ddata(z).LMCA.Adj(:,2) = C3Ddata(z).LMCA.Adj(:,2) - Dist.LeftX(z); % MCA marker
%                 C3Ddata(z).LMCA.step(:,2) = C3Ddata(z).LMCA.step(:,2) - Dist.LeftX(z);
%                 % RIGHT
%                 C3Ddata(z).RHEE.Adj(:,2) = C3Ddata(z).RHEE.Adj(:,2) - Dist.RightX(z); % heel points
%                 C3Ddata(z).RHEE.base(:,2) = C3Ddata(z).RHEE.base(:,2) - Dist.RightX(z);
%                 C3Ddata(z).RTOE.Adj(:,2) = C3Ddata(z).RTOE.Adj(:,2) - Dist.RightX(z); % toe points
%                 C3Ddata(z).RTOE.top(:,2) = C3Ddata(z).RTOE.top(:,2) - Dist.RightX(z);
%                 C3Ddata(z).RTOE.step(:,2) = C3Ddata(z).RTOE.step(:,2) - Dist.RightX(z);
%                 C3Ddata(z).RD1M.Adj(:,2) = C3Ddata(z).RD1M.Adj(:,2) - Dist.RightX(z); % D1M marker
%                 C3Ddata(z).RD1M.step(:,2) = C3Ddata(z).RD1M.step(:,2) - Dist.RightX(z);
%                 C3Ddata(z).RD5M.Adj(:,2) = C3Ddata(z).RD5M.Adj(:,2) - Dist.RightX(z); % D5M marker
%                 C3Ddata(z).RD5M.step(:,2) = C3Ddata(z).RD5M.step(:,2) - Dist.RightX(z);
%                 C3Ddata(z).RP1M.Adj(:,2) = C3Ddata(z).RP1M.Adj(:,2) - Dist.RightX(z); % P1M marker
%                 C3Ddata(z).RP1M.step(:,2) = C3Ddata(z).RP1M.step(:,2) - Dist.RightX(z);
%                 C3Ddata(z).RP5M.Adj(:,2) = C3Ddata(z).RP5M.Adj(:,2) - Dist.RightX(z); % P5M marker
%                 C3Ddata(z).RP5M.step(:,2) = C3Ddata(z).RP5M.step(:,2) - Dist.RightX(z);
%                 C3Ddata(z).RLCA.Adj(:,2) = C3Ddata(z).RLCA.Adj(:,2) - Dist.RightX(z); % LCA marker
%                 C3Ddata(z).RLCA.step(:,2) = C3Ddata(z).RLCA.step(:,2) - Dist.RightX(z);
%                 C3Ddata(z).RMCA.Adj(:,2) = C3Ddata(z).RMCA.Adj(:,2) - Dist.RightX(z); % MCA marker
%                 C3Ddata(z).RMCA.step(:,2) = C3Ddata(z).RMCA.step(:,2) - Dist.RightX(z);
%             end
%         end
%     end
%     % RIGHT adjustment if necessary
% %     if isnan(Dist.LeftX(z)) == 0
% %         if abs(mean(Dist.RightX(z))) > AdjThresh
% %             if strcmp(C3Ddata(z).Direction, 'North')
% %           
% %             else
% %                
% %             end
% %         end
% %     end
% % end
% clearvars DelThresh p i j z locs Ti v W

%% Quality check coordinate system alignment AP
% clearvars Dist
% % only checks first step on each side
% AdjThresh = 1; % threshold for adjusting C3D locations = 0.5 cm (1 half-cm)
% for z = 1:length(C3Ddata)
%     for j = 1:length(Regions{1,z})
%         if strcmp(Regions{1,z}(j).Side, 'Left') % identify feet
%             LFeet(j) = 1;
%             RFeet(j) = 0;
%         else
%             RFeet(j) = 1;
%             LFeet(j) = 0;
%         end
%     end
%     LFeets = find(LFeet);
%     RFeets = find(RFeet);
%     if isempty(LFeets) == 0
%         Dist.LeftX(z) =  C3Ddata(z).LCentB(1,1) - mean(Regions{1,z}(LFeets(1)).HT.General(2,:));
%     else
%         Dist.LeftX(z) = NaN;
%     end
%     if isempty(RFeets) == 0
%         Dist.RightX(z) = C3Ddata(z).RCentB(1,1) - mean(Regions{1,z}(RFeets(1)).HT.General(2,:));
%     else
%         Dist.RightX(z) = NaN;
%     end
%     clearvars LFeets RFeets LFeet RFeet
% end
% 
% for z = 1:length(C3Ddata)
%     % LEFT adjustment if necessary
%     if isnan(Dist.LeftX(z)) == 0
%         if abs(mean(Dist.LeftX(z))) > AdjThresh
%             if strcmp(C3Ddata(z).Direction, 'North')
%                 C3Ddata(z).LHEE.Adj(:,1) = C3Ddata(z).LHEE.Adj(:,1) - Dist.LeftX(z); % heel points
%                 C3Ddata(z).LHEE.base(:,1) = C3Ddata(z).LHEE.base(:,1) - Dist.LeftX(z);
%                 C3Ddata(z).LTOE.Adj(:,1) = C3Ddata(z).LTOE.Adj(:,1) - Dist.LeftX(z); % toe points
%                 C3Ddata(z).LTOE.top(:,1) = C3Ddata(z).LTOE.top(:,1) - Dist.LeftX(z);
%                 C3Ddata(z).LTOE.step(:,1) = C3Ddata(z).LTOE.step(:,1) - Dist.LeftX(z);
%                 C3Ddata(z).LD1M.Adj(:,1) = C3Ddata(z).LD1M.Adj(:,1) - Dist.LeftX(z); % D1M marker
%                 C3Ddata(z).LD1M.step(:,1) = C3Ddata(z).LD1M.step(:,1) - Dist.LeftX(z);
%                 C3Ddata(z).LD5M.Adj(:,1) = C3Ddata(z).LD5M.Adj(:,1) - Dist.LeftX(z); % D5M marker
%                 C3Ddata(z).LD5M.step(:,1) = C3Ddata(z).LD5M.step(:,1) - Dist.LeftX(z);
%                 C3Ddata(z).LP1M.Adj(:,1) = C3Ddata(z).LP1M.Adj(:,1) - Dist.LeftX(z); % P1M marker
%                 C3Ddata(z).LP1M.step(:,1) = C3Ddata(z).LP1M.step(:,1) - Dist.LeftX(z);
%                 C3Ddata(z).LP5M.Adj(:,1) = C3Ddata(z).LP5M.Adj(:,1) - Dist.LeftX(z); % P5M marker
%                 C3Ddata(z).LP5M.step(:,1) = C3Ddata(z).LP5M.step(:,1) - Dist.LeftX(z);
%                 C3Ddata(z).LLCA.Adj(:,1) = C3Ddata(z).LLCA.Adj(:,1) - Dist.LeftX(z); % LCA marker
%                 C3Ddata(z).LLCA.step(:,1) = C3Ddata(z).LLCA.step(:,1) - Dist.LeftX(z);
%                 C3Ddata(z).LMCA.Adj(:,1) = C3Ddata(z).LMCA.Adj(:,1) - Dist.LeftX(z); % MCA marker
%                 C3Ddata(z).LMCA.step(:,1) = C3Ddata(z).LMCA.step(:,1) - Dist.LeftX(z);
%             else
%                 C3Ddata(z).LHEE.Adj(:,1) = C3Ddata(z).LHEE.Adj(:,1) - Dist.LeftX(z); % heel points
%                 C3Ddata(z).LHEE.base(:,1) = C3Ddata(z).LHEE.base(:,1) - Dist.LeftX(z);
%                 C3Ddata(z).LTOE.Adj(:,1) = C3Ddata(z).LTOE.Adj(:,1) - Dist.LeftX(z); % toe points
%                 C3Ddata(z).LTOE.top(:,1) = C3Ddata(z).LTOE.top(:,1) - Dist.LeftX(z);
%                 C3Ddata(z).LTOE.step(:,1) = C3Ddata(z).LTOE.step(:,1) - Dist.LeftX(z);
%                 C3Ddata(z).LD1M.Adj(:,1) = C3Ddata(z).LD1M.Adj(:,1) - Dist.LeftX(z); % D1M marker
%                 C3Ddata(z).LD1M.step(:,1) = C3Ddata(z).LD1M.step(:,1) - Dist.LeftX(z);
%                 C3Ddata(z).LD5M.Adj(:,1) = C3Ddata(z).LD5M.Adj(:,1) - Dist.LeftX(z); % D5M marker
%                 C3Ddata(z).LD5M.step(:,1) = C3Ddata(z).LD5M.step(:,1) - Dist.LeftX(z);
%                 C3Ddata(z).LP1M.Adj(:,1) = C3Ddata(z).LP1M.Adj(:,1) - Dist.LeftX(z); % P1M marker
%                 C3Ddata(z).LP1M.step(:,1) = C3Ddata(z).LP1M.step(:,1) - Dist.LeftX(z);
%                 C3Ddata(z).LP5M.Adj(:,1) = C3Ddata(z).LP5M.Adj(:,1) - Dist.LeftX(z); % P5M marker
%                 C3Ddata(z).LP5M.step(:,1) = C3Ddata(z).LP5M.step(:,1) - Dist.LeftX(z);
%                 C3Ddata(z).LLCA.Adj(:,1) = C3Ddata(z).LLCA.Adj(:,1) - Dist.LeftX(z); % LCA marker
%                 C3Ddata(z).LLCA.step(:,1) = C3Ddata(z).LLCA.step(:,1) - Dist.LeftX(z);
%                 C3Ddata(z).LMCA.Adj(:,1) = C3Ddata(z).LMCA.Adj(:,1) - Dist.LeftX(z); % MCA marker
%                 C3Ddata(z).LMCA.step(:,1) = C3Ddata(z).LMCA.step(:,1) - Dist.LeftX(z);
%             end
%         end
%     end
%     % RIGHT adjustment if necessary
%     if isnan(Dist.LeftX(z)) == 0
%         if abs(mean(Dist.RightX(z))) > AdjThresh
%             if strcmp(C3Ddata(z).Direction, 'North')
%                 C3Ddata(z).RHEE.Adj(:,1) = C3Ddata(z).RHEE.Adj(:,1) - Dist.RightX(z); % heel points
%                 C3Ddata(z).RHEE.base(:,1) = C3Ddata(z).RHEE.base(:,1) - Dist.RightX(z);
%                 C3Ddata(z).RTOE.Adj(:,1) = C3Ddata(z).RTOE.Adj(:,1) - Dist.RightX(z); % toe points
%                 C3Ddata(z).RTOE.top(:,1) = C3Ddata(z).RTOE.top(:,1) - Dist.RightX(z);
%                 C3Ddata(z).RTOE.step(:,1) = C3Ddata(z).RTOE.step(:,1) - Dist.RightX(z);
%                 C3Ddata(z).RD1M.Adj(:,1) = C3Ddata(z).RD1M.Adj(:,1) - Dist.RightX(z); % D1M marker
%                 C3Ddata(z).RD1M.step(:,1) = C3Ddata(z).RD1M.step(:,1) - Dist.RightX(z);
%                 C3Ddata(z).RD5M.Adj(:,1) = C3Ddata(z).RD5M.Adj(:,1) - Dist.RightX(z); % D5M marker
%                 C3Ddata(z).RD5M.step(:,1) = C3Ddata(z).RD5M.step(:,1) - Dist.RightX(z);
%                 C3Ddata(z).RP1M.Adj(:,1) = C3Ddata(z).RP1M.Adj(:,1) - Dist.RightX(z); % P1M marker
%                 C3Ddata(z).RP1M.step(:,1) = C3Ddata(z).RP1M.step(:,1) - Dist.RightX(z);
%                 C3Ddata(z).RP5M.Adj(:,1) = C3Ddata(z).RP5M.Adj(:,1) - Dist.RightX(z); % P5M marker
%                 C3Ddata(z).RP5M.step(:,1) = C3Ddata(z).RP5M.step(:,1) - Dist.RightX(z);
%                 C3Ddata(z).RLCA.Adj(:,1) = C3Ddata(z).RLCA.Adj(:,1) - Dist.RightX(z); % LCA marker
%                 C3Ddata(z).RLCA.step(:,1) = C3Ddata(z).RLCA.step(:,1) - Dist.RightX(z);
%                 C3Ddata(z).RMCA.Adj(:,1) = C3Ddata(z).RMCA.Adj(:,1) - Dist.RightX(z); % MCA marker
%                 C3Ddata(z).RMCA.step(:,1) = C3Ddata(z).RMCA.step(:,1) - Dist.RightX(z);
%             else
%                 C3Ddata(z).RHEE.Adj(:,1) = C3Ddata(z).RHEE.Adj(:,1) - Dist.RightX(z); % heel points
%                 C3Ddata(z).RHEE.base(:,1) = C3Ddata(z).RHEE.base(:,1) - Dist.RightX(z);
%                 C3Ddata(z).RTOE.Adj(:,1) = C3Ddata(z).RTOE.Adj(:,1) - Dist.RightX(z); % toe points
%                 C3Ddata(z).RTOE.top(:,1) = C3Ddata(z).RTOE.top(:,1) - Dist.RightX(z);
%                 C3Ddata(z).RTOE.step(:,1) = C3Ddata(z).RTOE.step(:,1) - Dist.RightX(z);
%                 C3Ddata(z).RD1M.Adj(:,1) = C3Ddata(z).RD1M.Adj(:,1) - Dist.RightX(z); % D1M marker
%                 C3Ddata(z).RD1M.step(:,1) = C3Ddata(z).RD1M.step(:,1) - Dist.RightX(z);
%                 C3Ddata(z).RD5M.Adj(:,1) = C3Ddata(z).RD5M.Adj(:,1) - Dist.RightX(z); % D5M marker
%                 C3Ddata(z).RD5M.step(:,1) = C3Ddata(z).RD5M.step(:,1) - Dist.RightX(z);
%                 C3Ddata(z).RP1M.Adj(:,1) = C3Ddata(z).RP1M.Adj(:,1) - Dist.RightX(z); % P1M marker
%                 C3Ddata(z).RP1M.step(:,1) = C3Ddata(z).RP1M.step(:,1) - Dist.RightX(z);
%                 C3Ddata(z).RP5M.Adj(:,1) = C3Ddata(z).RP5M.Adj(:,1) - Dist.RightX(z); % P5M marker
%                 C3Ddata(z).RP5M.step(:,1) = C3Ddata(z).RP5M.step(:,1) - Dist.RightX(z);
%                 C3Ddata(z).RLCA.Adj(:,1) = C3Ddata(z).RLCA.Adj(:,1) - Dist.RightX(z); % LCA marker
%                 C3Ddata(z).RLCA.step(:,1) = C3Ddata(z).RLCA.step(:,1) - Dist.RightX(z);
%                 C3Ddata(z).RMCA.Adj(:,1) = C3Ddata(z).RMCA.Adj(:,1) - Dist.RightX(z); % MCA marker
%                 C3Ddata(z).RMCA.step(:,1) = C3Ddata(z).RMCA.step(:,1) - Dist.RightX(z);
%             end
%         end
%     end
% end
% clearvars DelThresh p i j z locs Ti v W

%% delete any extra steps that weren't identified on PP mat
% again using another method just to be sure
clearvars ToDel
for z = 1:length(C3Ddata)
    % identify # of steps
    for j = 1:length(Regions{1,z})
        FindLeft(j) = strcmp(Regions{1,z}(j).Side, 'Left');
        FindRight(j) = strcmp(Regions{1,z}(j).Side, 'Right');
    end
    Lefts(z) = sum(FindLeft);
    Rights(z) = sum(FindRight);
    clearvars FindLeft FindRight
    
    % LEFT
    i = 1;
    for j = 1:size(C3Ddata(z).LHA,1)
        if Lefts(z) == 0
            ToDel(i) = j;
            i = i + 1;
        else
            if C3Ddata(z).LHEE.base(j,1) < 0 || C3Ddata(z).LTOE.top(j,1) < 0
                ToDel(i) = j;
                i = i + 1;
            elseif C3Ddata(z).LHEE.base(j,1) > PPOutput.Length || C3Ddata(z).LTOE.top(j,1) > PPOutput.Length
                ToDel(i) = j;
                i = i + 1;
            end
        end
    end
    if exist('ToDel','var')
        C3Ddata(z).LeftContact(ToDel) = [];
        C3Ddata(z).LHA(ToDel,:) = [];
        C3Ddata(z).LFA(ToDel,:) = [];
        C3Ddata(z).LHAb(ToDel,:) = [];
        C3Ddata(z).LFAb(ToDel,:) = [];
        C3Ddata(z).LCent(ToDel,:) = [];
        C3Ddata(z).LCentB(ToDel,:) = [];
        C3Ddata(z).LHEE.base(ToDel,:) = [];
        C3Ddata(z).LTOE.top(ToDel,:) = [];
        C3Ddata(z).LD1M.step(ToDel,:) = [];
        C3Ddata(z).LD5M.step(ToDel,:) = [];
        C3Ddata(z).LP1M.step(ToDel,:) = [];
        C3Ddata(z).LP5M.step(ToDel,:) = [];
        C3Ddata(z).LLCA.step(ToDel,:) = [];
        C3Ddata(z).LMCA.step(ToDel,:) = [];
        clearvars i ToDel j
    end
    % RIGHT
    i = 1;
    for j = 1:size(C3Ddata(z).RHA,1)
        if Rights(z) == 0
            ToDel(i) = j;
            i = i + 1;
        else
            if C3Ddata(z).RHA(j,1) < 0 || C3Ddata(z).RFA(j,1) < 0
                ToDel(i) = j;
                i = i + 1;
            elseif C3Ddata(z).RHA(j,1) > PPOutput.Length || C3Ddata(z).RFA(j,1) > PPOutput.Length
                ToDel(i) = j;
                i = i + 1;
            end
        end
    end
    if exist('ToDel','var')
        C3Ddata(z).RightContact(ToDel) = [];
        C3Ddata(z).RHA(ToDel,:) = [];
        C3Ddata(z).RFA(ToDel,:) = [];
        C3Ddata(z).RHAb(ToDel,:) = [];
        C3Ddata(z).RFAb(ToDel,:) = [];
        C3Ddata(z).RCent(ToDel,:) = [];
        C3Ddata(z).RCentB(ToDel,:) = [];
        C3Ddata(z).RHEE.base(ToDel,:) = [];
        C3Ddata(z).RTOE.top(ToDel,:) = [];
        C3Ddata(z).RD1M.step(ToDel,:) = [];
        C3Ddata(z).RD5M.step(ToDel,:) = [];
        C3Ddata(z).RP1M.step(ToDel,:) = [];
        C3Ddata(z).RP5M.step(ToDel,:) = [];
        C3Ddata(z).RLCA.step(ToDel,:) = [];
        C3Ddata(z).RMCA.step(ToDel,:) = [];
        clearvars i ToDel j
    end
end

%% Plot Aligned coordinate system
% for z = 1:NumTrials
%     figure;
% %     subplot(121);
%     hold on; grid on;
%     contour(DynamicPPTrials(z).SumTM,25); axis equal; % plot PP pressures
%     % LEFT
%     if isempty(C3Ddata(z).LHA) == 0
%         for i = 1:length(C3Ddata(z).LeftContact)
%             %         plot(C3Ddata(z).LHEE.step(i,2), C3Ddata(z).LHEE.step(i,1), 'r.','Markersize', 20);
%             %         plot(C3Ddata(z).LTOE.step(i,2), C3Ddata(z).LTOE.step(i,1), 'r.','Markersize', 20);
%             plot(C3Ddata(z).LD1M.step(i,2), C3Ddata(z).LD1M.step(i,1), 'r*');
%             plot(C3Ddata(z).LD5M.step(i,2), C3Ddata(z).LD5M.step(i,1), 'r*');
%             plot(C3Ddata(z).LP1M.step(i,2), C3Ddata(z).LP1M.step(i,1), 'r*');
%             plot(C3Ddata(z).LP5M.step(i,2), C3Ddata(z).LP5M.step(i,1), 'r*');
%             plot(C3Ddata(z).LLCA.step(i,2), C3Ddata(z).LLCA.step(i,1), 'r*');
%             plot(C3Ddata(z).LMCA.step(i,2), C3Ddata(z).LMCA.step(i,1), 'r*');
%             plot(C3Ddata(z).LHEE.base(i,2), C3Ddata(z).LHEE.base(i,1), 'm.', 'Markersize', 20);
%             plot(C3Ddata(z).LHAb(i,2), C3Ddata(z).LHA(i,1), 'k*');
%             plot(C3Ddata(z).LFAb(i,2), C3Ddata(z).LFA(i,1), 'k*');
%             plot(C3Ddata(z).LTOE.top(i,2), C3Ddata(z).LTOE.top(i,1), 'ko');
%             %         plot(C3Ddata(z).LHA(i,2), C3Ddata(z).LHA(i,1), 'k*');
%             %         plot(C3Ddata(z).LFA(i,2), C3Ddata(z).LFA(i,1), 'k*');
%         end
%     end
%     % RIGHT
%     if isempty(C3Ddata(z).RHA) == 0
%         for i = 1:length(C3Ddata(z).RightContact)
%             %         plot(C3Ddata(z).RHEE.step(i,2), C3Ddata(z).RHEE.step(i,1), 'g.','Markersize', 16);
%             %         plot(C3Ddata(z).RTOE.step(i,2), C3Ddata(z).RTOE.step(i,1), 'g.','Markersize', 20);
%             plot(C3Ddata(z).RD1M.step(i,2), C3Ddata(z).RD1M.step(i,1), 'g*');
%             plot(C3Ddata(z).RD5M.step(i,2), C3Ddata(z).RD5M.step(i,1), 'g*');
%             plot(C3Ddata(z).RP1M.step(i,2), C3Ddata(z).RP1M.step(i,1), 'g*');
%             plot(C3Ddata(z).RP5M.step(i,2), C3Ddata(z).RP5M.step(i,1), 'g*');
%             plot(C3Ddata(z).RLCA.step(i,2), C3Ddata(z).RLCA.step(i,1), 'g*');
%             plot(C3Ddata(z).RMCA.step(i,2), C3Ddata(z).RMCA.step(i,1), 'g*');
%             plot(C3Ddata(z).RHEE.base(i,2), C3Ddata(z).RHEE.base(i,1), 'm.','MarkerSize', 20);
%             plot(C3Ddata(z).RHAb(i,2), C3Ddata(z).RHA(i,1), 'k*');
%             plot(C3Ddata(z).RFAb(i,2), C3Ddata(z).RFA(i,1), 'k*');
%             plot(C3Ddata(z).RTOE.top(i,2), C3Ddata(z).RTOE.top(i,1), 'ko');
%             %         plot(C3Ddata(z).RHA(i,2), C3Ddata(z).RHA(i,1), 'k*');
%             %         plot(C3Ddata(z).RFA(i,2), C3Ddata(z).RFA(i,1), 'k*');
%         end
%     end
%     rectangle('Position', [0, 0, 440 ./ 5, (3548-2103) ./ 5]);
%     axis equal; ylim([0 300]);
%     title('Top view');
% end

%% Trim a copy of c3d trials with gait events only on the pp mat
TrimC3Ddata = C3Ddata;
for z = 1:length(TrimC3Ddata)
    Trim(z).All = [Trim(z).Left; Trim(z).Right];
    Trim(z).All = Trim(z).All + TrimC3Ddata(z).HeaderInfo(4).data - 1;
    clearvars ToDel A a Offs OffMats
    for i = 1:size(TrimC3Ddata(z).Events.GaitEventsSimple,2)
        a =  double(cell2mat(TrimC3Ddata(z).Events.GaitEventsSimple(3,i))); %+ C3Ddata(z).HeaderInfo(4).data - 1;
        A = find(a == Trim(z).All,1);
        if isempty(A) == 1
            OffMats(i) = 1;
        else
            OffMats(i) = 0;
        end
        if strcmp(TrimC3Ddata(z).Events.GaitEventsSimple(1,i),'Foot Off')
            Offs(i) = 1;
        else
            Offs(i) = 0;
        end
    end
    ToDel = OffMats+Offs;
    ToDel = logical(ToDel>0);
    TrimC3Ddata(z).Events.GaitEventsSimple(:,ToDel) = [];
    TrimC3Ddata(z).Events.GaitEvents = [];
    TrimC3Ddata(z).TempSpat = [];
    % organize L and R trajectory data for temp spat calcs
    clearvars ToDel
    for i = 1:size(TrimC3Ddata(z).Events.GaitEventsSimple,2)
        if strcmp(TrimC3Ddata(z).Events.GaitEventsSimple(2,i), 'Left')
            Lefts(i) = 1;
        else
            Lefts(i) = 0;
        end
        Ind = double(cell2mat(TrimC3Ddata(z).Events.GaitEventsSimple(3,i))) - TrimC3Ddata(z).HeaderInfo(4).data +1;
        if Lefts(i) == 1
            TrimC3Ddata(z).Events.GaitEvents(2,i) = TrimC3Ddata(z).LHEE.raw(Ind,1);
            TrimC3Ddata(z).Events.GaitEvents(3,i) = TrimC3Ddata(z).LHEE.raw(Ind,2);
        else
            TrimC3Ddata(z).Events.GaitEvents(2,i) = TrimC3Ddata(z).RHEE.raw(Ind,1);
            TrimC3Ddata(z).Events.GaitEvents(3,i) = TrimC3Ddata(z).RHEE.raw(Ind,2);
        end
    end
    TrimC3Ddata(z).Events.GaitEvents(1,:) =  double(cell2mat(TrimC3Ddata(z).Events.GaitEventsSimple(3,:)));
    clearvars ToDel
    for i = 1:size(TrimC3Ddata(z).Events.GaitEvents,2) % delete out of range events
        if TrimC3Ddata(z).Events.GaitEvents(2,i) < 2090 ||  TrimC3Ddata(z).Events.GaitEvents(2,i) > 3560
            ToDel(i) = 1;
        else
            ToDel(i) = 0;
        end
    end
    ToDel = logical(ToDel);
    TrimC3Ddata(z).Events.GaitEvents(:,ToDel) = [];
    TrimC3Ddata(z).Events.GaitEvents(4,:) = TrimC3Ddata(z).Events.GaitEvents(1,:) / TrimC3Ddata(z).HeaderInfo(6).data;
    % Temporal spatial calculations
    % Calculate gait speed and cadence
    p = size(TrimC3Ddata(z).Events.GaitEvents,2);
    for j = 1:p-1
        MpS(j) = abs((TrimC3Ddata(z).Events.GaitEvents(2,j+1) - TrimC3Ddata(z).Events.GaitEvents(2,j))) /  (1000*(TrimC3Ddata(z).Events.GaitEvents(4,j+1) - TrimC3Ddata(z).Events.GaitEvents(4,j)));
        Cadence(j) = 60*(p-1) / (TrimC3Ddata(z).Events.GaitEvents(4,p) - TrimC3Ddata(z).Events.GaitEvents(4,1));
        StepW(j) = abs(TrimC3Ddata(z).Events.GaitEvents(3,j+1) - TrimC3Ddata(z).Events.GaitEvents(3,j)) / 1000;
        StrideL(j) = abs(TrimC3Ddata(z).Events.GaitEvents(2,j+1) - TrimC3Ddata(z).Events.GaitEvents(2,j)) / 1000;
    end
    
    TrimSpat(z).GaitSpeed_MpS_Avg = nanmean(MpS); % averages
    TrimSpat(z).GaitSpeed_MpM_Avg = TrimSpat(z).GaitSpeed_MpS_Avg * 60;
    TrimSpat(z).Cadence_Avg = nanmean(Cadence);
    TrimSpat(z).StepWidth_Avg = nanmean(StepW);
    TrimSpat(z).StrideLength_Avg = nanmean(StrideL);
    TrimSpat(z).GaitSpeed_MpS_Measures = MpS; % measures
    TrimSpat(z).GaitSpeed_MpM_Measures = TrimSpat(z).GaitSpeed_MpS_Measures * 60;
    TrimSpat(z).Cadence_Measures = Cadence;
    TrimSpat(z).StepWidth_Measures = StepW;
    TrimSpat(z).StrideLength_Measures = StrideL;
    TrimSpat(z).GaitSpeed_MpS_Std = nanstd(MpS); % standard deviations
    TrimSpat(z).GaitSpeed_MpM_Std = nanstd(TrimSpat(z).GaitSpeed_MpM_Measures);
    TrimSpat(z).Cadence_Std = nanstd(Cadence);
    TrimSpat(z).StepWidth_Std = nanstd(StepW);
    TrimSpat(z).StrideLength_Std = nanstd(StrideL);
    
    clearvars MpS Cadence StepW StrideL v W L i j p a A ans Lefts locs NextEnd OffMats Offs question Question ToDel Xoffset Xoffset
end

TrimSpat = struct2table(TrimSpat);

end