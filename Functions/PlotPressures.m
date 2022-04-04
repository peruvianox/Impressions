function [] = PlotPressures(DynamicPPTrials, C3Ddata, PPSettings, Chosen4)

%% set up figure dimensions
set(0,'units','pixels');
Pix_SS = get(0,'screensize');
Pix_SS(1:2) = [];
DisplayDim = [40, 40, Pix_SS(1)*.66, Pix_SS(2)*.8];
DisplayTrials = figure;
set(DisplayTrials, 'Position', DisplayDim);

%% set number of plots
nPlots = length(Chosen4);
if length(Chosen4) > 4
    nPlots = 4;
    Chosen4(5:end) = [];
end

%% Plot contours
Green = rgb('ForestGreen');
[m,n,~] = size(DynamicPPTrials(1).SumTM);
sp = 0;
for q = Chosen4
    sp = sp+1;
    subplot(1, nPlots, sp); hold on;
    contour(DynamicPPTrials(q).SumTM, 100);
    hline(m, 'k');
    vline(n, 'k');
    if isfield(DynamicPPTrials(q).Zones,'L')
        for j = 1:length(DynamicPPTrials(q).Zones.L) % Left mask regions
            %             pgon = polyshape(DynamicPPTrials(q).Zones.L(j).LatHeel(:,1)', DynamicPPTrials(q).Zones.L(j).LatHeel(:,2)');
            %             plot(pgon,'EdgeColor', 'r', 'FaceColor', 'none');
            %             pgon = polyshape(DynamicPPTrials(q).Zones.L(j).MedHeel(:,1)', DynamicPPTrials(q).Zones.L(j).MedHeel(:,2)');
            %             plot(pgon,'EdgeColor', 'r', 'FaceColor', 'none');
            %             pgon = polyshape(DynamicPPTrials(q).Zones.L(j).LatArch(:,1)', DynamicPPTrials(q).Zones.L(j).LatArch(:,2)');
            %             plot(pgon,'EdgeColor', 'r', 'FaceColor', 'none');
            %             pgon = polyshape(DynamicPPTrials(q).Zones.L(j).MedArch(:,1)', DynamicPPTrials(q).Zones.L(j).MedArch(:,2)');
            %             plot(pgon,'EdgeColor', 'r', 'FaceColor', 'none');
            %             pgon = polyshape(DynamicPPTrials(q).Zones.L(j).LatFore(:,1)', DynamicPPTrials(q).Zones.L(j).LatFore(:,2)');
            %             plot(pgon,'EdgeColor', 'r', 'FaceColor', 'none');
            %             pgon = polyshape(DynamicPPTrials(q).Zones.L(j).MedFore(:,1)', DynamicPPTrials(q).Zones.L(j).MedFore(:,2)');
            %             plot(pgon,'EdgeColor', 'r', 'FaceColor', 'none');
            
            pgon = [DynamicPPTrials(q).Zones.L(j).LatHeel(:,1), DynamicPPTrials(q).Zones.L(j).LatHeel(:,2)];
            pgon(end+1, :) = pgon(1, :);
            line(pgon(:,1), pgon(:,2),'Color', 'r');
            clearvars pgon
            pgon = [DynamicPPTrials(q).Zones.L(j).MedHeel(:,1), DynamicPPTrials(q).Zones.L(j).MedHeel(:,2)];
            pgon(end+1, :) = pgon(1, :);
            line(pgon(:,1), pgon(:,2),'Color', 'r');
            clearvars pgon
            pgon = [DynamicPPTrials(q).Zones.L(j).LatArch(:,1), DynamicPPTrials(q).Zones.L(j).LatArch(:,2)];
            pgon(end+1, :) = pgon(1, :);
            line(pgon(:,1), pgon(:,2),'Color', 'r');
            clearvars pgon
            pgon = [DynamicPPTrials(q).Zones.L(j).MedArch(:,1), DynamicPPTrials(q).Zones.L(j).MedArch(:,2)];
            pgon(end+1, :) = pgon(1, :);
            line(pgon(:,1), pgon(:,2),'Color', 'r');
            clearvars pgon
            pgon = [DynamicPPTrials(q).Zones.L(j).LatFore(:,1), DynamicPPTrials(q).Zones.L(j).LatFore(:,2)];
            pgon(end+1, :) = pgon(1, :);
            line(pgon(:,1), pgon(:,2),'Color', 'r');
            clearvars pgon
            pgon = [DynamicPPTrials(q).Zones.L(j).MedFore(:,1), DynamicPPTrials(q).Zones.L(j).MedFore(:,2)];
            pgon(end+1, :) = pgon(1, :);
            line(pgon(:,1), pgon(:,2),'Color', 'r');
            clearvars pgon
            
        end
    end
    if isfield(DynamicPPTrials(q).Zones,'R')
        for j = 1:length(DynamicPPTrials(q).Zones.R) % right mask regions
            %             pgon = polyshape(DynamicPPTrials(q).Zones.R(j).LatHeel(:,1)', DynamicPPTrials(q).Zones.R(j).LatHeel(:,2)');
            %             plot(pgon,'EdgeColor', Green, 'FaceColor', 'none');
            %             pgon = polyshape(DynamicPPTrials(q).Zones.R(j).MedHeel(:,1)', DynamicPPTrials(q).Zones.R(j).MedHeel(:,2)');
            %             plot(pgon,'EdgeColor', Green, 'FaceColor', 'none');
            %             pgon = polyshape(DynamicPPTrials(q).Zones.R(j).LatArch(:,1)', DynamicPPTrials(q).Zones.R(j).LatArch(:,2)');
            %             plot(pgon,'EdgeColor', Green, 'FaceColor', 'none');
            %             pgon = polyshape(DynamicPPTrials(q).Zones.R(j).MedArch(:,1)', DynamicPPTrials(q).Zones.R(j).MedArch(:,2)');
            %             plot(pgon,'EdgeColor', Green, 'FaceColor', 'none');
            %             pgon = polyshape(DynamicPPTrials(q).Zones.R(j).LatFore(:,1)', DynamicPPTrials(q).Zones.R(j).LatFore(:,2)');
            %             plot(pgon,'EdgeColor', Green, 'FaceColor', 'none');
            %             pgon = polyshape(DynamicPPTrials(q).Zones.R(j).MedFore(:,1)', DynamicPPTrials(q).Zones.R(j).MedFore(:,2)');
            %             plot(pgon,'EdgeColor', Green, 'FaceColor', 'none');
            
            pgon = [DynamicPPTrials(q).Zones.R(j).LatHeel(:,1), DynamicPPTrials(q).Zones.R(j).LatHeel(:,2)];
            pgon(end+1, :) = pgon(1, :);
            line(pgon(:,1), pgon(:,2),'Color', 'Green');
            clearvars pgon
            pgon = [DynamicPPTrials(q).Zones.R(j).MedHeel(:,1), DynamicPPTrials(q).Zones.R(j).MedHeel(:,2)];
            pgon(end+1, :) = pgon(1, :);
            line(pgon(:,1), pgon(:,2),'Color', 'Green');
            clearvars pgon
            pgon = [DynamicPPTrials(q).Zones.R(j).LatArch(:,1), DynamicPPTrials(q).Zones.R(j).LatArch(:,2)];
            pgon(end+1, :) = pgon(1, :);
            line(pgon(:,1), pgon(:,2),'Color', 'Green');
            clearvars pgon
            pgon = [DynamicPPTrials(q).Zones.R(j).MedArch(:,1), DynamicPPTrials(q).Zones.R(j).MedArch(:,2)];
            pgon(end+1, :) = pgon(1, :);
            line(pgon(:,1), pgon(:,2),'Color', 'Green');
            clearvars pgon
            pgon = [DynamicPPTrials(q).Zones.R(j).LatFore(:,1), DynamicPPTrials(q).Zones.R(j).LatFore(:,2)];
            pgon(end+1, :) = pgon(1, :);
            line(pgon(:,1), pgon(:,2),'Color', 'Green');
            clearvars pgon
            pgon = [DynamicPPTrials(q).Zones.R(j).MedFore(:,1), DynamicPPTrials(q).Zones.R(j).MedFore(:,2)];
            pgon(end+1, :) = pgon(1, :);
            line(pgon(:,1), pgon(:,2),'Color', 'Green');
            clearvars pgon
        end
    end
end

%% add C3D data if present
if strcmp(PPSettings.C3DInput, 'Yes')
    MkrSz = 5;
    sp = 0;
    for q = Chosen4
        sp = sp+1;
        subplot(1, nPlots, sp); hold on;
        plot(C3Ddata(q).LHEE.base(:,2), C3Ddata(q).LHEE.base(:,1),'k*', 'MarkerSize', MkrSz);
        plot(C3Ddata(q).LTOE.top(:,2), C3Ddata(q).LTOE.top(:,1),'ko', 'MarkerSize', MkrSz);
        plot(C3Ddata(q).LTOE.step(:,2), C3Ddata(q).LTOE.step(:,1),'k*', 'MarkerSize', MkrSz);
        plot(C3Ddata(q).LD1M.step(:,2), C3Ddata(q).LD1M.step(:,1),'k*', 'MarkerSize', MkrSz);
        plot(C3Ddata(q).LD5M.step(:,2), C3Ddata(q).LD5M.step(:,1),'k*', 'MarkerSize', MkrSz);
        plot(C3Ddata(q).LP1M.step(:,2), C3Ddata(q).LP1M.step(:,1),'k*', 'MarkerSize', MkrSz);
        plot(C3Ddata(q).LP5M.step(:,2), C3Ddata(q).LP5M.step(:,1),'k*', 'MarkerSize', MkrSz);
        plot(C3Ddata(q).LLCA.step(:,2), C3Ddata(q).LLCA.step(:,1),'k*', 'MarkerSize', MkrSz);
        plot(C3Ddata(q).LMCA.step(:,2), C3Ddata(q).LMCA.step(:,1),'k*', 'MarkerSize', MkrSz);
        plot(C3Ddata(q).LHAb(:,2), C3Ddata(q).LHAb(:,1),'ko', 'MarkerSize', MkrSz);
        plot(C3Ddata(q).LFAb(:,2), C3Ddata(q).LFAb(:,1),'ko', 'MarkerSize', MkrSz);
        
        plot(C3Ddata(q).RHEE.base(:,2), C3Ddata(q).RHEE.base(:,1),'k*', 'MarkerSize', MkrSz);
        plot(C3Ddata(q).RTOE.top(:,2), C3Ddata(q).RTOE.top(:,1),'ko', 'MarkerSize', MkrSz);
        plot(C3Ddata(q).RTOE.step(:,2), C3Ddata(q).RTOE.step(:,1),'k*', 'MarkerSize', MkrSz);
        plot(C3Ddata(q).RD1M.step(:,2), C3Ddata(q).RD1M.step(:,1),'k*', 'MarkerSize', MkrSz);
        plot(C3Ddata(q).RD5M.step(:,2), C3Ddata(q).RD5M.step(:,1),'k*', 'MarkerSize', MkrSz);
        plot(C3Ddata(q).RP1M.step(:,2), C3Ddata(q).RP1M.step(:,1),'k*', 'MarkerSize', MkrSz);
        plot(C3Ddata(q).RP5M.step(:,2), C3Ddata(q).RP5M.step(:,1),'k*', 'MarkerSize', MkrSz);
        plot(C3Ddata(q).RLCA.step(:,2), C3Ddata(q).RLCA.step(:,1),'k*', 'MarkerSize', MkrSz);
        plot(C3Ddata(q).RMCA.step(:,2), C3Ddata(q).RMCA.step(:,1),'k*', 'MarkerSize', MkrSz);
        plot(C3Ddata(q).RHAb(:,2), C3Ddata(q).RHAb(:,1),'ko', 'MarkerSize', MkrSz);
        plot(C3Ddata(q).RFAb(:,2), C3Ddata(q).RFAb(:,1),'ko', 'MarkerSize', MkrSz);
    end
end

%% apply settings for saving
sp = 0;
for q = Chosen4
    sp = sp+1;
    subplot(1, nPlots, sp);
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
    axis equal;
end

end
