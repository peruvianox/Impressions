function [DynamicPPTrials] = SearchParallel(DynamicPPTrials, Selection, Block, PPSettings)

%% initializing
% define search "area"
Srch = 2;
% set loop distances
Loop = 1:0.1:20;

%% loop through trials to adjust points
for q = 1:length(Selection)
    % define matrix size
    [m,n] = size(DynamicPPTrials(Selection(q)).MaskLog);
    %% LEFT
    if DynamicPPTrials(Selection(q)).NumLeft > 0
        for k = 1: DynamicPPTrials(Selection(q)).NumLeft
            % Manual points
            if isfield(DynamicPPTrials(q).Lht, 'Manual') == 1
                if strcmp(PPSettings.MaskChoice, 'Manual') || strcmp(PPSettings.MaskChoice, 'Manual (pink)') || strcmp(PPSettings.MaskChoice, 'Validation')
                    for i = Loop
                        xAdj = DynamicPPTrials(Selection(q)).Lht.ManualInV{k}(:,1);
                        yAdj = DynamicPPTrials(Selection(q)).Lht.ManualInV{k}(:,2);
                        a = DynamicPPTrials(Selection(q)).Lht.Manual{k}(:,1)+(i*xAdj);
                        b = DynamicPPTrials(Selection(q)).Lht.Manual{k}(:,2)+(i*yAdj);
                        SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                        SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                        if sum(sum(isnan([SrchAreaX, SrchAreaY]))) == 0
                            Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                            Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                            if sum(sum(Region)) == 0
                                %                         subplot(1,NumDynTrials, Selection(q));
                                %                         line(a,b,'Color','k');
                                % uncomment to check search points
                                %                                 subplot(1, NumDynTrials, q);
                                %                                 hold on;
                                %                                 for ii = 1:4
                                %                                     plot(SrchAreaX(ii), SrchAreaY(ii),'k.')
                                %                                 end
                                %                                 plot(a(:), b(:), 'r')
                                break
                            end
                        end
                    end
                    if a(1) < DynamicPPTrials(Selection(q)).Lht.Manual{k}(1,1)
                        DynamicPPTrials(Selection(q)).L.LatManual{k} = [a b];
                    else
                        DynamicPPTrials(Selection(q)).L.MedManual{k} = [a b];
                    end
                    % search in the other direction
                    for i = Loop
                        xAdj = DynamicPPTrials(Selection(q)).Lht.ManualInV{k}(:,1);
                        yAdj = DynamicPPTrials(Selection(q)).Lht.ManualInV{k}(:,2);
                        a = DynamicPPTrials(Selection(q)).Lht.Manual{k}(:,1)-(i*xAdj);
                        b = DynamicPPTrials(Selection(q)).Lht.Manual{k}(:,2)-(i*yAdj);
                        SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                        SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                        if sum(sum(isnan([SrchAreaX, SrchAreaY]))) == 0
                            Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                            Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                            if sum(sum(Region)) == 0
                                %                         subplot(1,NumDynTrials, Selection(q));
                                %                         line(a,b,'Color','k');
                                % uncomment to check search points
                                %                 subplot(1, NumDynTrials, q);
                                %                 hold on;
                                %                 for ii = 1:4
                                %                     plot(SrchAreaX(ii), SrchAreaY(ii),'k.')
                                %                 end
                                %                 plot(a(:), b(:), 'r')
                                break
                            end
                        end
                    end
                    if a(1) < DynamicPPTrials(Selection(q)).Lht.Manual{k}(1,1)
                        DynamicPPTrials(Selection(q)).L.LatManual{k} = [a b];
                    else
                        DynamicPPTrials(Selection(q)).L.MedManual{k} = [a b];
                    end
                    clearvars xAdj yAdj a b SrchAreaX SrchAreaY Zone Region
                end
            end
            
            % General image recognition points
            for i = Loop
                xAdj = DynamicPPTrials(Selection(q)).Lht.GeneralInV{k}(:,1);
                yAdj = DynamicPPTrials(Selection(q)).Lht.GeneralInV{k}(:,2);
                a = DynamicPPTrials(Selection(q)).Lht.General{k}(:,1)+(i*xAdj);
                b = DynamicPPTrials(Selection(q)).Lht.General{k}(:,2)+(i*yAdj);
                SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                if sum(sum(isnan([SrchAreaX, SrchAreaY]))) == 0
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) == 0
                        %                         subplot(1,NumDynTrials, Selection(q));
                        %                         line(a,b,'Color','k');
                        % uncomment to check search points
                        %                                 subplot(1, NumDynTrials, q);
                        %                                 hold on;
                        %                                 for ii = 1:4
                        %                                     plot(SrchAreaX(ii), SrchAreaY(ii),'k.')
                        %                                 end
                        %                                 plot(a(:), b(:), 'r')
                        break
                    end
                end
            end
            if a(1) < DynamicPPTrials(Selection(q)).Lht.General{k}(1,1)
                DynamicPPTrials(Selection(q)).L.LatGeneral{k} = [a b];
            else
                DynamicPPTrials(Selection(q)).L.MedGeneral{k} = [a b];
            end
            % search in the other direction
            for i = Loop
                xAdj = DynamicPPTrials(Selection(q)).Lht.GeneralInV{k}(:,1);
                yAdj = DynamicPPTrials(Selection(q)).Lht.GeneralInV{k}(:,2);
                a = DynamicPPTrials(Selection(q)).Lht.General{k}(:,1)-(i*xAdj);
                b = DynamicPPTrials(Selection(q)).Lht.General{k}(:,2)-(i*yAdj);
                SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                if sum(sum(isnan([SrchAreaX, SrchAreaY]))) == 0
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) == 0
                        %                         subplot(1,NumDynTrials, Selection(q));
                        %                         line(a,b,'Color','k');
                        % uncomment to check search points
                        %                 subplot(1, NumDynTrials, q);
                        %                 hold on;
                        %                 for ii = 1:4
                        %                     plot(SrchAreaX(ii), SrchAreaY(ii),'k.')
                        %                 end
                        %                 plot(a(:), b(:), 'r')
                        break
                    end
                end
            end
            if a(1) < DynamicPPTrials(Selection(q)).Lht.General{k}(1,1)
                DynamicPPTrials(Selection(q)).L.LatGeneral{k} = [a b];
            else
                DynamicPPTrials(Selection(q)).L.MedGeneral{k} = [a b];
            end
            clearvars xAdj yAdj a b SrchAreaX SrchAreaY Zone Region
            
            % CoP points
            for i = Loop
                xAdj = DynamicPPTrials(Selection(q)).Lht.InVCoP{k}(:,1);
                yAdj = DynamicPPTrials(Selection(q)).Lht.InVCoP{k}(:,2);
                a = DynamicPPTrials(Selection(q)).Lht.CoP{k}(:,1)+(i*xAdj);
                b = DynamicPPTrials(Selection(q)).Lht.CoP{k}(:,2)+(i*yAdj);
                SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                if sum(sum(isnan([SrchAreaX, SrchAreaY]))) == 0
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) == 0
                        %                 subplot(1,NumDynTrials, Selection(q));
                        %                 line(a,b,'Color','k');
                        break
                    end
                end
                if a(1) < DynamicPPTrials(Selection(q)).Lht.CoP{k}(1,1)
                    DynamicPPTrials(Selection(q)).L.LatCoP{k} = [a b];
                else
                    DynamicPPTrials(Selection(q)).L.MedCoP{k} = [a b];
                end
            end
            % search in the other direction
            for i = Loop
                xAdj = DynamicPPTrials(Selection(q)).Lht.InVCoP{k}(:,1);
                yAdj = DynamicPPTrials(Selection(q)).Lht.InVCoP{k}(:,2);
                a = DynamicPPTrials(Selection(q)).Lht.CoP{k}(:,1)-(i*xAdj);
                b = DynamicPPTrials(Selection(q)).Lht.CoP{k}(:,2)-(i*yAdj);
                SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                if sum(sum(isnan([SrchAreaX, SrchAreaY]))) == 0
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) == 0
                        %                 subplot(1,NumDynTrials, Selection(q));
                        %                 line(a,b,'Color','k');
                        break
                    end
                end
            end
            if a(1) < DynamicPPTrials(Selection(q)).Lht.CoP{k}(1,1)
                DynamicPPTrials(Selection(q)).L.LatCoP{k} = [a b];
            else
                DynamicPPTrials(Selection(q)).L.MedCoP{k} = [a b];
            end
            clearvars xAdj yAdj a b SrchAreaX SrchAreaY Zone Region
            
            % 66% CoP points
            for i = Loop
                xAdj = DynamicPPTrials(Selection(q)).Lht.CoPInV66{k}(:,1);
                yAdj = DynamicPPTrials(Selection(q)).Lht.CoPInV66{k}(:,2);
                a = DynamicPPTrials(Selection(q)).Lht.CoP66{k}(:,1)+(i*xAdj);
                b = DynamicPPTrials(Selection(q)).Lht.CoP66{k}(:,2)+(i*yAdj);
                SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                if sum(sum(isnan([SrchAreaX, SrchAreaY]))) == 0
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) == 0
                        %                 subplot(1,NumDynTrials, Selection(q));
                        %                 line(a,b,'Color','k');
                        break
                    end
                end
            end
            if a(1) < DynamicPPTrials(Selection(q)).Lht.CoP66{k}(1,1)
                DynamicPPTrials(Selection(q)).L.Lat66{k} = [a b];
            else
                DynamicPPTrials(Selection(q)).L.Med66{k} = [a b];
            end
            % search in the other direction
            for i = Loop
                xAdj = DynamicPPTrials(Selection(q)).Lht.CoPInV66{k}(:,1);
                yAdj = DynamicPPTrials(Selection(q)).Lht.CoPInV66{k}(:,2);
                a = DynamicPPTrials(Selection(q)).Lht.CoP66{k}(:,1)-(i*xAdj);
                b = DynamicPPTrials(Selection(q)).Lht.CoP66{k}(:,2)-(i*yAdj);
                SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                if sum(sum(isnan([SrchAreaX, SrchAreaY]))) == 0
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) == 0
                        %                 subplot(1,NumDynTrials, Selection(q));
                        %                 line(a,b,'Color','k');
                        break
                    end
                end
            end
            if a(1) < DynamicPPTrials(Selection(q)).Lht.CoP66{k}(1,1)
                DynamicPPTrials(Selection(q)).L.Lat66{k} = [a b];
            else
                DynamicPPTrials(Selection(q)).L.Med66{k} = [a b];
            end
            clearvars xAdj yAdj a b SrchAreaX SrchAreaY Zone Region
            
            % inter peak points
            if Block.IP == 0
                for i = Loop
                    xAdj = DynamicPPTrials(Selection(q)).Lht.InVIP{k}(:,1);
                    yAdj = DynamicPPTrials(Selection(q)).Lht.InVIP{k}(:,2);
                    a = DynamicPPTrials(Selection(q)).Lht.IP{k}(:,1)+(i*xAdj);
                    b = DynamicPPTrials(Selection(q)).Lht.IP{k}(:,2)+(i*yAdj);
                    SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                    SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                    if sum(sum(isnan([SrchAreaX, SrchAreaY]))) == 0
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) == 0
                            %                 subplot(1,NumDynTrials, Selection(q));
                            %                 line(a,b,'Color','k');
                            break
                        end
                    end
                end
                if a(1) < DynamicPPTrials(Selection(q)).Lht.IP{k}(1,1)
                    DynamicPPTrials(Selection(q)).L.LatIP{k} = [a b];
                else
                    DynamicPPTrials(Selection(q)).L.MedIP{k} = [a b];
                end
                clearvars xAdj yAdj a b SrchAreaX SrchAreaY Zone Region
                % search in the other direction
                for i = Loop
                    xAdj = DynamicPPTrials(Selection(q)).Lht.InVIP{k}(:,1);
                    yAdj = DynamicPPTrials(Selection(q)).Lht.InVIP{k}(:,2);
                    a = DynamicPPTrials(Selection(q)).Lht.IP{k}(:,1)-(i*xAdj);
                    b = DynamicPPTrials(Selection(q)).Lht.IP{k}(:,2)-(i*yAdj);
                    SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                    SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                    if sum(sum(isnan([SrchAreaX, SrchAreaY]))) == 0
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) == 0
                            %                 subplot(1,NumDynTrials, Selection(q));
                            %                 line(a,b,'Color','k');
                            break
                        end
                    end
                end
                if a(1) < DynamicPPTrials(Selection(q)).Lht.IP{k}(1,1)
                    DynamicPPTrials(Selection(q)).L.LatIP{k} = [a b];
                else
                    DynamicPPTrials(Selection(q)).L.MedIP{k} = [a b];
                end
                clearvars xAdj yAdj a b SrchAreaX SrchAreaY Zone Region
            end
            
            % heel to centroid points
            if Block.HC == 0
                for i = Loop
                    xAdj = DynamicPPTrials(Selection(q)).Lht.InVHeelCent{k}(:,1);
                    yAdj = DynamicPPTrials(Selection(q)).Lht.InVHeelCent{k}(:,2);
                    a = DynamicPPTrials(Selection(q)).Lht.HeelCent{k}(:,1)+(i*xAdj);
                    b = DynamicPPTrials(Selection(q)).Lht.HeelCent{k}(:,2)+(i*yAdj);
                    SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                    SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                    if sum(sum(isnan([SrchAreaX, SrchAreaY]))) == 0
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) == 0
                            %                 subplot(1,NumDynTrials, Selection(q));
                            %                 line(a,b,'Color','k');
                            break
                        end
                    end
                end
                if a(1) < DynamicPPTrials(Selection(q)).Lht.HeelCent{k}(1,1)
                    DynamicPPTrials(Selection(q)).L.LatHeelCent{k} = [a b];
                else
                    DynamicPPTrials(Selection(q)).L.MedHeelCent{k} = [a b];
                end
                % search in the other direction
                for i = Loop
                    xAdj = DynamicPPTrials(Selection(q)).Lht.InVHeelCent{k}(:,1);
                    yAdj = DynamicPPTrials(Selection(q)).Lht.InVHeelCent{k}(:,2);
                    a = DynamicPPTrials(Selection(q)).Lht.HeelCent{k}(:,1)-(i*xAdj);
                    b = DynamicPPTrials(Selection(q)).Lht.HeelCent{k}(:,2)-(i*yAdj);
                    SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                    SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                    if sum(sum(isnan([SrchAreaX, SrchAreaY]))) == 0
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) == 0
                            %                 subplot(1,NumDynTrials, Selection(q));
                            %                 line(a,b,'Color','k');
                            break
                        end
                    end
                end
                if a(1) < DynamicPPTrials(Selection(q)).Lht.HeelCent{k}(1,1)
                    DynamicPPTrials(Selection(q)).L.LatHeelCent{k} = [a b];
                else
                    DynamicPPTrials(Selection(q)).L.MedHeelCent{k} = [a b];
                end
                clearvars xAdj yAdj a b SrchAreaX SrchAreaY Zone Region
            end
        end
    end
    
    %% RIGHT
    if DynamicPPTrials(Selection(q)).NumRight > 0
        for k  = 1:DynamicPPTrials(Selection(q)).NumRight
            
            % manual points
            if isfield(DynamicPPTrials(q).Rht, 'Manual') == 1
                if strcmp(PPSettings.MaskChoice, 'Manual') || strcmp(PPSettings.MaskChoice, 'Manual (pink)') || strcmp(PPSettings.MaskChoice, 'Validation')
                    for i = Loop % increment search additions
                        xAdj = DynamicPPTrials(Selection(q)).Rht.ManualInV{k}(:,1);
                        yAdj = DynamicPPTrials(Selection(q)).Rht.ManualInV{k}(:,2);
                        a = DynamicPPTrials(Selection(q)).Rht.Manual{k}(:,1)+(i*xAdj);
                        b = DynamicPPTrials(Selection(q)).Rht.Manual{k}(:,2)+(i*yAdj);
                        SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj]);
                        SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj]);
                        if sum(sum(isnan([SrchAreaX, SrchAreaY]))) == 0
                            Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                            Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                            if sum(sum(Region)) == 0
                                %                         subplot(1,NumDynTrials, Selection(q));
                                %                         line(a,b,'Color','k');
                                % uncomment to check search points
                                %                 subplot(1, NumDynTrials, q);
                                %                 hold on;
                                %                 for ii = 1:4
                                %                     plot(SrchAreaX(ii), SrchAreaY(ii),'k.')
                                %                 end
                                %                 plot(a(:), b(:), 'r')
                                break
                            end
                        end
                    end
                    if a(1) > DynamicPPTrials(Selection(q)).Rht.Manual{k}(1,1)
                        DynamicPPTrials(Selection(q)).R.LatManual{k} = [a b];
                    else
                        DynamicPPTrials(Selection(q)).R.MedManual{k} = [a b];
                    end
                    % search in the other direction
                    for i = Loop
                        xAdj = DynamicPPTrials(Selection(q)).Rht.ManualInV{k}(:,1);
                        yAdj = DynamicPPTrials(Selection(q)).Rht.ManualInV{k}(:,2);
                        a = DynamicPPTrials(Selection(q)).Rht.Manual{k}(:,1)-(i*xAdj);
                        b = DynamicPPTrials(Selection(q)).Rht.Manual{k}(:,2)-(i*yAdj);
                        SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                        SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                        if sum(sum(isnan([SrchAreaX, SrchAreaY]))) == 0
                            Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                            Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                            if sum(sum(Region)) == 0
                                %                         subplot(1,NumDynTrials, Selection(q));
                                %                         line(a,b,'Color','k');
                                % uncomment to check search points
                                %                 subplot(1, NumDynTrials, q);
                                %                 hold on;
                                %                 for ii = 1:4
                                %                     plot(SrchAreaX(ii), SrchAreaY(ii),'k.')
                                %                 end
                                %                 plot(a(:), b(:), 'r')
                                break
                            end
                        end
                    end
                    if a(1) > DynamicPPTrials(Selection(q)).Rht.Manual{k}(1,1)
                        DynamicPPTrials(Selection(q)).R.LatManual{k} = [a b];
                    else
                        DynamicPPTrials(Selection(q)).R.MedManual{k} = [a b];
                    end
                    clearvars xAdj yAdj a b SrchAreaX SrchAreaY Zone Region
                end
            end
            % general image recognition points
            for i = Loop % increment search additions
                xAdj = DynamicPPTrials(Selection(q)).Rht.GeneralInV{k}(:,1);
                yAdj = DynamicPPTrials(Selection(q)).Rht.GeneralInV{k}(:,2);
                a = DynamicPPTrials(Selection(q)).Rht.General{k}(:,1)+(i*xAdj);
                b = DynamicPPTrials(Selection(q)).Rht.General{k}(:,2)+(i*yAdj);
                SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj]);
                SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj]);
                if sum(sum(isnan([SrchAreaX, SrchAreaY]))) == 0
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) == 0
                        %                         subplot(1,NumDynTrials, Selection(q));
                        %                         line(a,b,'Color','k');
                        % uncomment to check search points
                        %                 subplot(1, NumDynTrials, q);
                        %                 hold on;
                        %                 for ii = 1:4
                        %                     plot(SrchAreaX(ii), SrchAreaY(ii),'k.')
                        %                 end
                        %                 plot(a(:), b(:), 'r')
                        break
                    end
                end
            end
            if a(1) > DynamicPPTrials(Selection(q)).Rht.General{k}(1,1)
                DynamicPPTrials(Selection(q)).R.LatGeneral{k} = [a b];
            else
                DynamicPPTrials(Selection(q)).R.MedGeneral{k} = [a b];
            end
            % search in the other direction
            for i = Loop
                xAdj = DynamicPPTrials(Selection(q)).Rht.GeneralInV{k}(:,1);
                yAdj = DynamicPPTrials(Selection(q)).Rht.GeneralInV{k}(:,2);
                a = DynamicPPTrials(Selection(q)).Rht.General{k}(:,1)-(i*xAdj);
                b = DynamicPPTrials(Selection(q)).Rht.General{k}(:,2)-(i*yAdj);
                SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                if sum(sum(isnan([SrchAreaX, SrchAreaY]))) == 0
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) == 0
                        %                         subplot(1,NumDynTrials, Selection(q));
                        %                         line(a,b,'Color','k');
                        % uncomment to check search points
                        %                 subplot(1, NumDynTrials, q);
                        %                 hold on;
                        %                 for ii = 1:4
                        %                     plot(SrchAreaX(ii), SrchAreaY(ii),'k.')
                        %                 end
                        %                 plot(a(:), b(:), 'r')
                        break
                    end
                end
            end
            if a(1) > DynamicPPTrials(Selection(q)).Rht.General{k}(1,1)
                DynamicPPTrials(Selection(q)).R.LatGeneral{k} = [a b];
            else
                DynamicPPTrials(Selection(q)).R.MedGeneral{k} = [a b];
            end
            clearvars xAdj yAdj a b SrchAreaX SrchAreaY Zone Region
            
            % CoP points
            for i = Loop
                xAdj = DynamicPPTrials(Selection(q)).Rht.InVCoP{k}(:,1);
                yAdj = DynamicPPTrials(Selection(q)).Rht.InVCoP{k}(:,2);
                a = DynamicPPTrials(Selection(q)).Rht.CoP{k}(:,1)+(i*xAdj);
                b = DynamicPPTrials(Selection(q)).Rht.CoP{k}(:,2)+(i*yAdj);
                SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                if sum(sum(isnan([SrchAreaX, SrchAreaY]))) == 0
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) == 0
                        %                 subplot(1,NumDynTrials, Selection(q));
                        %                 line(a,b,'Color','k');
                        break
                    end
                end
            end
            if a(1) > DynamicPPTrials(Selection(q)).Rht.CoP{k}(1,1)
                DynamicPPTrials(Selection(q)).R.LatCoP{k} = [a b];
            else
                DynamicPPTrials(Selection(q)).R.MedCoP{k} = [a b];
            end
            % search in the other direction
            for i = Loop
                xAdj = DynamicPPTrials(Selection(q)).Rht.InVCoP{k}(:,1);
                yAdj = DynamicPPTrials(Selection(q)).Rht.InVCoP{k}(:,2);
                a = DynamicPPTrials(Selection(q)).Rht.CoP{k}(:,1)-(i*xAdj);
                b = DynamicPPTrials(Selection(q)).Rht.CoP{k}(:,2)-(i*yAdj);
                SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                if sum(sum(isnan([SrchAreaX, SrchAreaY]))) == 0
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) == 0
                        %                 subplot(1,NumDynTrials, Selection(q));
                        %                 line(a,b,'Color','k');
                        break
                    end
                end
            end
            if a(1) > DynamicPPTrials(Selection(q)).Rht.CoP{k}(1,1)
                DynamicPPTrials(Selection(q)).R.LatCoP{k} = [a b];
            else
                DynamicPPTrials(Selection(q)).R.MedCoP{k} = [a b];
            end
            clearvars xAdj yAdj a b SrchAreaX SrchAreaY Zone Region
            
            % 66% CoP points
            for i = Loop
                xAdj = DynamicPPTrials(Selection(q)).Rht.CoPInV66{k}(:,1);
                yAdj = DynamicPPTrials(Selection(q)).Rht.CoPInV66{k}(:,2);
                a = DynamicPPTrials(Selection(q)).Rht.CoP66{k}(:,1)+(i*xAdj);
                b = DynamicPPTrials(Selection(q)).Rht.CoP66{k}(:,2)+(i*yAdj);
                SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                if sum(sum(isnan([SrchAreaX, SrchAreaY]))) == 0
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) == 0
                        %                 subplot(1,NumDynTrials, Selection(q));
                        %                 line(a,b,'Color','k');
                        break
                    end
                end
            end
            if a(1) > DynamicPPTrials(Selection(q)).Rht.CoP66{k}(1,1)
                DynamicPPTrials(Selection(q)).R.Lat66{k} = [a b];
            else
                DynamicPPTrials(Selection(q)).R.Med66{k} = [a b];
            end
            % search in the other direction
            for i = Loop
                xAdj = DynamicPPTrials(Selection(q)).Rht.CoPInV66{k}(:,1);
                yAdj = DynamicPPTrials(Selection(q)).Rht.CoPInV66{k}(:,2);
                a = DynamicPPTrials(Selection(q)).Rht.CoP66{k}(:,1)-(i*xAdj);
                b = DynamicPPTrials(Selection(q)).Rht.CoP66{k}(:,2)-(i*yAdj);
                SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                if sum(sum(isnan([SrchAreaX, SrchAreaY]))) == 0
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) == 0
                        %                 subplot(1,NumDynTrials, Selection(q));
                        %                 line(a,b,'Color','k');
                        break
                    end
                end
            end
            if a(1) > DynamicPPTrials(Selection(q)).Rht.CoP66{k}(1,1)
                DynamicPPTrials(Selection(q)).R.Lat66{k} = [a b];
            else
                DynamicPPTrials(Selection(q)).R.Med66{k} = [a b];
            end
            clearvars xAdj yAdj a b SrchAreaX SrchAreaY Zone Region
            
            % inter peak points
            if Block.IP == 0
                for i = Loop
                    xAdj = DynamicPPTrials(Selection(q)).Rht.InVIP{k}(:,1);
                    yAdj = DynamicPPTrials(Selection(q)).Rht.InVIP{k}(:,2);
                    a = DynamicPPTrials(Selection(q)).Rht.IP{k}(:,1)+(i*xAdj);
                    b = DynamicPPTrials(Selection(q)).Rht.IP{k}(:,2)+(i*yAdj);
                    SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                    SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                    if sum(sum(isnan([SrchAreaX, SrchAreaY]))) == 0
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) == 0
                            %                 subplot(1,NumDynTrials, Selection(q));
                            %                 line(a,b,'Color','k');
                            break
                        end
                    end
                end
                if a(1) > DynamicPPTrials(Selection(q)).Rht.IP{k}(1,1)
                    DynamicPPTrials(Selection(q)).R.LatIP{k} = [a b];
                else
                    DynamicPPTrials(Selection(q)).R.MedIP{k} = [a b];
                end
                % search in the other direction
                for i = Loop
                    xAdj = DynamicPPTrials(Selection(q)).Rht.InVIP{k}(:,1);
                    yAdj = DynamicPPTrials(Selection(q)).Rht.InVIP{k}(:,2);
                    a = DynamicPPTrials(Selection(q)).Rht.IP{k}(:,1)-(i*xAdj);
                    b = DynamicPPTrials(Selection(q)).Rht.IP{k}(:,2)-(i*yAdj);
                    SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                    SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                    if sum(sum(isnan([SrchAreaX, SrchAreaY]))) == 0
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) == 0
                            %                 subplot(1,NumDynTrials, Selection(q));
                            %                 line(a,b,'Color','k');
                            break
                        end
                    end
                end
                if a(1) > DynamicPPTrials(Selection(q)).Rht.IP{k}(1,1)
                    DynamicPPTrials(Selection(q)).R.LatIP{k} = [a b];
                else
                    DynamicPPTrials(Selection(q)).R.MedIP{k} = [a b];
                end
                clearvars xAdj yAdj a b SrchAreaX SrchAreaY Zone Region
            end
            
            % heel to centroid points
            if Block.HC == 0
                for i = Loop
                    xAdj = DynamicPPTrials(Selection(q)).Rht.InVHeelCent{k}(:,1);
                    yAdj = DynamicPPTrials(Selection(q)).Rht.InVHeelCent{k}(:,2);
                    a = DynamicPPTrials(Selection(q)).Rht.HeelCent{k}(:,1)+(i*xAdj);
                    b = DynamicPPTrials(Selection(q)).Rht.HeelCent{k}(:,2)+(i*yAdj);
                    SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                    SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                    if sum(sum(isnan([SrchAreaX, SrchAreaY]))) == 0
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) == 0
                            %                 subplot(1,NumDynTrials, Selection(q));
                            %                 line(a,b,'Color','k');
                            break
                        end
                    end
                end
                if a(1) > DynamicPPTrials(Selection(q)).Rht.HeelCent{k}(1,1)
                    DynamicPPTrials(Selection(q)).R.LatHeelCent{k} = [a b];
                else
                    DynamicPPTrials(Selection(q)).R.MedHeelCent{k} = [a b];
                end
                % search in the other direction
                for i = Loop
                    xAdj = DynamicPPTrials(Selection(q)).Rht.InVHeelCent{k}(:,1);
                    yAdj = DynamicPPTrials(Selection(q)).Rht.InVHeelCent{k}(:,2);
                    a = DynamicPPTrials(Selection(q)).Rht.HeelCent{k}(:,1)-(i*xAdj);
                    b = DynamicPPTrials(Selection(q)).Rht.HeelCent{k}(:,2)-(i*yAdj);
                    SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                    SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                    if sum(sum(isnan([SrchAreaX, SrchAreaY]))) == 0
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) == 0
                            %                 subplot(1,NumDynTrials, Selection(q));
                            %                 line(a,b,'Color','k');
                            break
                        end
                    end
                end
                if a(1) > DynamicPPTrials(Selection(q)).Rht.HeelCent{k}(1,1)
                    DynamicPPTrials(Selection(q)).R.LatHeelCent{k} = [a b];
                else
                    DynamicPPTrials(Selection(q)).R.MedHeelCent{k} = [a b];
                end
                clearvars xAdj yAdj a b SrchAreaX SrchAreaY Zone Region
            end
        end
    end
end

end