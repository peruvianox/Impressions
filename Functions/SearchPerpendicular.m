function [DynamicPPTrials] = SearchPerpendicular(DynamicPPTrials, Selection, Block, PPSettings)

%% initialize
% set search area
Srch = 2;
% set loop distances
Loop = 0.1:0.1:20;

%% loop through trials to adjust points
for q = 1:length(Selection)
    [m,n] = size(DynamicPPTrials(Selection(q)).MaskLog); % define matrix size
    %% LEFT foot length
    for j = 1: DynamicPPTrials(Selection(q)).NumLeft % loop through left steps
        if DynamicPPTrials(Selection(q)).PoorPressure.Left(j) == 0
            %% manual points
            if isfield(DynamicPPTrials(q).Lht, 'Manual') == 1
                if strcmp(PPSettings.MaskChoice, 'Manual') || strcmp(PPSettings.MaskChoice, 'Manual (pink)') || strcmp(PPSettings.MaskChoice, 'Validation')
                    % determine if the current TOE point is on a pressure
                    xAdj = DynamicPPTrials(Selection(q)).Lht.ManualV{j}(:,1);
                    yAdj = DynamicPPTrials(Selection(q)).Lht.ManualV{j}(:,2);
                    a = [DynamicPPTrials(Selection(q)).L.LatManual{j}(2,1), DynamicPPTrials(Selection(q)).L.MedManual{j}(2,1)];
                    b = [DynamicPPTrials(Selection(q)).L.LatManual{j}(2,2), DynamicPPTrials(Selection(q)).L.MedManual{j}(2,2)];
                    SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                    SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                    %                subplot(1,NumDynTrials, Selection(q)); hold on; plot(SrchAreaX, SrchAreaY, '.g');
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) > 0
                        PtOnOff = 'On';
                    else
                        PtOnOff = 'Off';
                    end
                    if strcmp(PtOnOff,'On') == 1 % if the FP points are within the border of the foot, search from the toe forward
                        for i = Loop % increment search points out from toe
                            a = [DynamicPPTrials(Selection(q)).L.LatManual{j}(2,1) + (i*xAdj), DynamicPPTrials(Selection(q)).L.MedManual{j}(2,1) + (i*xAdj)];
                            b = [DynamicPPTrials(Selection(q)).L.LatManual{j}(2,2) + (i*yAdj), DynamicPPTrials(Selection(q)).L.MedManual{j}(2,2) + (i*yAdj)];
                            SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                            SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                            Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                            Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                            if sum(sum(Region)) == 0
                                %                         subplot(1,NumDynTrials, Selection(q));
                                %                         line(a,b,'Color','k');
                                DynamicPPTrials(Selection(q)).L.LatManual{j}(2,:) = [a(1), b(1)];
                                DynamicPPTrials(Selection(q)).L.MedManual{j}(2,:) = [a(2), b(2)];
                                % uncomment to check search points
                                %                 subplot(1, NumDynTrials, q);  hold on;
                                %                 for ii = 1:4
                                %                     plot(SrchAreaX(ii), SrchAreaY(ii),'k.')
                                %                 end
                                %                 plot(a(:), b(:), 'r')
                                break
                            end
                        end
                    else % if the FP points are outside the border of the foot, search from the toe back towards the heel
                        for i = Loop % search forward to make sure there arent more pressures forward
                            a = [DynamicPPTrials(Selection(q)).L.LatManual{j}(2,1) + (i*xAdj), DynamicPPTrials(Selection(q)).L.MedManual{j}(2,1) + (i*xAdj)];
                            b = [DynamicPPTrials(Selection(q)).L.LatManual{j}(2,2) + (i*yAdj), DynamicPPTrials(Selection(q)).L.MedManual{j}(2,2) + (i*yAdj)];
                            SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                            SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                            Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                            Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                            if sum(sum(Region)) > 0
                                PressForward = 'Yes';
                                ForwardRegion = [a b];
                                break
                            else
                                PressForward = 'No';
                            end
                        end
                        if strcmp(PressForward,'No') == 1
                            for i = Loop % increment search points back from toe
                                a = [DynamicPPTrials(Selection(q)).L.LatManual{j}(2,1) - (i*xAdj), DynamicPPTrials(Selection(q)).L.MedManual{j}(2,1) - (i*xAdj)];
                                b = [DynamicPPTrials(Selection(q)).L.LatManual{j}(2,2) - (i*yAdj), DynamicPPTrials(Selection(q)).L.MedManual{j}(2,2) - (i*yAdj)];
                                SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                                SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                                Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                                Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                                if sum(sum(Region)) > 0
                                    %                             subplot(1,NumDynTrials, Selection(q));
                                    %                             line(a,b,'Color','k');
                                    DynamicPPTrials(Selection(q)).L.LatManual{j}(2,:) = [a(1), b(1)];
                                    DynamicPPTrials(Selection(q)).L.MedManual{j}(2,:) = [a(2), b(2)];
                                    % uncomment to check search points
                                    %                                      subplot(1, NumDynTrials, q);  hold on;
                                    %                                      for ii = 1:4
                                    %                                          plot(SrchAreaX(ii), SrchAreaY(ii),'k.')
                                    %                                      end
                                    %                                      plot(a(:), b(:), 'r')
                                    break
                                end
                            end
                        else % pressures are in front of the original toe point, keep searching forward
                            for i = Loop
                                a = [ForwardRegion(1) + (i*xAdj), ForwardRegion(2) + (i*xAdj)];
                                b = [ForwardRegion(3) + (i*yAdj), ForwardRegion(4) + (i*yAdj)];
                                SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                                SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                                Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                                Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                                if sum(sum(Region)) > 0
                                    %                             subplot(1,NumDynTrials, Selection(q));
                                    %                             line(a,b,'Color','k');
                                    DynamicPPTrials(Selection(q)).L.LatManual{j}(2,:) = [a(1), b(1)];
                                    DynamicPPTrials(Selection(q)).L.MedManual{j}(2,:) = [a(2), b(2)];
                                    % uncomment to check search points
                                    %                 subplot(1, NumDynTrials, q);  hold on;
                                    %                 for ii = 1:4
                                    %                     plot(SrchAreaX(ii), SrchAreaY(ii),'k.')
                                    %                 end
                                    %                 plot(a(:), b(:), 'r')
                                    break
                                end
                            end
                        end
                    end
                    
                    % determine if the current HEEL point is on a pressure
                    a = [DynamicPPTrials(Selection(q)).L.LatManual{j}(1,1), DynamicPPTrials(Selection(q)).L.MedManual{j}(1,1)];
                    b = [DynamicPPTrials(Selection(q)).L.LatManual{j}(1,2), DynamicPPTrials(Selection(q)).L.MedManual{j}(1,2)];
                    SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                    SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) > 0
                        PtOnOff = 'On';
                    else
                        PtOnOff = 'Off';
                    end
                    if strcmp(PtOnOff,'On') == 1 % if the FP points are within the border of the foot, search from the Heel backward
                        for i = Loop % increment search points out from heel
                            a = [DynamicPPTrials(Selection(q)).L.LatManual{j}(1,1) - (i*xAdj), DynamicPPTrials(Selection(q)).L.MedManual{j}(1,1) - (i*xAdj)];
                            b = [DynamicPPTrials(Selection(q)).L.LatManual{j}(1,2) - (i*yAdj), DynamicPPTrials(Selection(q)).L.MedManual{j}(1,2) - (i*yAdj)];
                            SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                            SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                            Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                            Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                            if sum(sum(Region)) == 0
                                %                         subplot(1,NumDynTrials, Selection(q));
                                %                         line(a,b,'Color','k');
                                DynamicPPTrials(Selection(q)).L.LatManual{j}(1,:) = [a(1), b(1)];
                                DynamicPPTrials(Selection(q)).L.MedManual{j}(1,:) = [a(2), b(2)];
                                % uncomment to check search points
                                %                 subplot(1, NumDynTrials, q); hold on;
                                %                 for ii = 1:4
                                %                     plot(SrchAreaX(ii), SrchAreaY(ii),'k.')
                                %                 end
                                %                 plot(a(:), b(:), 'r')
                                break
                            end
                        end
                    else % if the FP points are outside the border of the foot, search from the heel towards the toe
                        for i = Loop % increment search points forward from heel
                            a = [DynamicPPTrials(Selection(q)).L.LatManual{j}(1,1) + (i*xAdj), DynamicPPTrials(Selection(q)).L.MedManual{j}(1,1) + (i*xAdj)];
                            b = [DynamicPPTrials(Selection(q)).L.LatManual{j}(1,2) + (i*yAdj), DynamicPPTrials(Selection(q)).L.MedManual{j}(1,2) + (i*yAdj)];
                            SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                            SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                            Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                            Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                            if sum(sum(Region)) > 0
                                %                         subplot(1,NumDynTrials, Selection(q));
                                %                         line(a,b,'Color','k');
                                DynamicPPTrials(Selection(q)).L.LatManual{j}(1,:) = [a(1), b(1)];
                                DynamicPPTrials(Selection(q)).L.MedManual{j}(1,:) = [a(2), b(2)];
                                % uncomment to check search points
                                %                                      subplot(1, NumDynTrials, q); hold on;
                                %                                      for ii = 1:4
                                %                                          plot(SrchAreaX(ii), SrchAreaY(ii),'k.')
                                %                                      end
                                %                                      plot(a(:), b(:), 'r')
                                break
                            end
                        end
                    end
                end
            end
            
            %% general image processing points
            % determine if the current TOE point is on a pressure
            xAdj = DynamicPPTrials(Selection(q)).Lht.GeneralV{j}(:,1);
            yAdj = DynamicPPTrials(Selection(q)).Lht.GeneralV{j}(:,2);
            a = [DynamicPPTrials(Selection(q)).L.LatGeneral{j}(2,1), DynamicPPTrials(Selection(q)).L.MedGeneral{j}(2,1)];
            b = [DynamicPPTrials(Selection(q)).L.LatGeneral{j}(2,2), DynamicPPTrials(Selection(q)).L.MedGeneral{j}(2,2)];
            SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
            SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
            %                subplot(1,NumDynTrials, Selection(q)); hold on; plot(SrchAreaX, SrchAreaY, '.g');
            Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
            Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
            if sum(sum(Region)) > 0
                PtOnOff = 'On';
            else
                PtOnOff = 'Off';
            end
            if strcmp(PtOnOff,'On') == 1 % if the FP points are within the border of the foot, search from the toe forward
                for i = Loop % increment search points out from toe
                    a = [DynamicPPTrials(Selection(q)).L.LatGeneral{j}(2,1) + (i*xAdj), DynamicPPTrials(Selection(q)).L.MedGeneral{j}(2,1) + (i*xAdj)];
                    b = [DynamicPPTrials(Selection(q)).L.LatGeneral{j}(2,2) + (i*yAdj), DynamicPPTrials(Selection(q)).L.MedGeneral{j}(2,2) + (i*yAdj)];
                    SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                    SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) == 0
                        %                         subplot(1,NumDynTrials, Selection(q));
                        %                         line(a,b,'Color','k');
                        DynamicPPTrials(Selection(q)).L.LatGeneral{j}(2,:) = [a(1), b(1)];
                        DynamicPPTrials(Selection(q)).L.MedGeneral{j}(2,:) = [a(2), b(2)];
                        % uncomment to check search points
                        %                 subplot(1, NumDynTrials, q);  hold on;
                        %                 for ii = 1:4
                        %                     plot(SrchAreaX(ii), SrchAreaY(ii),'k.')
                        %                 end
                        %                 plot(a(:), b(:), 'r')
                        break
                    end
                end
            else % if the FP points are outside the border of the foot, search from the toe back towards the heel
                for i = Loop % search forward to make sure there arent more pressures forward
                    a = [DynamicPPTrials(Selection(q)).L.LatGeneral{j}(2,1) + (i*xAdj), DynamicPPTrials(Selection(q)).L.MedGeneral{j}(2,1) + (i*xAdj)];
                    b = [DynamicPPTrials(Selection(q)).L.LatGeneral{j}(2,2) + (i*yAdj), DynamicPPTrials(Selection(q)).L.MedGeneral{j}(2,2) + (i*yAdj)];
                    SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                    SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) > 0
                        PressForward = 'Yes';
                        ForwardRegion = [a b];
                        break
                    else
                        PressForward = 'No';
                    end
                end
                if strcmp(PressForward,'No') == 1
                    for i = Loop % increment search points back from toe
                        a = [DynamicPPTrials(Selection(q)).L.LatGeneral{j}(2,1) - (i*xAdj), DynamicPPTrials(Selection(q)).L.MedGeneral{j}(2,1) - (i*xAdj)];
                        b = [DynamicPPTrials(Selection(q)).L.LatGeneral{j}(2,2) - (i*yAdj), DynamicPPTrials(Selection(q)).L.MedGeneral{j}(2,2) - (i*yAdj)];
                        SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                        SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) > 0
                            %                             subplot(1,NumDynTrials, Selection(q));
                            %                             line(a,b,'Color','k');
                            DynamicPPTrials(Selection(q)).L.LatGeneral{j}(2,:) = [a(1), b(1)];
                            DynamicPPTrials(Selection(q)).L.MedGeneral{j}(2,:) = [a(2), b(2)];
                            % uncomment to check search points
                            %                                      subplot(1, NumDynTrials, q);  hold on;
                            %                                      for ii = 1:4
                            %                                          plot(SrchAreaX(ii), SrchAreaY(ii),'k.')
                            %                                      end
                            %                                      plot(a(:), b(:), 'r')
                            break
                        end
                    end
                else % pressures are in front of the original toe point, keep searching forward
                    for i = Loop
                        a = [ForwardRegion(1) + (i*xAdj), ForwardRegion(2) + (i*xAdj)];
                        b = [ForwardRegion(3) + (i*yAdj), ForwardRegion(4) + (i*yAdj)];
                        SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                        SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) > 0
                            %                             subplot(1,NumDynTrials, Selection(q));
                            %                             line(a,b,'Color','k');
                            DynamicPPTrials(Selection(q)).L.LatGeneral{j}(2,:) = [a(1), b(1)];
                            DynamicPPTrials(Selection(q)).L.MedGeneral{j}(2,:) = [a(2), b(2)];
                            % uncomment to check search points
                            %                 subplot(1, NumDynTrials, q);  hold on;
                            %                 for ii = 1:4
                            %                     plot(SrchAreaX(ii), SrchAreaY(ii),'k.')
                            %                 end
                            %                 plot(a(:), b(:), 'r')
                            break
                        end
                    end
                end
            end
            
            % determine if the current HEEL point is on a pressure
            a = [DynamicPPTrials(Selection(q)).L.LatGeneral{j}(1,1), DynamicPPTrials(Selection(q)).L.MedGeneral{j}(1,1)];
            b = [DynamicPPTrials(Selection(q)).L.LatGeneral{j}(1,2), DynamicPPTrials(Selection(q)).L.MedGeneral{j}(1,2)];
            SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
            SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
            Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
            Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
            if sum(sum(Region)) > 0
                PtOnOff = 'On';
            else
                PtOnOff = 'Off';
            end
            if strcmp(PtOnOff,'On') == 1 % if the FP points are within the border of the foot, search from the Heel backward
                for i = Loop % increment search points out from heel
                    a = [DynamicPPTrials(Selection(q)).L.LatGeneral{j}(1,1) - (i*xAdj), DynamicPPTrials(Selection(q)).L.MedGeneral{j}(1,1) - (i*xAdj)];
                    b = [DynamicPPTrials(Selection(q)).L.LatGeneral{j}(1,2) - (i*yAdj), DynamicPPTrials(Selection(q)).L.MedGeneral{j}(1,2) - (i*yAdj)];
                    SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                    SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) == 0
                        %                         subplot(1,NumDynTrials, Selection(q));
                        %                         line(a,b,'Color','k');
                        DynamicPPTrials(Selection(q)).L.LatGeneral{j}(1,:) = [a(1), b(1)];
                        DynamicPPTrials(Selection(q)).L.MedGeneral{j}(1,:) = [a(2), b(2)];
                        % uncomment to check search points
                        %                 subplot(1, NumDynTrials, q); hold on;
                        %                 for ii = 1:4
                        %                     plot(SrchAreaX(ii), SrchAreaY(ii),'k.')
                        %                 end
                        %                 plot(a(:), b(:), 'r')
                        break
                    end
                end
            else % if the FP points are outside the border of the foot, search from the heel towards the toe
                for i = Loop % increment search points forward from heel
                    a = [DynamicPPTrials(Selection(q)).L.LatGeneral{j}(1,1) + (i*xAdj), DynamicPPTrials(Selection(q)).L.MedGeneral{j}(1,1) + (i*xAdj)];
                    b = [DynamicPPTrials(Selection(q)).L.LatGeneral{j}(1,2) + (i*yAdj), DynamicPPTrials(Selection(q)).L.MedGeneral{j}(1,2) + (i*yAdj)];
                    SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                    SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) > 0
                        %                         subplot(1,NumDynTrials, Selection(q));
                        %                         line(a,b,'Color','k');
                        DynamicPPTrials(Selection(q)).L.LatGeneral{j}(1,:) = [a(1), b(1)];
                        DynamicPPTrials(Selection(q)).L.MedGeneral{j}(1,:) = [a(2), b(2)];
                        % uncomment to check search points
                        %                                      subplot(1, NumDynTrials, q); hold on;
                        %                                      for ii = 1:4
                        %                                          plot(SrchAreaX(ii), SrchAreaY(ii),'k.')
                        %                                      end
                        %                                      plot(a(:), b(:), 'r')
                        break
                    end
                end
            end
            
            
            %% CoP points
            % determine if the current TOE point is on a pressure
            xAdj = DynamicPPTrials(Selection(q)).Lht.VCoP{j}(:,1);
            yAdj = DynamicPPTrials(Selection(q)).Lht.VCoP{j}(:,2);
            a = [DynamicPPTrials(Selection(q)).L.LatCoP{j}(2,1), DynamicPPTrials(Selection(q)).L.MedCoP{j}(2,1)];
            b = [DynamicPPTrials(Selection(q)).L.LatCoP{j}(2,2), DynamicPPTrials(Selection(q)).L.MedCoP{j}(2,2)];
            SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
            SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
            Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
            Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
            if sum(sum(Region)) > 0
                PtOnOff = 'On';
            else
                PtOnOff = 'Off';
            end
            if strcmp(PtOnOff,'On') == 1 % if the FP points are within the border of the foot, search from the toe forward
                for i = Loop % increment search points out from toe
                    a = [DynamicPPTrials(Selection(q)).L.LatCoP{j}(2,1) + (i*xAdj), DynamicPPTrials(Selection(q)).L.MedCoP{j}(2,1) + (i*xAdj)];
                    b = [DynamicPPTrials(Selection(q)).L.LatCoP{j}(2,2) + (i*yAdj), DynamicPPTrials(Selection(q)).L.MedCoP{j}(2,2) + (i*yAdj)];
                    SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                    SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) == 0
                        %                     subplot(1,NumDynTrials, Selection(q));
                        %                     line(a,b,'Color','k');
                        DynamicPPTrials(Selection(q)).L.LatCoP{j}(2,:) = [a(1), b(1)];
                        DynamicPPTrials(Selection(q)).L.MedCoP{j}(2,:) = [a(2), b(2)];
                        break
                    end
                end
            else % if the FP points are outside the border of the foot, search from the toe back towards the heel
                for i = Loop % search forward to make sure there arent more pressures forward
                    a = [DynamicPPTrials(Selection(q)).L.LatCoP{j}(2,1) + (i*xAdj), DynamicPPTrials(Selection(q)).L.MedCoP{j}(2,1) + (i*xAdj)];
                    b = [DynamicPPTrials(Selection(q)).L.LatCoP{j}(2,2) + (i*yAdj), DynamicPPTrials(Selection(q)).L.MedCoP{j}(2,2) + (i*yAdj)];
                    SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                    SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) > 0
                        PressForward = 'Yes';
                        ForwardRegion = [a b];
                        break
                    else
                        PressForward = 'No';
                    end
                end
                if strcmp(PressForward,'No') == 1
                    for i = Loop % increment search points back from toe
                        a = [DynamicPPTrials(Selection(q)).L.LatCoP{j}(2,1) - (i*xAdj), DynamicPPTrials(Selection(q)).L.MedCoP{j}(2,1) - (i*xAdj)];
                        b = [DynamicPPTrials(Selection(q)).L.LatCoP{j}(2,2) - (i*yAdj), DynamicPPTrials(Selection(q)).L.MedCoP{j}(2,2) - (i*yAdj)];
                        SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                        SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) > 0
                            %                         subplot(1,NumDynTrials, Selection(q));
                            %                         line(a,b,'Color','k');
                            DynamicPPTrials(Selection(q)).L.LatCoP{j}(2,:) = [a(1), b(1)];
                            DynamicPPTrials(Selection(q)).L.MedCoP{j}(2,:) = [a(2), b(2)];
                            break
                        end
                    end
                else % pressures are in front of the original toe point, keep searching forward
                    for i = Loop
                        a = [ForwardRegion(1) + (i*xAdj), ForwardRegion(2) + (i*xAdj)];
                        b = [ForwardRegion(3) + (i*yAdj), ForwardRegion(4) + (i*yAdj)];
                        SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                        SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) > 0
                            %                         subplot(1,NumDynTrials, Selection(q));
                            %                         line(a,b,'Color','k');
                            DynamicPPTrials(Selection(q)).L.LatCoP{j}(2,:) = [a(1), b(1)];
                            DynamicPPTrials(Selection(q)).L.MedCoP{j}(2,:) = [a(2), b(2)];
                            break
                        end
                    end
                end
            end
            % determine if the current HEEL point is on a pressure
            a = [DynamicPPTrials(Selection(q)).L.LatCoP{j}(1,1), DynamicPPTrials(Selection(q)).L.MedCoP{j}(1,1)];
            b = [DynamicPPTrials(Selection(q)).L.LatCoP{j}(1,2), DynamicPPTrials(Selection(q)).L.MedCoP{j}(1,2)];
            SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
            SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
            Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
            Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
            if sum(sum(Region)) > 0
                PtOnOff = 'On';
            else
                PtOnOff = 'Off';
            end
            if strcmp(PtOnOff,'On') == 1 % if the FP points are within the border of the foot, search from the Heel backward
                for i = Loop % increment search points out from heel
                    a = [DynamicPPTrials(Selection(q)).L.LatCoP{j}(1,1) - (i*xAdj), DynamicPPTrials(Selection(q)).L.MedCoP{j}(1,1) - (i*xAdj)];
                    b = [DynamicPPTrials(Selection(q)).L.LatCoP{j}(1,2) - (i*yAdj), DynamicPPTrials(Selection(q)).L.MedCoP{j}(1,2) - (i*yAdj)];
                    SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                    SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) == 0
                        %                     subplot(1,NumDynTrials, Selection(q));
                        %                     line(a,b,'Color','k');
                        DynamicPPTrials(Selection(q)).L.LatCoP{j}(1,:) = [a(1), b(1)];
                        DynamicPPTrials(Selection(q)).L.MedCoP{j}(1,:) = [a(2), b(2)];
                        break
                    end
                end
            else % if the FP points are outside the border of the foot, search from the heel towards the toe
                for i = Loop % increment search points forward from heel
                    a = [DynamicPPTrials(Selection(q)).L.LatCoP{j}(1,1) + (i*xAdj), DynamicPPTrials(Selection(q)).L.MedCoP{j}(1,1) + (i*xAdj)];
                    b = [DynamicPPTrials(Selection(q)).L.LatCoP{j}(1,2) + (i*yAdj), DynamicPPTrials(Selection(q)).L.MedCoP{j}(1,2) + (i*yAdj)];
                    SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                    SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) > 0
                        %                     subplot(1,NumDynTrials, Selection(q));
                        %                     line(a,b,'Color','k');
                        DynamicPPTrials(Selection(q)).L.LatCoP{j}(1,:) = [a(1), b(1)];
                        DynamicPPTrials(Selection(q)).L.MedCoP{j}(1,:) = [a(2), b(2)];
                        break
                    end
                end
            end
            
            %% 66% CoP points
            % determine if the current TOE point is on a pressure
            xAdj = DynamicPPTrials(Selection(q)).Lht.CoPV66{j}(:,1);
            yAdj = DynamicPPTrials(Selection(q)).Lht.CoPV66{j}(:,2);
            a = [DynamicPPTrials(Selection(q)).L.Lat66{j}(2,1), DynamicPPTrials(Selection(q)).L.Med66{j}(2,1)];
            b = [DynamicPPTrials(Selection(q)).L.Lat66{j}(2,2), DynamicPPTrials(Selection(q)).L.Med66{j}(2,2)];
            SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
            SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
            Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
            Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
            if sum(sum(Region)) > 0
                PtOnOff = 'On';
            else
                PtOnOff = 'Off';
            end
            if strcmp(PtOnOff,'On') == 1 % if the FP points are within the border of the foot, search from the toe forward
                for i = Loop % increment search points out from toe
                    a = [DynamicPPTrials(Selection(q)).L.Lat66{j}(2,1) + (i*xAdj), DynamicPPTrials(Selection(q)).L.Med66{j}(2,1) + (i*xAdj)];
                    b = [DynamicPPTrials(Selection(q)).L.Lat66{j}(2,2) + (i*yAdj), DynamicPPTrials(Selection(q)).L.Med66{j}(2,2) + (i*yAdj)];
                    SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                    SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) == 0
                        %                     subplot(1,NumDynTrials, Selection(q));
                        %                     line(a,b,'Color','k');
                        DynamicPPTrials(Selection(q)).L.Lat66{j}(2,:) = [a(1), b(1)];
                        DynamicPPTrials(Selection(q)).L.Med66{j}(2,:) = [a(2), b(2)];
                        break
                    end
                end
            else % if the FP points are outside the border of the foot, search from the toe back towards the heel
                for i = Loop % search forward to make sure there arent more pressures forward
                    a = [DynamicPPTrials(Selection(q)).L.Lat66{j}(2,1) + (i*xAdj), DynamicPPTrials(Selection(q)).L.Med66{j}(2,1) + (i*xAdj)];
                    b = [DynamicPPTrials(Selection(q)).L.Lat66{j}(2,2) + (i*yAdj), DynamicPPTrials(Selection(q)).L.Med66{j}(2,2) + (i*yAdj)];
                    SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                    SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) > 0
                        PressForward = 'Yes';
                        ForwardRegion = [a b];
                        break
                    else
                        PressForward = 'No';
                    end
                end
                if strcmp(PressForward,'No') == 1
                    for i = Loop % increment search points back from toe
                        a = [DynamicPPTrials(Selection(q)).L.Lat66{j}(2,1) - (i*xAdj), DynamicPPTrials(Selection(q)).L.Med66{j}(2,1) - (i*xAdj)];
                        b = [DynamicPPTrials(Selection(q)).L.Lat66{j}(2,2) - (i*yAdj), DynamicPPTrials(Selection(q)).L.Med66{j}(2,2) - (i*yAdj)];
                        SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                        SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) > 0
                            %                         subplot(1,NumDynTrials, Selection(q));
                            %                         line(a,b,'Color','k');
                            DynamicPPTrials(Selection(q)).L.Lat66{j}(2,:) = [a(1), b(1)];
                            DynamicPPTrials(Selection(q)).L.Med66{j}(2,:) = [a(2), b(2)];
                            break
                        end
                    end
                else % pressures are in front of the original toe point, keep searching forward
                    for i = Loop
                        a = [ForwardRegion(1) + (i*xAdj), ForwardRegion(2) + (i*xAdj)];
                        b = [ForwardRegion(3) + (i*yAdj), ForwardRegion(4) + (i*yAdj)];
                        SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                        SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) > 0
                            %                         subplot(1,NumDynTrials, Selection(q));
                            %                         line(a,b,'Color','k');
                            DynamicPPTrials(Selection(q)).L.Lat66{j}(2,:) = [a(1), b(1)];
                            DynamicPPTrials(Selection(q)).L.Med66{j}(2,:) = [a(2), b(2)];
                            break
                        end
                    end
                end
            end
            % determine if the current HEEL point is on a pressure
            a = [DynamicPPTrials(Selection(q)).L.Lat66{j}(1,1), DynamicPPTrials(Selection(q)).L.Med66{j}(1,1)];
            b = [DynamicPPTrials(Selection(q)).L.Lat66{j}(1,2), DynamicPPTrials(Selection(q)).L.Med66{j}(1,2)];
            SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
            SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
            Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
            Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
            if sum(sum(Region)) > 0
                PtOnOff = 'On';
            else
                PtOnOff = 'Off';
            end
            if strcmp(PtOnOff,'On') == 1 % if the FP points are within the border of the foot, search from the Heel backward
                for i = Loop % increment search points out from heel
                    a = [DynamicPPTrials(Selection(q)).L.Lat66{j}(1,1) - (i*xAdj), DynamicPPTrials(Selection(q)).L.Med66{j}(1,1) - (i*xAdj)];
                    b = [DynamicPPTrials(Selection(q)).L.Lat66{j}(1,2) - (i*yAdj), DynamicPPTrials(Selection(q)).L.Med66{j}(1,2) - (i*yAdj)];
                    SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                    SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) == 0
                        %                     subplot(1,NumDynTrials, Selection(q));
                        %                     line(a,b,'Color','k');
                        DynamicPPTrials(Selection(q)).L.Lat66{j}(1,:) = [a(1), b(1)];
                        DynamicPPTrials(Selection(q)).L.Med66{j}(1,:) = [a(2), b(2)];
                        break
                    end
                end
            else % if the FP points are outside the border of the foot, search from the heel towards the toe
                for i = Loop % increment search points forward from heel
                    a = [DynamicPPTrials(Selection(q)).L.Lat66{j}(1,1) + (i*xAdj), DynamicPPTrials(Selection(q)).L.Med66{j}(1,1) + (i*xAdj)];
                    b = [DynamicPPTrials(Selection(q)).L.Lat66{j}(1,2) + (i*yAdj), DynamicPPTrials(Selection(q)).L.Med66{j}(1,2) + (i*yAdj)];
                    SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                    SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) > 0
                        %                     subplot(1,NumDynTrials, Selection(q));
                        %                     line(a,b,'Color','k');
                        DynamicPPTrials(Selection(q)).L.Lat66{j}(1,:) = [a(1), b(1)];
                        DynamicPPTrials(Selection(q)).L.Med66{j}(1,:) = [a(2), b(2)];
                        break
                    end
                end
            end
            
            %% IP points
            if Block.IP == 0
                % determine if the current TOE point is on a pressure
                xAdj = DynamicPPTrials(Selection(q)).Lht.VIP{j}(:,1);
                yAdj = DynamicPPTrials(Selection(q)).Lht.VIP{j}(:,2);
                if isempty(DynamicPPTrials(Selection(q)).L.LatIP)
                    break
                end
                a = [DynamicPPTrials(Selection(q)).L.LatIP{j}(2,1), DynamicPPTrials(Selection(q)).L.MedIP{j}(2,1)];
                b = [DynamicPPTrials(Selection(q)).L.LatIP{j}(2,2), DynamicPPTrials(Selection(q)).L.MedIP{j}(2,2)];
                SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                if sum(sum(Region)) > 0
                    PtOnOff = 'On';
                else
                    PtOnOff = 'Off';
                end
                if strcmp(PtOnOff,'On') == 1 % if the FP points are within the border of the foot, search from the toe forward
                    for i = Loop % increment search points out from toe
                        a = [DynamicPPTrials(Selection(q)).L.LatIP{j}(2,1) + (i*xAdj), DynamicPPTrials(Selection(q)).L.MedIP{j}(2,1) + (i*xAdj)];
                        b = [DynamicPPTrials(Selection(q)).L.LatIP{j}(2,2) + (i*yAdj), DynamicPPTrials(Selection(q)).L.MedIP{j}(2,2) + (i*yAdj)];
                        SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                        SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) == 0
                            %                     subplot(1,NumDynTrials, Selection(q));
                            %                     line(a,b,'Color','k');
                            DynamicPPTrials(Selection(q)).L.LatIP{j}(2,:) = [a(1), b(1)];
                            DynamicPPTrials(Selection(q)).L.MedIP{j}(2,:) = [a(2), b(2)];
                            break
                        end
                    end
                else % if the FP points are outside the border of the foot, search from the toe back towards the heel
                    for i = Loop % search forward to make sure there arent more pressures forward
                        a = [DynamicPPTrials(Selection(q)).L.LatIP{j}(2,1) + (i*xAdj), DynamicPPTrials(Selection(q)).L.MedIP{j}(2,1) + (i*xAdj)];
                        b = [DynamicPPTrials(Selection(q)).L.LatIP{j}(2,2) + (i*yAdj), DynamicPPTrials(Selection(q)).L.MedIP{j}(2,2) + (i*yAdj)];
                        SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                        SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) > 0
                            PressForward = 'Yes';
                            ForwardRegion = [a b];
                            break
                        else
                            PressForward = 'No';
                        end
                    end
                    if strcmp(PressForward,'No') == 1
                        for i = Loop % increment search points back from toe
                            a = [DynamicPPTrials(Selection(q)).L.LatIP{j}(2,1) - (i*xAdj), DynamicPPTrials(Selection(q)).L.MedIP{j}(2,1) - (i*xAdj)];
                            b = [DynamicPPTrials(Selection(q)).L.LatIP{j}(2,2) - (i*yAdj), DynamicPPTrials(Selection(q)).L.MedIP{j}(2,2) - (i*yAdj)];
                            SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                            SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                            Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                            Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                            if sum(sum(Region)) > 0
                                %                         subplot(1,NumDynTrials, Selection(q));
                                %                         line(a,b,'Color','k');
                                DynamicPPTrials(Selection(q)).L.LatIP{j}(2,:) = [a(1), b(1)];
                                DynamicPPTrials(Selection(q)).L.MedIP{j}(2,:) = [a(2), b(2)];
                                break
                            end
                        end
                    else % pressures are in front of the original toe point, keep searching forward
                        for i = Loop
                            a = [ForwardRegion(1) + (i*xAdj), ForwardRegion(2) + (i*xAdj)];
                            b = [ForwardRegion(3) + (i*yAdj), ForwardRegion(4) + (i*yAdj)];
                            SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                            SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                            Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                            Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                            if sum(sum(Region)) > 0
                                %                         subplot(1,NumDynTrials, Selection(q));
                                %                         line(a,b,'Color','k');
                                DynamicPPTrials(Selection(q)).L.LatIP{j}(2,:) = [a(1), b(1)];
                                DynamicPPTrials(Selection(q)).L.MedIP{j}(2,:) = [a(2), b(2)];
                                break
                            end
                        end
                    end
                end
                % determine if the current HEEL point is on a pressure
                a = [DynamicPPTrials(Selection(q)).L.LatIP{j}(1,1), DynamicPPTrials(Selection(q)).L.MedIP{j}(1,1)];
                b = [DynamicPPTrials(Selection(q)).L.LatIP{j}(1,2), DynamicPPTrials(Selection(q)).L.MedIP{j}(1,2)];
                SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                if sum(sum(Region)) > 0
                    PtOnOff = 'On';
                else
                    PtOnOff = 'Off';
                end
                if strcmp(PtOnOff,'On') == 1 % if the FP points are within the border of the foot, search from the Heel backward
                    for i = Loop % increment search points out from heel
                        a = [DynamicPPTrials(Selection(q)).L.LatIP{j}(1,1) - (i*xAdj), DynamicPPTrials(Selection(q)).L.MedIP{j}(1,1) - (i*xAdj)];
                        b = [DynamicPPTrials(Selection(q)).L.LatIP{j}(1,2) - (i*yAdj), DynamicPPTrials(Selection(q)).L.MedIP{j}(1,2) - (i*yAdj)];
                        SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                        SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) == 0
                            %                     subplot(1,NumDynTrials, Selection(q));
                            %                     line(a,b,'Color','k');
                            DynamicPPTrials(Selection(q)).L.LatIP{j}(1,:) = [a(1), b(1)];
                            DynamicPPTrials(Selection(q)).L.MedIP{j}(1,:) = [a(2), b(2)];
                            break
                        end
                    end
                else % if the FP points are outside the border of the foot, search from the heel towards the toe
                    for i = Loop % increment search points forward from heel
                        a = [DynamicPPTrials(Selection(q)).L.LatIP{j}(1,1) + (i*xAdj), DynamicPPTrials(Selection(q)).L.MedIP{j}(1,1) + (i*xAdj)];
                        b = [DynamicPPTrials(Selection(q)).L.LatIP{j}(1,2) + (i*yAdj), DynamicPPTrials(Selection(q)).L.MedIP{j}(1,2) + (i*yAdj)];
                        SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                        SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) > 0
                            %                     subplot(1,NumDynTrials, Selection(q));
                            %                     line(a,b,'Color','k');
                            DynamicPPTrials(Selection(q)).L.LatIP{j}(1,:) = [a(1), b(1)];
                            DynamicPPTrials(Selection(q)).L.MedIP{j}(1,:) = [a(2), b(2)];
                            break
                        end
                    end
                end
            end
            
            %% Heel to centroid points
            if Block.HC == 0
                % determine if the current TOE point is on a pressure
                xAdj = DynamicPPTrials(Selection(q)).Lht.VHeelCent{j}(:,1);
                yAdj = DynamicPPTrials(Selection(q)).Lht.VHeelCent{j}(:,2);
                a = [DynamicPPTrials(Selection(q)).L.LatHeelCent{j}(2,1), DynamicPPTrials(Selection(q)).L.MedHeelCent{j}(2,1)];
                b = [DynamicPPTrials(Selection(q)).L.LatHeelCent{j}(2,2), DynamicPPTrials(Selection(q)).L.MedHeelCent{j}(2,2)];
                SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                if sum(sum(Region)) > 0
                    PtOnOff = 'On';
                else
                    PtOnOff = 'Off';
                end
                if strcmp(PtOnOff,'On') == 1 % if the FP points are within the border of the foot, search from the toe forward
                    for i = Loop % increment search points out from toe
                        a = [DynamicPPTrials(Selection(q)).L.LatHeelCent{j}(2,1) + (i*xAdj), DynamicPPTrials(Selection(q)).L.MedHeelCent{j}(2,1) + (i*xAdj)];
                        b = [DynamicPPTrials(Selection(q)).L.LatHeelCent{j}(2,2) + (i*yAdj), DynamicPPTrials(Selection(q)).L.MedHeelCent{j}(2,2) + (i*yAdj)];
                        SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                        SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) == 0
                            %                     subplot(1,NumDynTrials, Selection(q));
                            %                     line(a,b,'Color','k');
                            DynamicPPTrials(Selection(q)).L.LatHeelCent{j}(2,:) = [a(1), b(1)];
                            DynamicPPTrials(Selection(q)).L.MedHeelCent{j}(2,:) = [a(2), b(2)];
                            break
                        end
                    end
                else % if the FP points are outside the border of the foot, search from the toe back towards the heel
                    for i = Loop % search forward to make sure there arent more pressures forward
                        a = [DynamicPPTrials(Selection(q)).L.LatHeelCent{j}(2,1) + (i*xAdj), DynamicPPTrials(Selection(q)).L.MedHeelCent{j}(2,1) + (i*xAdj)];
                        b = [DynamicPPTrials(Selection(q)).L.LatHeelCent{j}(2,2) + (i*yAdj), DynamicPPTrials(Selection(q)).L.MedHeelCent{j}(2,2) + (i*yAdj)];
                        SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                        SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) > 0
                            PressForward = 'Yes';
                            ForwardRegion = [a b];
                            break
                        else
                            PressForward = 'No';
                        end
                    end
                    if strcmp(PressForward,'No') == 1
                        for i = Loop % increment search points back from toe
                            a = [DynamicPPTrials(Selection(q)).L.LatHeelCent{j}(2,1) - (i*xAdj), DynamicPPTrials(Selection(q)).L.MedHeelCent{j}(2,1) - (i*xAdj)];
                            b = [DynamicPPTrials(Selection(q)).L.LatHeelCent{j}(2,2) - (i*yAdj), DynamicPPTrials(Selection(q)).L.MedHeelCent{j}(2,2) - (i*yAdj)];
                            SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                            SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                            Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                            Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                            if sum(sum(Region)) > 0
                                %                         subplot(1,NumDynTrials, Selection(q));
                                %                         line(a,b,'Color','k');
                                DynamicPPTrials(Selection(q)).L.LatHeelCent{j}(2,:) = [a(1), b(1)];
                                DynamicPPTrials(Selection(q)).L.MedHeelCent{j}(2,:) = [a(2), b(2)];
                                break
                            end
                        end
                    else % pressures are in front of the original toe point, keep searching forward
                        for i = Loop
                            a = [ForwardRegion(1) + (i*xAdj), ForwardRegion(2) + (i*xAdj)];
                            b = [ForwardRegion(3) + (i*yAdj), ForwardRegion(4) + (i*yAdj)];
                            SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                            SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                            Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                            Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                            if sum(sum(Region)) > 0
                                %                         subplot(1,NumDynTrials, Selection(q));
                                %                         line(a,b,'Color','k');
                                DynamicPPTrials(Selection(q)).L.LatHeelCent{j}(2,:) = [a(1), b(1)];
                                DynamicPPTrials(Selection(q)).L.MedHeelCent{j}(2,:) = [a(2), b(2)];
                                break
                            end
                        end
                    end
                end
                % determine if the current HEEL point is on a pressure
                a = [DynamicPPTrials(Selection(q)).L.LatHeelCent{j}(1,1), DynamicPPTrials(Selection(q)).L.MedHeelCent{j}(1,1)];
                b = [DynamicPPTrials(Selection(q)).L.LatHeelCent{j}(1,2), DynamicPPTrials(Selection(q)).L.MedHeelCent{j}(1,2)];
                SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                if sum(sum(Region)) > 0
                    PtOnOff = 'On';
                else
                    PtOnOff = 'Off';
                end
                if strcmp(PtOnOff,'On') == 1 % if the FP points are within the border of the foot, search from the Heel backward
                    for i = Loop % increment search points out from heel
                        a = [DynamicPPTrials(Selection(q)).L.LatHeelCent{j}(1,1) - (i*xAdj), DynamicPPTrials(Selection(q)).L.MedHeelCent{j}(1,1) - (i*xAdj)];
                        b = [DynamicPPTrials(Selection(q)).L.LatHeelCent{j}(1,2) - (i*yAdj), DynamicPPTrials(Selection(q)).L.MedHeelCent{j}(1,2) - (i*yAdj)];
                        SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                        SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) == 0
                            %                     subplot(1,NumDynTrials, Selection(q));
                            %                     line(a,b,'Color','k');
                            DynamicPPTrials(Selection(q)).L.LatHeelCent{j}(1,:) = [a(1), b(1)];
                            DynamicPPTrials(Selection(q)).L.MedHeelCent{j}(1,:) = [a(2), b(2)];
                            break
                        end
                    end
                else % if the FP points are outside the border of the foot, search from the heel towards the toe
                    for i = Loop % increment search points forward from heel
                        a = [DynamicPPTrials(Selection(q)).L.LatHeelCent{j}(1,1) + (i*xAdj), DynamicPPTrials(Selection(q)).L.MedHeelCent{j}(1,1) + (i*xAdj)];
                        b = [DynamicPPTrials(Selection(q)).L.LatHeelCent{j}(1,2) + (i*yAdj), DynamicPPTrials(Selection(q)).L.MedHeelCent{j}(1,2) + (i*yAdj)];
                        SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                        SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) > 0
                            %                     subplot(1,NumDynTrials, Selection(q));
                            %                     line(a,b,'Color','k');
                            DynamicPPTrials(Selection(q)).L.LatHeelCent{j}(1,:) = [a(1), b(1)];
                            DynamicPPTrials(Selection(q)).L.MedHeelCent{j}(1,:) = [a(2), b(2)];
                            break
                        end
                    end
                end
            end
        end
    end
end

% RIGHT foot length
for q = 1:length(Selection)
    for j = 1: DynamicPPTrials(Selection(q)).NumRight
        % determine if the current TOE point is on a pressure
        if DynamicPPTrials(Selection(q)).PoorPressure.Right(j) == 0
            %% manual points
            if isfield(DynamicPPTrials(q).Rht, 'Manual') == 1
                if strcmp(PPSettings.MaskChoice, 'Manual') || strcmp(PPSettings.MaskChoice, 'Manual (pink)') || strcmp(PPSettings.MaskChoice, 'Validation')
                    xAdj = DynamicPPTrials(Selection(q)).Rht.ManualV{j}(:,1);
                    yAdj = DynamicPPTrials(Selection(q)).Rht.ManualV{j}(:,2);
                    a = [DynamicPPTrials(Selection(q)).R.LatManual{j}(2,1), DynamicPPTrials(Selection(q)).R.MedManual{j}(2,1)];
                    b = [DynamicPPTrials(Selection(q)).R.LatManual{j}(2,2), DynamicPPTrials(Selection(q)).R.MedManual{j}(2,2)];
                    SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                    SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                    %             subplot(1,NumDynTrials, Selection(q)); hold on;  plot(SrchAreaX, SrchAreaY, '.g');
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) > 0
                        PtOnOff = 'On';
                    else
                        PtOnOff = 'Off';
                    end
                    if strcmp(PtOnOff,'On') == 1 % if the FP points are within the border of the foot, search from the toe forward
                        for i = Loop % increment search points out from toe
                            a = [DynamicPPTrials(Selection(q)).R.LatManual{j}(2,1) + (i*xAdj), DynamicPPTrials(Selection(q)).R.MedManual{j}(2,1) + (i*xAdj)];
                            b = [DynamicPPTrials(Selection(q)).R.LatManual{j}(2,2) + (i*yAdj), DynamicPPTrials(Selection(q)).R.MedManual{j}(2,2) + (i*yAdj)];
                            SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                            SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                            Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                            Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                            if sum(sum(Region)) == 0
                                %                         subplot(1,NumDynTrials, Selection(q));
                                %                         line(a,b,'Color','k');
                                DynamicPPTrials(Selection(q)).R.LatManual{j}(2,:) = [a(1), b(1)];
                                DynamicPPTrials(Selection(q)).R.MedManual{j}(2,:) = [a(2), b(2)];
                                % uncomment to check search points
                                %                 subplot(1, NumDynTrials, q);  hold on;
                                %                 for ii = 1:4
                                %                     plot(SrchAreaX(ii), SrchAreaY(ii),'k.')
                                %                 end
                                %                 plot(a(:), b(:), 'r')
                                break
                            end
                        end
                    else % if the FP points are outside the border of the foot, search from the toe back towards the heel
                        for i = Loop % search forward to make sure there arent more pressures forward
                            a = [DynamicPPTrials(Selection(q)).R.LatManual{j}(2,1) + (i*xAdj), DynamicPPTrials(Selection(q)).R.MedManual{j}(2,1) + (i*xAdj)];
                            b = [DynamicPPTrials(Selection(q)).R.LatManual{j}(2,2) + (i*yAdj), DynamicPPTrials(Selection(q)).R.MedManual{j}(2,2) + (i*yAdj)];
                            SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                            SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                            Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                            Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                            if sum(sum(Region)) > 0
                                PressForward = 'Yes';
                                ForwardRegion = [a b];
                                break
                            else
                                PressForward = 'No';
                            end
                        end
                        if strcmp(PressForward,'No') == 1
                            for i = Loop % increment search points back from toe
                                a = [DynamicPPTrials(Selection(q)).R.LatManual{j}(2,1) - (i*xAdj), DynamicPPTrials(Selection(q)).R.MedManual{j}(2,1) - (i*xAdj)];
                                b = [DynamicPPTrials(Selection(q)).R.LatManual{j}(2,2) - (i*yAdj), DynamicPPTrials(Selection(q)).R.MedManual{j}(2,2) - (i*yAdj)];
                                SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                                SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                                Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                                Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                                if sum(sum(Region)) > 0
                                    %                             subplot(1,NumDynTrials, Selection(q));
                                    %                             line(a,b,'Color','k');
                                    DynamicPPTrials(Selection(q)).R.LatManual{j}(2,:) = [a(1), b(1)];
                                    DynamicPPTrials(Selection(q)).R.MedManual{j}(2,:) = [a(2), b(2)];
                                    % uncomment to check search points
                                    %                                      subplot(1, NumDynTrials, q);  hold on;
                                    %                                      for ii = 1:4
                                    %                                          plot(SrchAreaX(ii), SrchAreaY(ii),'k.')
                                    %                                      end
                                    %                                      plot(a(:), b(:), 'r')
                                    break
                                end
                            end
                        else
                            for i = Loop % increment search points out from toe
                                a = [ForwardRegion(1) + (i*xAdj), ForwardRegion(2) + (i*xAdj)];
                                b = [ForwardRegion(3) + (i*yAdj), ForwardRegion(4) + (i*yAdj)];
                                SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                                SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                                Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                                Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                                if sum(sum(Region)) == 0
                                    %                             subplot(1,NumDynTrials, Selection(q));
                                    %                             line(a,b,'Color','k');
                                    DynamicPPTrials(Selection(q)).R.LatManual{j}(2,:) = [a(1), b(1)];
                                    DynamicPPTrials(Selection(q)).R.MedManual{j}(2,:) = [a(2), b(2)];
                                    % uncomment to check search points
                                    %                 subplot(1, NumDynTrials, q); hold on;
                                    %                 for ii = 1:4
                                    %                     plot(SrchAreaX(ii), SrchAreaY(ii),'k.')
                                    %                 end
                                    %                 plot(a(:), b(:), 'r')
                                    break
                                end
                            end
                        end
                    end
                    
                    % determine if the current HEEL point is on a pressure
                    a = [DynamicPPTrials(Selection(q)).R.LatManual{j}(1,1), DynamicPPTrials(Selection(q)).R.MedManual{j}(1,1)];
                    b = [DynamicPPTrials(Selection(q)).R.LatManual{j}(1,2), DynamicPPTrials(Selection(q)).R.MedManual{j}(1,2)];
                    SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                    SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                    %             subplot(1,NumDynTrials, Selection(q)); hold on;  plot(SrchAreaX, SrchAreaY, '.g');
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) > 0
                        PtOnOff = 'On';
                    else
                        PtOnOff = 'Off';
                    end
                    if strcmp(PtOnOff,'On') == 1 % if the FP points are within the border of the foot, search from the Heel backward
                        for i = Loop % increment search points out from heel
                            a = [DynamicPPTrials(Selection(q)).R.LatManual{j}(1,1) - (i*xAdj), DynamicPPTrials(Selection(q)).R.MedManual{j}(1,1) - (i*xAdj)];
                            b = [DynamicPPTrials(Selection(q)).R.LatManual{j}(1,2) - (i*yAdj), DynamicPPTrials(Selection(q)).R.MedManual{j}(1,2) - (i*yAdj)];
                            SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                            SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                            Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                            Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                            if sum(sum(Region)) == 0
                                %                         subplot(1,NumDynTrials, Selection(q));
                                %                         line(a,b,'Color','k');
                                DynamicPPTrials(Selection(q)).R.LatManual{j}(1,:) = [a(1), b(1)];
                                DynamicPPTrials(Selection(q)).R.MedManual{j}(1,:) = [a(2), b(2)];
                                % uncomment to check search points
                                %                 subplot(1, NumDynTrials, q); hold on;
                                %                 for ii = 1:4
                                %                     plot(SrchAreaX(ii), SrchAreaY(ii),'k.')
                                %                 end
                                %                 plot(a(:), b(:), 'r')
                                break
                            end
                        end
                    else % if the FP points are outside the border of the foot, search from the heel towards the toe
                        for i = Loop % increment search points forward from heel
                            a = [DynamicPPTrials(Selection(q)).R.LatManual{j}(1,1) + (i*xAdj), DynamicPPTrials(Selection(q)).R.MedManual{j}(1,1) + (i*xAdj)];
                            b = [DynamicPPTrials(Selection(q)).R.LatManual{j}(1,2) + (i*yAdj), DynamicPPTrials(Selection(q)).R.MedManual{j}(1,2) + (i*yAdj)];
                            SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                            SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                            Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                            Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                            if sum(sum(Region)) > 0
                                %                         subplot(1,NumDynTrials, Selection(q));
                                %                         line(a,b,'Color','k');
                                DynamicPPTrials(Selection(q)).R.LatManual{j}(1,:) = [a(1), b(1)];
                                DynamicPPTrials(Selection(q)).R.MedManual{j}(1,:) = [a(2), b(2)];
                                % uncomment to check search points
                                %                                      subplot(1, NumDynTrials, q); hold on;
                                %                                      for ii = 1:4
                                %                                          plot(SrchAreaX(ii), SrchAreaY(ii),'k.')
                                %                                      end
                                %                                      plot(a(:), b(:), 'r')
                                break
                            end
                        end
                    end
                end
                
                %% General image processing points
                xAdj = DynamicPPTrials(Selection(q)).Rht.GeneralV{j}(:,1);
                yAdj = DynamicPPTrials(Selection(q)).Rht.GeneralV{j}(:,2);
                a = [DynamicPPTrials(Selection(q)).R.LatGeneral{j}(2,1), DynamicPPTrials(Selection(q)).R.MedGeneral{j}(2,1)];
                b = [DynamicPPTrials(Selection(q)).R.LatGeneral{j}(2,2), DynamicPPTrials(Selection(q)).R.MedGeneral{j}(2,2)];
                SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                %             subplot(1,NumDynTrials, Selection(q)); hold on;  plot(SrchAreaX, SrchAreaY, '.g');
                Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                if sum(sum(Region)) > 0
                    PtOnOff = 'On';
                else
                    PtOnOff = 'Off';
                end
                if strcmp(PtOnOff,'On') == 1 % if the FP points are within the border of the foot, search from the toe forward
                    for i = Loop % increment search points out from toe
                        a = [DynamicPPTrials(Selection(q)).R.LatGeneral{j}(2,1) + (i*xAdj), DynamicPPTrials(Selection(q)).R.MedGeneral{j}(2,1) + (i*xAdj)];
                        b = [DynamicPPTrials(Selection(q)).R.LatGeneral{j}(2,2) + (i*yAdj), DynamicPPTrials(Selection(q)).R.MedGeneral{j}(2,2) + (i*yAdj)];
                        SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                        SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) == 0
                            %                         subplot(1,NumDynTrials, Selection(q));
                            %                         line(a,b,'Color','k');
                            DynamicPPTrials(Selection(q)).R.LatGeneral{j}(2,:) = [a(1), b(1)];
                            DynamicPPTrials(Selection(q)).R.MedGeneral{j}(2,:) = [a(2), b(2)];
                            % uncomment to check search points
                            %                 subplot(1, NumDynTrials, q);  hold on;
                            %                 for ii = 1:4
                            %                     plot(SrchAreaX(ii), SrchAreaY(ii),'k.')
                            %                 end
                            %                 plot(a(:), b(:), 'r')
                            break
                        end
                    end
                else % if the FP points are outside the border of the foot, search from the toe back towards the heel
                    for i = Loop % search forward to make sure there arent more pressures forward
                        a = [DynamicPPTrials(Selection(q)).R.LatGeneral{j}(2,1) + (i*xAdj), DynamicPPTrials(Selection(q)).R.MedGeneral{j}(2,1) + (i*xAdj)];
                        b = [DynamicPPTrials(Selection(q)).R.LatGeneral{j}(2,2) + (i*yAdj), DynamicPPTrials(Selection(q)).R.MedGeneral{j}(2,2) + (i*yAdj)];
                        SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                        SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) > 0
                            PressForward = 'Yes';
                            ForwardRegion = [a b];
                            break
                        else
                            PressForward = 'No';
                        end
                    end
                    if strcmp(PressForward,'No') == 1
                        for i = Loop % increment search points back from toe
                            a = [DynamicPPTrials(Selection(q)).R.LatGeneral{j}(2,1) - (i*xAdj), DynamicPPTrials(Selection(q)).R.MedGeneral{j}(2,1) - (i*xAdj)];
                            b = [DynamicPPTrials(Selection(q)).R.LatGeneral{j}(2,2) - (i*yAdj), DynamicPPTrials(Selection(q)).R.MedGeneral{j}(2,2) - (i*yAdj)];
                            SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                            SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                            Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                            Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                            if sum(sum(Region)) > 0
                                %                             subplot(1,NumDynTrials, Selection(q));
                                %                             line(a,b,'Color','k');
                                DynamicPPTrials(Selection(q)).R.LatGeneral{j}(2,:) = [a(1), b(1)];
                                DynamicPPTrials(Selection(q)).R.MedGeneral{j}(2,:) = [a(2), b(2)];
                                % uncomment to check search points
                                %                                      subplot(1, NumDynTrials, q);  hold on;
                                %                                      for ii = 1:4
                                %                                          plot(SrchAreaX(ii), SrchAreaY(ii),'k.')
                                %                                      end
                                %                                      plot(a(:), b(:), 'r')
                                break
                            end
                        end
                    else
                        for i = Loop % increment search points out from toe
                            a = [ForwardRegion(1) + (i*xAdj), ForwardRegion(2) + (i*xAdj)];
                            b = [ForwardRegion(3) + (i*yAdj), ForwardRegion(4) + (i*yAdj)];
                            SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                            SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                            Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                            Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                            if sum(sum(Region)) == 0
                                %                             subplot(1,NumDynTrials, Selection(q));
                                %                             line(a,b,'Color','k');
                                DynamicPPTrials(Selection(q)).R.LatGeneral{j}(2,:) = [a(1), b(1)];
                                DynamicPPTrials(Selection(q)).R.MedGeneral{j}(2,:) = [a(2), b(2)];
                                % uncomment to check search points
                                %                 subplot(1, NumDynTrials, q); hold on;
                                %                 for ii = 1:4
                                %                     plot(SrchAreaX(ii), SrchAreaY(ii),'k.')
                                %                 end
                                %                 plot(a(:), b(:), 'r')
                                break
                            end
                        end
                    end
                end
                
                % determine if the current HEEL point is on a pressure
                a = [DynamicPPTrials(Selection(q)).R.LatGeneral{j}(1,1), DynamicPPTrials(Selection(q)).R.MedGeneral{j}(1,1)];
                b = [DynamicPPTrials(Selection(q)).R.LatGeneral{j}(1,2), DynamicPPTrials(Selection(q)).R.MedGeneral{j}(1,2)];
                SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                %             subplot(1,NumDynTrials, Selection(q)); hold on;  plot(SrchAreaX, SrchAreaY, '.g');
                Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                if sum(sum(Region)) > 0
                    PtOnOff = 'On';
                else
                    PtOnOff = 'Off';
                end
                if strcmp(PtOnOff,'On') == 1 % if the FP points are within the border of the foot, search from the Heel backward
                    for i = Loop % increment search points out from heel
                        a = [DynamicPPTrials(Selection(q)).R.LatGeneral{j}(1,1) - (i*xAdj), DynamicPPTrials(Selection(q)).R.MedGeneral{j}(1,1) - (i*xAdj)];
                        b = [DynamicPPTrials(Selection(q)).R.LatGeneral{j}(1,2) - (i*yAdj), DynamicPPTrials(Selection(q)).R.MedGeneral{j}(1,2) - (i*yAdj)];
                        SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                        SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) == 0
                            %                         subplot(1,NumDynTrials, Selection(q));
                            %                         line(a,b,'Color','k');
                            DynamicPPTrials(Selection(q)).R.LatGeneral{j}(1,:) = [a(1), b(1)];
                            DynamicPPTrials(Selection(q)).R.MedGeneral{j}(1,:) = [a(2), b(2)];
                            % uncomment to check search points
                            %                 subplot(1, NumDynTrials, q); hold on;
                            %                 for ii = 1:4
                            %                     plot(SrchAreaX(ii), SrchAreaY(ii),'k.')
                            %                 end
                            %                 plot(a(:), b(:), 'r')
                            break
                        end
                    end
                else % if the FP points are outside the border of the foot, search from the heel towards the toe
                    for i = Loop % increment search points forward from heel
                        a = [DynamicPPTrials(Selection(q)).R.LatGeneral{j}(1,1) + (i*xAdj), DynamicPPTrials(Selection(q)).R.MedGeneral{j}(1,1) + (i*xAdj)];
                        b = [DynamicPPTrials(Selection(q)).R.LatGeneral{j}(1,2) + (i*yAdj), DynamicPPTrials(Selection(q)).R.MedGeneral{j}(1,2) + (i*yAdj)];
                        SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                        SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) > 0
                            %                         subplot(1,NumDynTrials, Selection(q));
                            %                         line(a,b,'Color','k');
                            DynamicPPTrials(Selection(q)).R.LatGeneral{j}(1,:) = [a(1), b(1)];
                            DynamicPPTrials(Selection(q)).R.MedGeneral{j}(1,:) = [a(2), b(2)];
                            % uncomment to check search points
                            %                                      subplot(1, NumDynTrials, q); hold on;
                            %                                      for ii = 1:4
                            %                                          plot(SrchAreaX(ii), SrchAreaY(ii),'k.')
                            %                                      end
                            %                                      plot(a(:), b(:), 'r')
                            break
                        end
                    end
                end
            end
            
            %% CoP points
            % determine if the current TOE point is on a pressure
            xAdj = DynamicPPTrials(Selection(q)).Rht.VCoP{j}(:,1);
            yAdj = DynamicPPTrials(Selection(q)).Rht.VCoP{j}(:,2);
            a = [DynamicPPTrials(Selection(q)).R.LatCoP{j}(2,1), DynamicPPTrials(Selection(q)).R.MedCoP{j}(2,1)];
            b = [DynamicPPTrials(Selection(q)).R.LatCoP{j}(2,2), DynamicPPTrials(Selection(q)).R.MedCoP{j}(2,2)];
            SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
            SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
            Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
            Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
            if sum(sum(Region)) > 0
                PtOnOff = 'On';
            else
                PtOnOff = 'Off';
            end
            if strcmp(PtOnOff,'On') == 1 % if the FP points are within the border of the foot, search from the toe forward
                for i = Loop % increment search points out from toe
                    a = [DynamicPPTrials(Selection(q)).R.LatCoP{j}(2,1) + (i*xAdj), DynamicPPTrials(Selection(q)).R.MedCoP{j}(2,1) + (i*xAdj)];
                    b = [DynamicPPTrials(Selection(q)).R.LatCoP{j}(2,2) + (i*yAdj), DynamicPPTrials(Selection(q)).R.MedCoP{j}(2,2) + (i*yAdj)];
                    SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                    SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) == 0
                        %                     subplot(1,NumDynTrials, Selection(q));
                        %                     line(a,b,'Color','k');
                        DynamicPPTrials(Selection(q)).R.LatCoP{j}(2,:) = [a(1), b(1)];
                        DynamicPPTrials(Selection(q)).R.MedCoP{j}(2,:) = [a(2), b(2)];
                        break
                    end
                end
            else % if the FP points are outside the border of the foot, search from the toe back towards the heel
                for i = Loop % search forward to make sure there arent more pressures forward
                    a = [DynamicPPTrials(Selection(q)).R.LatCoP{j}(2,1) + (i*xAdj), DynamicPPTrials(Selection(q)).R.MedCoP{j}(2,1) + (i*xAdj)];
                    b = [DynamicPPTrials(Selection(q)).R.LatCoP{j}(2,2) + (i*yAdj), DynamicPPTrials(Selection(q)).R.MedCoP{j}(2,2) + (i*yAdj)];
                    SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                    SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) > 0
                        PressForward = 'Yes';
                        ForwardRegion = [a b];
                        break
                    else
                        PressForward = 'No';
                    end
                end
                if strcmp(PressForward,'No') == 1
                    for i = Loop % increment search points back from toe
                        a = [DynamicPPTrials(Selection(q)).R.LatCoP{j}(2,1) - (i*xAdj), DynamicPPTrials(Selection(q)).R.MedCoP{j}(2,1) - (i*xAdj)];
                        b = [DynamicPPTrials(Selection(q)).R.LatCoP{j}(2,2) - (i*yAdj), DynamicPPTrials(Selection(q)).R.MedCoP{j}(2,2) - (i*yAdj)];
                        SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                        SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) > 0
                            %                         subplot(1,NumDynTrials, Selection(q));
                            %                         line(a,b,'Color','k');
                            DynamicPPTrials(Selection(q)).R.LatCoP{j}(2,:) = [a(1), b(1)];
                            DynamicPPTrials(Selection(q)).R.MedCoP{j}(2,:) = [a(2), b(2)];
                            break
                        end
                    end
                else % pressures are in front of the original toe point, keep searching forward
                    for i = Loop
                        a = [ForwardRegion(1) + (i*xAdj), ForwardRegion(2) + (i*xAdj)];
                        b = [ForwardRegion(3) + (i*yAdj), ForwardRegion(4) + (i*yAdj)];
                        SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                        SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) > 0
                            %                         subplot(1,NumDynTrials, Selection(q));
                            %                         line(a,b,'Color','k');
                            DynamicPPTrials(Selection(q)).R.LatCoP{j}(2,:) = [a(1), b(1)];
                            DynamicPPTrials(Selection(q)).R.MedCoP{j}(2,:) = [a(2), b(2)];
                            break
                        end
                    end
                end
            end
            % determine if the current HEEL point is on a pressure
            a = [DynamicPPTrials(Selection(q)).R.LatCoP{j}(1,1), DynamicPPTrials(Selection(q)).R.MedCoP{j}(1,1)];
            b = [DynamicPPTrials(Selection(q)).R.LatCoP{j}(1,2), DynamicPPTrials(Selection(q)).R.MedCoP{j}(1,2)];
            SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
            SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
            Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
            Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
            if sum(sum(Region)) > 0
                PtOnOff = 'On';
            else
                PtOnOff = 'Off';
            end
            if strcmp(PtOnOff,'On') == 1 % if the FP points are within the border of the foot, search from the Heel backward
                for i = Loop % increment search points out from heel
                    a = [DynamicPPTrials(Selection(q)).R.LatCoP{j}(1,1) - (i*xAdj), DynamicPPTrials(Selection(q)).R.MedCoP{j}(1,1) - (i*xAdj)];
                    b = [DynamicPPTrials(Selection(q)).R.LatCoP{j}(1,2) - (i*yAdj), DynamicPPTrials(Selection(q)).R.MedCoP{j}(1,2) - (i*yAdj)];
                    SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                    SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) == 0
                        %                     subplot(1,NumDynTrials, Selection(q));
                        %                     line(a,b,'Color','k');
                        DynamicPPTrials(Selection(q)).R.LatCoP{j}(1,:) = [a(1), b(1)];
                        DynamicPPTrials(Selection(q)).R.MedCoP{j}(1,:) = [a(2), b(2)];
                        break
                    end
                end
            else % if the FP points are outside the border of the foot, search from the heel towards the toe
                for i = Loop % increment search points forward from heel
                    a = [DynamicPPTrials(Selection(q)).R.LatCoP{j}(1,1) + (i*xAdj), DynamicPPTrials(Selection(q)).R.MedCoP{j}(1,1) + (i*xAdj)];
                    b = [DynamicPPTrials(Selection(q)).R.LatCoP{j}(1,2) + (i*yAdj), DynamicPPTrials(Selection(q)).R.MedCoP{j}(1,2) + (i*yAdj)];
                    SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                    SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) > 0
                        %                     subplot(1,NumDynTrials, Selection(q));
                        %                     line(a,b,'Color','k');
                        DynamicPPTrials(Selection(q)).R.LatCoP{j}(1,:) = [a(1), b(1)];
                        DynamicPPTrials(Selection(q)).R.MedCoP{j}(1,:) = [a(2), b(2)];
                        break
                    end
                end
            end
            
            %% 66% CoP points
            % determine if the current TOE point is on a pressure
            xAdj = DynamicPPTrials(Selection(q)).Rht.CoPV66{j}(:,1);
            yAdj = DynamicPPTrials(Selection(q)).Rht.CoPV66{j}(:,2);
            a = [DynamicPPTrials(Selection(q)).R.Lat66{j}(2,1), DynamicPPTrials(Selection(q)).R.Med66{j}(2,1)];
            b = [DynamicPPTrials(Selection(q)).R.Lat66{j}(2,2), DynamicPPTrials(Selection(q)).R.Med66{j}(2,2)];
            SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
            SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
            Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
            Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
            if sum(sum(Region)) > 0
                PtOnOff = 'On';
            else
                PtOnOff = 'Off';
            end
            if strcmp(PtOnOff,'On') == 1 % if the FP points are within the border of the foot, search from the toe forward
                for i = Loop % increment search points out from toe
                    a = [DynamicPPTrials(Selection(q)).R.Lat66{j}(2,1) + (i*xAdj), DynamicPPTrials(Selection(q)).R.Med66{j}(2,1) + (i*xAdj)];
                    b = [DynamicPPTrials(Selection(q)).R.Lat66{j}(2,2) + (i*yAdj), DynamicPPTrials(Selection(q)).R.Med66{j}(2,2) + (i*yAdj)];
                    SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                    SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) == 0
                        %                     subplot(1,NumDynTrials, Selection(q));
                        %                     line(a,b,'Color','k');
                        DynamicPPTrials(Selection(q)).R.Lat66{j}(2,:) = [a(1), b(1)];
                        DynamicPPTrials(Selection(q)).R.Med66{j}(2,:) = [a(2), b(2)];
                        break
                    end
                end
            else % if the FP points are outside the border of the foot, search from the toe back towards the heel
                for i = Loop % search forward to make sure there arent more pressures forward
                    a = [DynamicPPTrials(Selection(q)).R.Lat66{j}(2,1) + (i*xAdj), DynamicPPTrials(Selection(q)).R.Med66{j}(2,1) + (i*xAdj)];
                    b = [DynamicPPTrials(Selection(q)).R.Lat66{j}(2,2) + (i*yAdj), DynamicPPTrials(Selection(q)).R.Med66{j}(2,2) + (i*yAdj)];
                    SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                    SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) > 0
                        PressForward = 'Yes';
                        ForwardRegion = [a b];
                        break
                    else
                        PressForward = 'No';
                    end
                end
                if strcmp(PressForward,'No') == 1
                    for i = Loop % increment search points back from toe
                        a = [DynamicPPTrials(Selection(q)).R.Lat66{j}(2,1) - (i*xAdj), DynamicPPTrials(Selection(q)).R.Med66{j}(2,1) - (i*xAdj)];
                        b = [DynamicPPTrials(Selection(q)).R.Lat66{j}(2,2) - (i*yAdj), DynamicPPTrials(Selection(q)).R.Med66{j}(2,2) - (i*yAdj)];
                        SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                        SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) > 0
                            %                         subplot(1,NumDynTrials, Selection(q));
                            %                         line(a,b,'Color','k');
                            DynamicPPTrials(Selection(q)).R.Lat66{j}(2,:) = [a(1), b(1)];
                            DynamicPPTrials(Selection(q)).R.Med66{j}(2,:) = [a(2), b(2)];
                            break
                        end
                    end
                else % pressures are in front of the original toe point, keep searching forward
                    for i = Loop
                        a = [ForwardRegion(1) + (i*xAdj), ForwardRegion(2) + (i*xAdj)];
                        b = [ForwardRegion(3) + (i*yAdj), ForwardRegion(4) + (i*yAdj)];
                        SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                        SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) > 0
                            %                         subplot(1,NumDynTrials, Selection(q));
                            %                         line(a,b,'Color','k');
                            DynamicPPTrials(Selection(q)).R.Lat66{j}(2,:) = [a(1), b(1)];
                            DynamicPPTrials(Selection(q)).R.Med66{j}(2,:) = [a(2), b(2)];
                            break
                        end
                    end
                end
            end
            % determine if the current HEEL point is on a pressure
            a = [DynamicPPTrials(Selection(q)).R.Lat66{j}(1,1), DynamicPPTrials(Selection(q)).R.Med66{j}(1,1)];
            b = [DynamicPPTrials(Selection(q)).R.Lat66{j}(1,2), DynamicPPTrials(Selection(q)).R.Med66{j}(1,2)];
            SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
            SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
            Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
            Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
            if sum(sum(Region)) > 0
                PtOnOff = 'On';
            else
                PtOnOff = 'Off';
            end
            if strcmp(PtOnOff,'On') == 1 % if the FP points are within the border of the foot, search from the Heel backward
                for i = Loop % increment search points out from heel
                    a = [DynamicPPTrials(Selection(q)).R.Lat66{j}(1,1) - (i*xAdj), DynamicPPTrials(Selection(q)).R.Med66{j}(1,1) - (i*xAdj)];
                    b = [DynamicPPTrials(Selection(q)).R.Lat66{j}(1,2) - (i*yAdj), DynamicPPTrials(Selection(q)).R.Med66{j}(1,2) - (i*yAdj)];
                    SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                    SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) == 0
                        %                     subplot(1,NumDynTrials, Selection(q));
                        %                     line(a,b,'Color','k');
                        DynamicPPTrials(Selection(q)).R.Lat66{j}(1,:) = [a(1), b(1)];
                        DynamicPPTrials(Selection(q)).R.Med66{j}(1,:) = [a(2), b(2)];
                        break
                    end
                end
            else % if the FP points are outside the border of the foot, search from the heel towards the toe
                for i = Loop % increment search points forward from heel
                    a = [DynamicPPTrials(Selection(q)).R.Lat66{j}(1,1) + (i*xAdj), DynamicPPTrials(Selection(q)).R.Med66{j}(1,1) + (i*xAdj)];
                    b = [DynamicPPTrials(Selection(q)).R.Lat66{j}(1,2) + (i*yAdj), DynamicPPTrials(Selection(q)).R.Med66{j}(1,2) + (i*yAdj)];
                    SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                    SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                    Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                    Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                    if sum(sum(Region)) > 0
                        %                     subplot(1,NumDynTrials, Selection(q));
                        %                     line(a,b,'Color','k');
                        DynamicPPTrials(Selection(q)).R.Lat66{j}(1,:) = [a(1), b(1)];
                        DynamicPPTrials(Selection(q)).R.Med66{j}(1,:) = [a(2), b(2)];
                        break
                    end
                end
            end
            
            %% IP points
            if Block.IP == 0
                % determine if the current TOE point is on a pressure
                xAdj = DynamicPPTrials(Selection(q)).Rht.VIP{j}(:,1);
                yAdj = DynamicPPTrials(Selection(q)).Rht.VIP{j}(:,2);
                if isempty(DynamicPPTrials(Selection(q)).R.LatIP)
                    break
                end
                a = [DynamicPPTrials(Selection(q)).R.LatIP{j}(2,1), DynamicPPTrials(Selection(q)).R.MedIP{j}(2,1)];
                b = [DynamicPPTrials(Selection(q)).R.LatIP{j}(2,2), DynamicPPTrials(Selection(q)).R.MedIP{j}(2,2)];
                SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                if sum(sum(Region)) > 0
                    PtOnOff = 'On';
                else
                    PtOnOff = 'Off';
                end
                if strcmp(PtOnOff,'On') == 1 % if the FP points are within the border of the foot, search from the toe forward
                    for i = Loop % increment search points out from toe
                        a = [DynamicPPTrials(Selection(q)).R.LatIP{j}(2,1) + (i*xAdj), DynamicPPTrials(Selection(q)).R.MedIP{j}(2,1) + (i*xAdj)];
                        b = [DynamicPPTrials(Selection(q)).R.LatIP{j}(2,2) + (i*yAdj), DynamicPPTrials(Selection(q)).R.MedIP{j}(2,2) + (i*yAdj)];
                        SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                        SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) == 0
                            %                     subplot(1,NumDynTrials, Selection(q));
                            %                     line(a,b,'Color','k');
                            DynamicPPTrials(Selection(q)).R.LatIP{j}(2,:) = [a(1), b(1)];
                            DynamicPPTrials(Selection(q)).R.MedIP{j}(2,:) = [a(2), b(2)];
                            break
                        end
                    end
                else % if the FP points are outside the border of the foot, search from the toe back towards the heel
                    for i = Loop % search forward to make sure there arent more pressures forward
                        a = [DynamicPPTrials(Selection(q)).R.LatIP{j}(2,1) + (i*xAdj), DynamicPPTrials(Selection(q)).R.MedIP{j}(2,1) + (i*xAdj)];
                        b = [DynamicPPTrials(Selection(q)).R.LatIP{j}(2,2) + (i*yAdj), DynamicPPTrials(Selection(q)).R.MedIP{j}(2,2) + (i*yAdj)];
                        SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                        SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) > 0
                            PressForward = 'Yes';
                            ForwardRegion = [a b];
                            break
                        else
                            PressForward = 'No';
                        end
                    end
                    if strcmp(PressForward,'No') == 1
                        for i = Loop % increment search points back from toe
                            a = [DynamicPPTrials(Selection(q)).R.LatIP{j}(2,1) - (i*xAdj), DynamicPPTrials(Selection(q)).R.MedIP{j}(2,1) - (i*xAdj)];
                            b = [DynamicPPTrials(Selection(q)).R.LatIP{j}(2,2) - (i*yAdj), DynamicPPTrials(Selection(q)).R.MedIP{j}(2,2) - (i*yAdj)];
                            SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                            SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                            Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                            Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                            if sum(sum(Region)) > 0
                                %                         subplot(1,NumDynTrials, Selection(q));
                                %                         line(a,b,'Color','k');
                                DynamicPPTrials(Selection(q)).R.LatIP{j}(2,:) = [a(1), b(1)];
                                DynamicPPTrials(Selection(q)).R.MedIP{j}(2,:) = [a(2), b(2)];
                                break
                            end
                        end
                    else % pressures are in front of the original toe point, keep searching forward
                        for i = Loop
                            a = [ForwardRegion(1) + (i*xAdj), ForwardRegion(2) + (i*xAdj)];
                            b = [ForwardRegion(3) + (i*yAdj), ForwardRegion(4) + (i*yAdj)];
                            SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                            SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                            Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                            Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                            if sum(sum(Region)) > 0
                                %                         subplot(1,NumDynTrials, Selection(q));
                                %                         line(a,b,'Color','k');
                                DynamicPPTrials(Selection(q)).R.LatIP{j}(2,:) = [a(1), b(1)];
                                DynamicPPTrials(Selection(q)).R.MedIP{j}(2,:) = [a(2), b(2)];
                                break
                            end
                        end
                    end
                end
                % determine if the current HEEL point is on a pressure
                a = [DynamicPPTrials(Selection(q)).R.LatIP{j}(1,1), DynamicPPTrials(Selection(q)).R.MedIP{j}(1,1)];
                b = [DynamicPPTrials(Selection(q)).R.LatIP{j}(1,2), DynamicPPTrials(Selection(q)).R.MedIP{j}(1,2)];
                SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                if sum(sum(Region)) > 0
                    PtOnOff = 'On';
                else
                    PtOnOff = 'Off';
                end
                if strcmp(PtOnOff,'On') == 1 % if the FP points are within the border of the foot, search from the Heel backward
                    for i = Loop % increment search points out from heel
                        a = [DynamicPPTrials(Selection(q)).R.LatIP{j}(1,1) - (i*xAdj), DynamicPPTrials(Selection(q)).R.MedIP{j}(1,1) - (i*xAdj)];
                        b = [DynamicPPTrials(Selection(q)).R.LatIP{j}(1,2) - (i*yAdj), DynamicPPTrials(Selection(q)).R.MedIP{j}(1,2) - (i*yAdj)];
                        SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                        SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) == 0
                            %                     subplot(1,NumDynTrials, Selection(q));
                            %                     line(a,b,'Color','k');
                            DynamicPPTrials(Selection(q)).R.LatIP{j}(1,:) = [a(1), b(1)];
                            DynamicPPTrials(Selection(q)).R.MedIP{j}(1,:) = [a(2), b(2)];
                            break
                        end
                    end
                else % if the FP points are outside the border of the foot, search from the heel towards the toe
                    for i = Loop % increment search points forward from heel
                        a = [DynamicPPTrials(Selection(q)).R.LatIP{j}(1,1) + (i*xAdj), DynamicPPTrials(Selection(q)).R.MedIP{j}(1,1) + (i*xAdj)];
                        b = [DynamicPPTrials(Selection(q)).R.LatIP{j}(1,2) + (i*yAdj), DynamicPPTrials(Selection(q)).R.MedIP{j}(1,2) + (i*yAdj)];
                        SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                        SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) > 0
                            %                     subplot(1,NumDynTrials, Selection(q));
                            %                     line(a,b,'Color','k');
                            DynamicPPTrials(Selection(q)).R.LatIP{j}(1,:) = [a(1), b(1)];
                            DynamicPPTrials(Selection(q)).R.MedIP{j}(1,:) = [a(2), b(2)];
                            break
                        end
                    end
                end
            end
            
            %% Heel to centroid points
            if Block.HC == 0
                % determine if the current TOE point is on a pressure
                xAdj = DynamicPPTrials(Selection(q)).Rht.VHeelCent{j}(:,1);
                yAdj = DynamicPPTrials(Selection(q)).Rht.VHeelCent{j}(:,2);
                a = [DynamicPPTrials(Selection(q)).R.LatHeelCent{j}(2,1), DynamicPPTrials(Selection(q)).R.MedHeelCent{j}(2,1)];
                b = [DynamicPPTrials(Selection(q)).R.LatHeelCent{j}(2,2), DynamicPPTrials(Selection(q)).R.MedHeelCent{j}(2,2)];
                SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                if sum(sum(Region)) > 0
                    PtOnOff = 'On';
                else
                    PtOnOff = 'Off';
                end
                if strcmp(PtOnOff,'On') == 1 % if the FP points are within the border of the foot, search from the toe forward
                    for i = Loop % increment search points out from toe
                        a = [DynamicPPTrials(Selection(q)).R.LatHeelCent{j}(2,1) + (i*xAdj), DynamicPPTrials(Selection(q)).R.MedHeelCent{j}(2,1) + (i*xAdj)];
                        b = [DynamicPPTrials(Selection(q)).R.LatHeelCent{j}(2,2) + (i*yAdj), DynamicPPTrials(Selection(q)).R.MedHeelCent{j}(2,2) + (i*yAdj)];
                        SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                        SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) == 0
                            %                     subplot(1,NumDynTrials, Selection(q));
                            %                     line(a,b,'Color','k');
                            DynamicPPTrials(Selection(q)).R.LatHeelCent{j}(2,:) = [a(1), b(1)];
                            DynamicPPTrials(Selection(q)).R.MedHeelCent{j}(2,:) = [a(2), b(2)];
                            break
                        end
                    end
                else % if the FP points are outside the border of the foot, search from the toe back towards the heel
                    for i = Loop % search forward to make sure there arent more pressures forward
                        a = [DynamicPPTrials(Selection(q)).R.LatHeelCent{j}(2,1) + (i*xAdj), DynamicPPTrials(Selection(q)).R.MedHeelCent{j}(2,1) + (i*xAdj)];
                        b = [DynamicPPTrials(Selection(q)).R.LatHeelCent{j}(2,2) + (i*yAdj), DynamicPPTrials(Selection(q)).R.MedHeelCent{j}(2,2) + (i*yAdj)];
                        SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                        SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) > 0
                            PressForward = 'Yes';
                            ForwardRegion = [a b];
                            break
                        else
                            PressForward = 'No';
                        end
                    end
                    if strcmp(PressForward,'No') == 1
                        for i = Loop % increment search points back from toe
                            a = [DynamicPPTrials(Selection(q)).R.LatHeelCent{j}(2,1) - (i*xAdj), DynamicPPTrials(Selection(q)).R.MedHeelCent{j}(2,1) - (i*xAdj)];
                            b = [DynamicPPTrials(Selection(q)).R.LatHeelCent{j}(2,2) - (i*yAdj), DynamicPPTrials(Selection(q)).R.MedHeelCent{j}(2,2) - (i*yAdj)];
                            SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                            SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                            Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                            Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                            if sum(sum(Region)) > 0
                                %                         subplot(1,NumDynTrials, Selection(q));
                                %                         line(a,b,'Color','k');
                                DynamicPPTrials(Selection(q)).R.LatHeelCent{j}(2,:) = [a(1), b(1)];
                                DynamicPPTrials(Selection(q)).R.MedHeelCent{j}(2,:) = [a(2), b(2)];
                                break
                            end
                        end
                    else % pressures are in front of the original toe point, keep searching forward
                        for i = Loop
                            a = [ForwardRegion(1) + (i*xAdj), ForwardRegion(2) + (i*xAdj)];
                            b = [ForwardRegion(3) + (i*yAdj), ForwardRegion(4) + (i*yAdj)];
                            SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                            SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                            Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                            Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                            if sum(sum(Region)) > 0
                                %                         subplot(1,NumDynTrials, Selection(q));
                                %                         line(a,b,'Color','k');
                                DynamicPPTrials(Selection(q)).R.LatHeelCent{j}(2,:) = [a(1), b(1)];
                                DynamicPPTrials(Selection(q)).R.MedHeelCent{j}(2,:) = [a(2), b(2)];
                                break
                            end
                        end
                    end
                end
                % determine if the current HEEL point is on a pressure
                a = [DynamicPPTrials(Selection(q)).R.LatHeelCent{j}(1,1), DynamicPPTrials(Selection(q)).R.MedHeelCent{j}(1,1)];
                b = [DynamicPPTrials(Selection(q)).R.LatHeelCent{j}(1,2), DynamicPPTrials(Selection(q)).R.MedHeelCent{j}(1,2)];
                SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                if sum(sum(Region)) > 0
                    PtOnOff = 'On';
                else
                    PtOnOff = 'Off';
                end
                if strcmp(PtOnOff,'On') == 1 % if the FP points are within the border of the foot, search from the Heel backward
                    for i = Loop % increment search points out from heel
                        a = [DynamicPPTrials(Selection(q)).R.LatHeelCent{j}(1,1) - (i*xAdj), DynamicPPTrials(Selection(q)).R.MedHeelCent{j}(1,1) - (i*xAdj)];
                        b = [DynamicPPTrials(Selection(q)).R.LatHeelCent{j}(1,2) - (i*yAdj), DynamicPPTrials(Selection(q)).R.MedHeelCent{j}(1,2) - (i*yAdj)];
                        SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                        SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) == 0
                            %                     subplot(1,NumDynTrials, Selection(q));
                            %                     line(a,b,'Color','k');
                            DynamicPPTrials(Selection(q)).R.LatHeelCent{j}(1,:) = [a(1), b(1)];
                            DynamicPPTrials(Selection(q)).R.MedHeelCent{j}(1,:) = [a(2), b(2)];
                            break
                        end
                    end
                else % if the FP points are outside the border of the foot, search from the heel towards the toe
                    for i = Loop % increment search points forward from heel
                        a = [DynamicPPTrials(Selection(q)).R.LatHeelCent{j}(1,1) + (i*xAdj), DynamicPPTrials(Selection(q)).R.MedHeelCent{j}(1,1) + (i*xAdj)];
                        b = [DynamicPPTrials(Selection(q)).R.LatHeelCent{j}(1,2) + (i*yAdj), DynamicPPTrials(Selection(q)).R.MedHeelCent{j}(1,2) + (i*yAdj)];
                        SrchAreaX = round([a(1)+Srch*xAdj; a(1)-Srch*xAdj; a(2)+Srch*xAdj; a(2)-Srch*xAdj],1);
                        SrchAreaY = round([b(1)+Srch*yAdj; b(1)-Srch*yAdj; b(2)+Srch*yAdj; b(2)-Srch*yAdj],1);
                        Zone = poly2mask(SrchAreaX,SrchAreaY,m,n);
                        Region = DynamicPPTrials(Selection(q)).MainMask .* Zone;
                        if sum(sum(Region)) > 0
                            %                     subplot(1,NumDynTrials, Selection(q));
                            %                     line(a,b,'Color','k');
                            DynamicPPTrials(Selection(q)).R.LatHeelCent{j}(1,:) = [a(1), b(1)];
                            DynamicPPTrials(Selection(q)).R.MedHeelCent{j}(1,:) = [a(2), b(2)];
                            break
                        end
                    end
                end
            end
        end
    end
end

end
