function [HTout] = FindIntersect(HTin, BoundingBox, ProgAng)
% Finds intersection points for heel and toe ends of the foot progression line

% Heel Line
yH = BoundingBox(2);
% Toe Line
yT = BoundingBox(2) + BoundingBox(4);
% Foot Prog Line Points
xFP(1) = HTin(1,1);
xFP(2) = HTin(2,1);
yFP(1) = HTin(1,2);
yFP(2) = HTin(2,2);
if abs(ProgAng) < 5 % if prog angle is steep, define points just as guestimate of intersection
    HTout(1,1) = xFP(1);
    HTout(1,2) = yH;
    HTout(2,1) = xFP(2);
    HTout(2,2) = yT;
else
    % Find starting heel point by line intersections
    pFP = polyfit(xFP, yFP,1);% define lines as linear function
    x = linspace(BoundingBox(1)-10,...  % create a series of x values from 10 below the bounding box min to 10 above the bounding box max
        BoundingBox(1)+BoundingBox(3)+10);
    yHval = yH * ones(1,100); % series of points at the same value = horizontal line
    yFPval = polyval(pFP, x); % evaluate the functions at the x values
    % Find heel points
    d = yHval ./ yFPval;  % divide each element of one function by the element of the other
    %         figure; subplot(1, NumDynTrials, Selection(q)); hold on;
    %         plot(x, yHval, 'm');
    %         plot(x, yFPval, 'm');
    ixH = x(d > 0.9 & d < 1.1); % when the quotient is near one, an intersection is present
    IntersectHFP(1) = mean(ixH);  % take the mean of the x intersection points
    iyFP = polyval(pFP,IntersectHFP(1));
    IntersectHFP(2) =  iyFP;  % average the y intersection points
    HTout(1,1) = IntersectHFP(1); % redefine the heel and toe points at the intersection
    HTout(1,2) = IntersectHFP(2);
    % find toe points
    yTval = yT * ones(1,100);
    d = yTval ./ yFPval;  % divide each element of one function by the element of the other
    %                 figure; subplot(1, NumDynTrials, Selection(q)); hold on;
    %                 plot(x, yTval, 'm');
    %                 plot(x, yFPval, 'm');
    ixT = x(d > 0.9 & d < 1.1); % when the quotient is near one, an intersection is present
    IntersectTFP(1) = mean(ixT);  % take the mean of the x intersection points
    iyFP = polyval(pFP,IntersectTFP(1));
    IntersectTFP(2) =  iyFP;  % average the y intersection points
    HTout(2,1) = IntersectTFP(1); % redefine the heel and toe points at the intersection
    HTout(2,2) = IntersectTFP(2);
end

if isnan(HTout(1,1)) ||  isnan(HTout(1,2)) 
    HTout(1,:) = HTin(1,:); 
end
if isnan(HTout(2,1)) ||  isnan(HTout(2,2)) 
    HTout(2,:) = HTin(2,:); 
end

end