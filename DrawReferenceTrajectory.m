function [x,y] = DrawReferenceTrajectory(axeshandle)

% This function allows user to draw a set of reference trajectory points in
% the figure window. If you do not draw all the way to the edge of the
% figure, this function continues the last line segment to the edge (same
% slope)
%
% 12/7/16 Kiersten Ruisard

% Prompt to draw line
[x,y] = getline(axeshandle);
xl = xlim(axeshandle);

% Continue left edge
if x(1) > xl(1)
    slope = (y(2)-y(1))/(x(2)-x(1));
    y0 = y(1)-slope*x(1);
    x = [xl(1);x]; 
    y = [y0;y];
end

% Continue right edge
if x(end) < xl(2)
    slope = (y(end)-y(end-1))/(x(end)-x(end-1));
    ye = y(end)+slope*(xl(2)-x(end));
    x = [x;xl(2)];
    y = [y;ye];
end

end % function end
