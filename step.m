function [x,y,xp,yp] = step(x0,y0,xp0,yp0)
global ds;
x = x0+xp0*ds;
y = y0+yp0*ds;
[xpp,ypp] = calc_prim2(x,y);
xp = xp0+xpp*ds;
yp = yp0+ypp*ds;


