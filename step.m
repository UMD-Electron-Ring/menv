function [x,y,xp,yp,D,Dp] = step(x0,y0,xp0,yp0,D0,Dp0)
global ds;
x = x0+xp0*ds;
y = y0+yp0*ds;
D = D0+Dp0*ds;
[xpp,ypp,Dpp] = calc_prim2(x,y,D);
xp = xp0+xpp*ds;
yp = yp0+ypp*ds;
Dp = Dp0+Dpp*ds;

