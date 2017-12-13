function f =  periodfunc( X )
global kx ky irho ds;
global KX KY IRHO nsteps;
global xw yw xpw ypw;
global Dw Dpw

% Evaluate x0 y0 xp0 yp0
x0 = X(1);
y0 = X(2);
xp0= X(3);
yp0= X(4);
D0 = X(5);
Dp0 = X(6);

% Leap frog half step
kx = KX(1); ky = KY(1); irho = IRHO(1);
[xpp,ypp,Dpp] = calc_prim2(x0,y0,D0);
xp = xp0+xpp*ds/2;
yp = yp0+ypp*ds/2;
Dp(1) = Dp0+Dpp*ds/2;
x = x0;
y = y0;
D = D0;

% Steps
for i=1:nsteps-1
   kx = KX(i+1); ky = KY(i+1); irho = IRHO(i+1);
   [x,y,xp,yp,D,Dp] = step(x,y,xp,yp,D,Dp);
end;

% xp and yp back half step
[xpp,ypp,Dpp] = calc_prim2(x,y,D);
xp = xp - xpp*ds/2;
yp = yp - ypp*ds/2;
Dp = Dp - Dpp*ds/2;
f = [x-x0,y-y0,xp-xp0,yp-yp0,D-D0,Dp-Dp0];
f = f.*[xw yw xpw ypw Dw Dpw];
