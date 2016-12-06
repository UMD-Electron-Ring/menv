function f =  periodfunc( X )
global kx ky ds;
global KX KY nsteps;
global xw yw xpw ypw;

% Evaluate x0 y0 xp0 yp0
x0 = X(1);
y0 = X(2);
xp0= X(3);
yp0= X(4);

% Leap frog half step
kx = KX(1); ky = KY(1);
[xpp,ypp] = calc_prim2(x0,y0);
xp = xp0+xpp*ds/2;
yp = yp0+ypp*ds/2;
x = x0;
y = y0;

% Steps
for i=1:nsteps-1
   kx = KX(i+1); ky = KY(i+1);
   [x,y,xp,yp] = step(x,y,xp,yp);
end;

% xp and yp back half step
[xpp,ypp] = calc_prim2(x,y);
xp = xp - xpp*ds/2;
yp = yp - ypp*ds/2;
f = [x-x0,y-y0,xp-xp0,yp-yp0];
f = f.*[xw yw xpw ypw];
