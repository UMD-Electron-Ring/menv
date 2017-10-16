function f =  stepfunc( X )
global kx ky irho ds;
global loc1 loc2 KX KY IRHO nsteps;
global x0 y0 xp0 yp0 D0 Dp0;
global x1 y1 xp1 yp1 D1 Dp1;
global xw yw xpw ypw Dw Dpw;
global OPT_ELE;

% Evaluate kappa
for i=1:length( X )
    if OPT_ELE(i)=='S'
        KX( loc1(i):loc2(i) ) = X(i); %0.96891*X(i);
        KY( loc1(i):loc2(i) ) = X(i); %-X(i);
    elseif OPT_ELE(i)=='Q'
        KX( loc1(i):loc2(i) ) = 0.955*X(i);
        KY( loc1(i):loc2(i) ) = -0.935*X(i);
    elseif ele(i)=='D'
        KX( loc1(i):loc2(i) ) = 1.9687*X(i)*(1-dipl_n(i));
        KY( loc1(i):loc2(i) ) = X(i)*dipl_n(i);
    end;
end;

% Leap-frog by half step
kx = KX(1); ky = KY(1); irho = IRHO(1);
[xpp,ypp,Dpp] = calc_prim2(x0,y0,D0);
xp = xp0+xpp*ds/2;
yp = yp0+ypp*ds/2;
Dp = Dp0+Dpp*ds/2;
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
f = [x-x1,y-y1,xp-xp1,yp-yp1,D-D1,Dp-Dp1];
f = f.*[xw yw xpw ypw Dw Dpw];
