function f =  stepfuncMO( X )

load 'runtmp';

KX = runtmp.KX;
KY = runtmp.KY;
nsteps = round(runtmp.distance/runtmp.stepsize) +1;

% Evaluate kappa
for i=1:length( X )
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
    end
end;

% -- integrate
ic = [runtmp.x0,runtmp.y0,runtmp.xp0,runtmp.yp0]; % initial conditions
[x,y,xp,yp] = envintegrator(KX,KY,ic,runtmp.stepsize,nsteps,1);

% -- calculate tune values
Lchan = 32;
d = 0:runtmp.stepsize:runtmp.distance;

betax = (x).^2/(runtmp.emitance);
betay = (y).^2/(runtmp.emitance);

% -- full lattice tune
psix = runtmp.stepsize*sum(1./betax);
psiy = runtmp.stepsize*sum(1./betay);
nux = psix/(2*pi);
nuy = psiy/(2*pi);

% condition for integer ring tune
nuxcontr = nux-runtmp.nuxt;
nuycontr = nuy-runtmp.nuyt;

% condition for envelope target
xcontr = x(end)-runtmp.x1;
ycontr = y(end)-runtmp.y1;
xpcontr = xp(end)-runtmp.xp1;
ypcontr = yp(end)-runtmp.yp1; 

% condition for symmetric beam
betacontr = max(abs(betax-betay));

f1 = [xcontr,xpcontr];
f2 = [ycontr,ypcontr];
f3 = [nuxcontr,nuycontr];
f1 = f1.*[runtmp.xw runtmp.xpw ];
f2 = f2.*[runtmp.yw runtmp.ypw ];
f3 = f3.*[runtmp.nuxw runtmp.nuyw];
f(1) = sum(abs(f1));
f(2) = sum(abs(f2));
f(3) = sum(abs(f3));
