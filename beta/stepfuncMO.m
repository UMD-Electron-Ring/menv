function f =  stepfuncMO( X )

% Update kappa for this step
menv_updatekappaarray(X)

% -- load lattice model parameters
load 'lattmp'
IRHO = lattmp.IRHO;
KX = lattmp.KX;
KY = lattmp.KY;
loc1 = lattmp.loc1;
loc2 = lattmp.loc2;
OPT_ELE = lattmp.OPT_ELE;
nsteps = lattmp.nsteps;

% -- load runtmp (just use emitance and stepsize, in SI units)
load 'runtmp'
K = runtmp.perveance; % Perveance
Ex = runtmp.emitance; % Emmitance x
Ey = Ex;
ds = runtmp.stepsize;

% -- integrate envelope
ic = runtmp.ic;
[x,y,xp,yp,D,Dp] = menv_integrator(Ex,Ey,K,KX,KY,IRHO,ic,ds,1);


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
