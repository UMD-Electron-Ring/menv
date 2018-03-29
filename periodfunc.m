function f =  periodfunc( X )

% Initial conditions based on optimization parameter X
ic = struct();
ic.x0 = X(1);
ic.y0 = X(2);
ic.xp0= X(3);
ic.yp0= X(4);
ic.D0 = X(5);
ic.Dp0 = X(6);

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
[x,y,xp,yp,D,Dp] = menv_integrator(Ex,Ey,K,KX,KY,IRHO,ic,ds,0);

% -- calculate min. function
f = [x-ic.x0,y-ic.y0,...
    xp-ic.xp0,yp-ic.yp0,...
    D-ic.D0,Dp-ic.Dp0];

% -- apply weight factors
try 
f = f.*[runtmp.weights.xw, runtmp.weights.yw,...
    runtmp.weights.xpw runtmp.weights.ypw,...
    runtmp.weights.Dw runtmp.weights.Dpw];
catch % -- if dispersion weights are not set (old version)
f = f.*[runtmp.weights.xw, runtmp.weights.yw,...
    runtmp.weights.xpw runtmp.weights.ypw,...
    0,0];
end
