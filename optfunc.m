function f =  optfunc( X )
global clm

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

% -- calculate other quantities

% full lattice tune
betax = (x).^2/(runtmp.emitance);
betay = (y).^2/(runtmp.emitance);
psix = runtmp.stepsize*sum(1./betax);
psiy = runtmp.stepsize*sum(1./betay);
nux = psix/(2*pi);
nuy = psiy/(2*pi);

% condition for symmetric beam
env_diff = rms(x-y);

% (need to add condition for ref. traj)

% -- calculate min. function
try
f = [x(end)-clm.usrdata.target.x1,...
    y(end)-clm.usrdata.target.y1,...
    xp(end)-clm.usrdata.target.xp1,...
    yp(end)-clm.usrdata.target.yp1,...
    D(end)-clm.usrdata.target.D1,...
    Dp(end)-clm.usrdata.target.Dp1,...
    nux-clm.usrdata.target.nuxt,...
    nuy-clm.usrdata.target.nuyt,...
    env_diff];
catch % -- backwards compatibility; unspecified targets set to 0
   f = [x(end)-clm.usrdata.target.x1,...
    y(end)-clm.usrdata.target.y1,...
    xp(end)-clm.usrdata.target.xp1,...
    yp(end)-clm.usrdata.target.yp1,...
    0,0,0,0,0];
end

% -- apply weight factors
try
f = f.*[clm.usrdata.weights.xw, clm.usrdata.weights.yw,...
    clm.usrdata.weights.xpw clm.usrdata.weights.ypw,...
    clm.usrdata.weights.Dw clm.usrdata.weights.Dpw,...
    clm.usrdata.weights.nuxw clm.usrdata.weights.nuyw,...
    clm.usrdata.weights.betaw];
catch % -- backwards compatibility; unspecified weights set to 0
    f = f.*[clm.usrdata.weights.xw, clm.usrdata.weights.yw,...
    clm.usrdata.weights.xpw clm.usrdata.weights.ypw,...
    0,0,0,0,0];
end
