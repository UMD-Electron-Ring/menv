function f =  optfunc( X )

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
% another option: max(abs(betax-betay));

% (need to add condition for ref. traj)

% -- calculate min. function
try
f = [x(end)-runtmp.target.x1,...
    y(end)-runtmp.target.y1,...
    xp(end)-runtmp.target.xp1,...
    yp(end)-runtmp.target.yp1,...
    D(end)-runtmp.target.D1,...
    Dp(end)-runtmp.target.Dp1,...
    nux-runtmp.target.nuxt,...
    nuy-runtmp.target.nuyt,...
    env_diff];
catch % -- backwards compatibility; unspecified targets set to 0
   f = [x(end)-runtmp.target.x1,...
    y(end)-runtmp.target.y1,...
    xp(end)-runtmp.target.xp1,...
    yp(end)-runtmp.target.yp1,...
    0,0,0,0,0];
end

% -- apply weight factors
try
f = f.*[runtmp.weights.xw, runtmp.weights.yw,...
    runtmp.weights.xpw runtmp.weights.ypw,...
    runtmp.weights.Dw runtmp.weights.Dpw,...
    runtmp.weights.nuxw runtmp.weights.nuyw,...
    runtmp.weights.betaw];
catch % -- backwards compatibility; unspecified weights set to 0
    f = f.*[runtmp.weights.xw, runtmp.weights.yw,...
    runtmp.weights.xpw runtmp.weights.ypw,...
    0,0,0,0,0];
end
