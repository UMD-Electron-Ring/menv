function f =  optfunc( X )
global clm

% Update kappa for this step
match_updatekappaarray(X)

% -- load lattice model parameters
load 'lattmp'
IRHO = lattmp.IRHO;
KX = lattmp.KX;
KY = lattmp.KY;
loc1 = lattmp.loc1;
loc2 = lattmp.loc2;
OPT_ELE = lattmp.OPT_ELE;
nsteps = lattmp.nsteps;
ds = lattmp.ds;
K = clm.usrdata.perveance; % Perveance
Ex = clm.usrdata.emitance; % Emmitance x
Ey = Ex;
ds = clm.usrdata.stepsize;

% -- integrate envelope
ic = [clm.usrdata.x0 clm.usrdata.y0 clm.usrdata.xp0 clm.usrdata.yp0 clm.usrdata.D0 clm.usrdata.Dp0];
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
env_diff = x-y;

% (need to add condition for ref. traj)

% -- calculate min. function
f = [x-clm.usrdata.x1,...
    y-clm.usrdata.y1,...
    xp-clm.usrdata.xp1,...
    yp-clm.usrdata.yp1,...
    D-clm.usrdata.D1,...
    Dp-clm.usrdata.Dp1
    nux-clm.usrdata.nuxt,...
    nuw-clm.usrdata.nuyt,...
    env_diff];


% -- apply weight factors
f = f.*[clm.usrdata.xw, clm.usrdata.yw,...
    clm.usrdata.xpw clm.usrdata.ypw,...
    clm.usrdata.Dw clm.usrdata.Dpw,...
    clm.usrdata.nuxw clm.usrdata.nuyw,...
    clm.usrdata.betaw];
