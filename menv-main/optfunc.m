function f =  optfunc( X )

% -- this is a stop-gap measure until I figure out what better way to
% manage passing structure variables without hard-disk writes
global lattmp runtmp

% Update kappa for this step
menv_updatekappaarray(X)

% -- load lattice model parameters
%load 'lattmp'
IRHO = lattmp.IRHO;
KX = lattmp.KX;
KY = lattmp.KY;
loc1 = lattmp.loc1;
loc2 = lattmp.loc2;
OPT_ELE = lattmp.OPT_ELE;
nsteps = lattmp.nsteps;

% -- load runtmp (just use emitance and stepsize, in SI units)
%load 'runtmp'
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

% -- condition for symmetric beam
env_diff = rms(x-y);
% another option: max(abs(betax-betay));

% -- condition for reference trajectory
try
    refx = mean(abs(x-runtmp.target.yref));
    refy = mean(abs(y-runtmp.target.yref));
catch
    [refx,refy] = deal(0);
end

% -- calculate min. function
try
f = [x(end)-runtmp.target.x1,...
    y(end)-runtmp.target.y1,...
    xp(end)-runtmp.target.xp1,...
    yp(end)-runtmp.target.yp1,...
    nux-runtmp.target.nux1,...
    nuy-runtmp.target.nuy1,...
    env_diff,refx,refy];
catch % -- backwards compatibility; Dispersion target set to 0
   f = [x(end)-runtmp.target.x1,...
    y(end)-runtmp.target.y1,...
    xp(end)-runtmp.target.xp1,...
    yp(end)-runtmp.target.yp1,...
    0,0,env_diff,refx,refy];
    %warning('No dispersion or tune targets set')
end

% -- apply weight factors to min. function vector
try
f = f.*[runtmp.weights.xw, runtmp.weights.yw,...
    runtmp.weights.xpw runtmp.weights.ypw,...
    runtmp.weights.nuxw runtmp.weights.nuyw,...
    runtmp.weights.betaw,runtmp.weights.refw,runtmp.weights.refw];
catch % -- backwards compatibility; unspecified weights set to 0
    f = f.*[runtmp.weights.xw, runtmp.weights.yw,...
    runtmp.weights.xpw runtmp.weights.ypw,...
    0,0,0,0,0];
    %warning('Missing some weight factors, only matching to x,y,x'',y'' target')
end

% -- save some vars in runtmp
runtmp.f = f;
runtmp.Xstr = X;
%save 'runtmp' runtmp;
