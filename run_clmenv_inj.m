global clm
clm = clmenv();

% -- load .spt file (optional)
%clm.open('FODO_0_7mA_dipl.spt')

%
% -- set params (over-writes spt file data if it's loaded)
%

% Beam parameters
emittance = 30; % mm-mrad
perveance = 0.0001045;
x0        = .3; % cm
y0        = .3; % cm
xp0       = 0.0;
yp0       = 0.0;

% Simulation parameters
stepsize = 0.005; % cm
distance = 32; % cm

% optimization parameters
iterations = 20;
tolerance  = 1e-6; 

%-- Lattice description
elements ='QDQ'; % element types
location = [8,16,24]; % element location
lengths = [3.6400,3.8500,3.6400]; % element length
str = [220,15.7000,-220]; % quadrupole strength (kappa)
did = [0 0.7200 0]; % dipole index
opt = [0,0,0];

% -- assign params to menv structures
clm.maketarget([1,1,1,1],[1,1,1,1])
clm.makeoptiset(iterations,tolerance)  % Makes file optiset
clm.defmatcher() % loads matcher settings

ic = struct();
ic.x0 = x0; ic.y0=y0; ic.xp0=xp0; ic.yp0=yp0; ic.D0 =0 ; ic.Dp0 = 0;
clm.makeparams(emittance,perveance,ic,stepsize,distance) % makes param file
clm.defparam() % loads params

clm.deflattice(elements,location,lengths,str,opt,did) % loads lattice

% -- look at initial solution
clm.solve()


%%
% -- solve for periodic solution
clm.periodicmatcher()

% -- save matching parameters in temp variable
FODOmatch = clm.soldata;

%%
% -- load injection .spt file (to get lattice settings)
clm.open('INJ_minsize_100mA_dipl.spt')

% -- beam parameters at aperture
x0aper = 0.3;
y0aper = 0.3;
xp0aper = 0;
yp0aper = 0;
icaper = struct();
icaper.x0 = x0aper; icaper.y0=y0aper; icaper.xp0=xp0aper; icaper.yp0=yp0aper;

inj_distance = clm.usrdata.distance; % length of injection line sim

clm.makeparams(emittance,perveance,icaper,stepsize,inj_distance) % makes param file
clm.defparam()

clm.makeoptiset(iterations,tolerance) 
clm.maketarget([FODOmatch.xf,FODOmatch.yf,FODOmatch.xpf,FODOmatch.ypf],[1,1,1,1])
clm.defmatcher()

% -- should also change lattice parameters here
elements = clm.usrdata.ele;
location = clm.usrdata.loc;
lengths = clm.usrdata.len;
str = clm.usrdata.str; % quadrupole strength (kappa)
did = clm.usrdata.did; % dipole index
opt = clm.usrdata.opt;



% -- look at initial solution
clm.solve()

%%

% -- match to target x,y,x',y'
clm.targetmatcher()


% -- save new magnet settings to variable
newlattice = clm.usrdata.str