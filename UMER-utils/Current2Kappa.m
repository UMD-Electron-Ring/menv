function kappa = Current2Kappa( I, ele )
N = length(I(:));
kappa = zeros(1,N);

% -- some constants
q = 1.602e-19;
M = 9.109e-31;
c = 2.9979e8;

E = 10.0e3;
T = q*E;
gamma = 1+T/(M*c^2);
beta = sqrt(1-gamma^-2);
ridg = (gamma*M*beta*c/q); % beam rigidity

% MODEL PARAMETERS -------------------------------------------------------
% -- quad model
g0 = 3.608e-2; % quad peak grad per amp [T/m/A]
%q_sbfact = 0.8354; % RK HE factor, to be used w/ eff length 4.475 cm
q_sbfact = 0.7208; % SB HE factor (2006 note), to be used w/ eff length 5.164 cm

% -- solenoid model
B0 = 17.597; % peak field per amp [Gauss]
B_off = 1.35; % offset [Gauss]
B_sbfact = sqrt(0.6945); % SB HEfactor

% -- dipole model
leff_dipo = 3.819e-2; % HE dipole length [m]
g0_dipox = 0.146e-4/leff_dipo; % quad grad By in dipo (due to sext.component) [T/A]
g0_dipoy = 0.565e-4/leff_dipo; % quad grad Bx in dipo (due to edge-focusing component) [T/A]
d_sbfact  = 1; % Has not yet been calculated for dipoles.


% -- YQ/QR1 model
g0_panofsky = 1.01e-2; % quad peak grad per amp [T/m/A]
panofsky_sbfact = 0.7224; 

% -- PD model
leff_pd  = 5.006e-2; %HE PD length [m] 
g0_pdx = 0./leff_pd; % HE quad grad in PD (due to sext.component) has not been calc.
g0_pdy = 0.041e-4/leff_pd; % HE quad grad in PD (due to sext.component) [T/A]
pd_sbfact  = 1; % Has not yet been calculated for dipoles.
% ------------------------------------------------------------------------

kappa = 0*I;
for i=1:N

if( nargin==1 || (nargin==2 && ele(i)=='Q') )
   kappa(i) = (I(i)*g0)/ridg*q_sbfact;
elseif( nargin==2 && ele(i)=='S' )
	kappa(i) = ((I(i)*B0+B_off)*B_sbfact/(1e4*2*ridg))^2;
elseif( nargin==2 && ele(i)=='Y' )
    kappa(i) = (I(i)*g0_panofsky)/ridg*panofsky_sbfact;
elseif( nargin==2 && ele(i)=='D' )
    kappa(i) = (I(i)*g0_dipoy)/ridg*d_sbfact;
elseif( nargin==2 && ele(i)=='P' )
    kappa(i) = (I(i)*g0_pdy)/ridg*pd_sbfact;
else
   kappa(i) = 0;
   
end
end


