function kappa = Current2Kappa( I, ele )
N = length(I);
kappa = zeros(1,N);

% -- some constants
q = 1.602e-19;
M = 9.109e-31;
c = 2.9979e8;

E = 10.0e3;
T = q*E;
gamma = 1+T/(M*c^2);
beta = sqrt(1-gamma^-2);

% MODEL PARAMETERS -------------------------------------------------------
% -- quad model
g0 = 3.608e-2; % quad peak grad per amp [T/m/A]
%q_sbfact = 0.8354; % RK HE factor, to be used w/ eff length 4.475 cm
q_sbfact = 0.7208; % SB HE factor (2006 note), to be used w/ eff length 5.164 cm

% -- solenoid model
B0 = 17.597; % peak field per amp [Gauss]
B_off = 1.35; % offset [Gauss]
B_sbfact = sqrt(0.6945); % SB HE factor

% -- dipole model
gint_dipox = 0.146e-4; % int. quad grad By in dipo (due to sext.component) [T/A]
gint_dipoy = 0.565e-4; % int. quad grad Bx in dipo (due to edge-focusing component) [T/A]
d_sbfact  = 1; % Has not yet been calculated for dipoles.

% -- YQ/QR1 model
g0_panofsky = 1.01e-2; % quad peak grad per amp [T/m/A]
panofsky_sbfact = 0.7224; 

% -- PD model
gint_pdx = 0.; % int. quad grad in PD (due to sext.component) has not been calc.
gint_pdy = 0.041e-4; % int. quad grad in PD (due to sext.component) [T/A]
pd_sbfact  = 1; % Has not yet been calculated for dipoles.
% ------------------------------------------------------------------------

for i=1:N

if( nargin==1 || (nargin==2 && ele(i)=='Q') )
   kappa(i) = (I(i)*g0)/(gamma*M*beta*c/q)*q_sbfact;
elseif( nargin==2 && ele(i)=='S' )
	kappa(i) = ((I(i)*B0+B_off)*B_sbfact/(1e4*2*(gamma*M*beta*c/q)))^2;
elseif( nargin==2 && ele(i)=='Y' )
    kappa(i) = (I(i)*g0_panofsky)/(gamma*M*beta*c/q)*panofsky_sbfact;
elseif( nargin==2 && ele(i)=='D' )
    kappa(i) = (I(i)*gint_dipoy)/(gamma*M*beta*c/q)*d_sbfact;
elseif( nargin==2 && ele(i)=='P' )
    kappa(i) = (I(i)*gint_pdy)/(gamma*M*beta*c/q)*pd_sbfact;
else
   kappa(i) = 0;
   
end
end


