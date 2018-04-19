function I = Kappa2Current( kappa, ele )

N = length(kappa(:));
I = zeros(1,N);

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
B_sbfact = sqrt(0.6945); % SB HE factor

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

I = kappa*0;
% -- loop over all inputs
for i=1:N
  
    if( nargin==1 || (nargin==2 && ele(i)=='Q') )
        G = ridg*kappa(i);
        I(i) = G/g0/q_sbfact;
        %(gamma*M*beta*c/q)/g0; % not sure what this is here for...
    elseif( nargin==2 && ele(i)=='S' )
        B = 2*ridg*sqrt(kappa(i))/B_sbfact; %Tesla
        B = B*1e4; %Gauss
        B = B-B_off;
        I(i) = B/B0;
    elseif( nargin==2 && ele(i)=='Y' )
        G = ridg*kappa(i);
        I(i) = G/g0_panofsky/panofsky_sbfact;
    elseif( nargin==2 && ele(i)=='D' )
        G = ridg*kappa(i);
        I(i) = G/g0_dipoy/d_sbfact;
    elseif( nargin==2 && ele(i)=='P' )
        G = ridg*kappa(i);
        I(i) = G/g0_pdy/pd_sbfact;
    end
    
end


