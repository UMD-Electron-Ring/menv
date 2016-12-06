function G = Kappa2Gauss( kappa, ele )
q = 1.602176462e-19;
M = 9.10938188e-31;
c = 2.99792458e8;

E = 10.0e3;
T = q*E;
gamma = 1+T/(M*c^2);
beta = sqrt(1-gamma^-2);

if( nargin==1 | (nargin==2 & ele=='Q') )
   G = (gamma*M*beta*c/q)*kappa; %Tesla/m
   G = G*1e2;                    %T/m --> G/cm
elseif( nargin==2 & ele=='S' )
   G = 2*(gamma*M*beta*c/q)*sqrt(kappa); %Tesla
   G = G*1e4;                            % T --> G
else
   G = 0;
end;
