function kappa = Gauss2Kappa( G, ele )
q = 1.602e-19;
M = 9.109e-31;
c = 2.9979e8;

E = 10.0e3;
T = q*E;
gamma = 1+T/(M*c^2);
beta = sqrt(1-gamma^-2);

if( nargin==1 | (nargin==2 & ele=='Q') )
   G = G*1e-2;    % G/cm --> T/m
   kappa = G / (gamma*M*beta*c/q);   
elseif( nargin==2 & ele=='S' )
   G = G*1e-4;    % G --> T
   kappa = ( G/(2*gamma*M*beta*c/q) )^2; 
else
   kappa = 0;
end;
