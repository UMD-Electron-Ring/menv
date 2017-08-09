function I = Kappa2Current( kappa, ele )
q = 1.602e-19;
M = 9.109e-31;
c = 2.9979e8;

E = 10.0e3;
T = q*E;
gamma = 1+T/(M*c^2);
beta = sqrt(1-gamma^-2);

if( nargin==1 | (nargin==2 & ele=='Q') )
   g0 = 4.14e-2;
   G = (gamma*M*beta*c/q)*kappa;
   I = G/g0;
   (gamma*M*beta*c/q)/g0
elseif( nargin==2 & ele=='S' )
   B = 2*(gamma*M*beta*c/q)*sqrt(kappa); %Tesla
   B = B*1e4; %Gauss
   B = B-1.35;
   B0 = 17.6; %Gauss/Amp
   I = B/B0;
else
   I = 0;
end;


