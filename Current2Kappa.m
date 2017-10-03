function kappa = Current2Kappa( I, ele )
q = 1.602e-19;
M = 9.109e-31;
c = 2.9979e8;

E = 10.0e3;
T = q*E;
gamma = 1+T/(M*c^2);
beta = sqrt(1-gamma^-2);

if( nargin==1 | (nargin==2 & ele=='Q') )
   g0 = 3.608e-2;
   kappa = (I*g0)/(gamma*M*beta*c/q)
elseif( nargin==2 & ele=='S' )
	B0 = 17.6; %Gauss/Amp
	kappa = ((I*B0+1.35)/(1e4*2*(gamma*M*beta*c/q)))^2

else
   I = 0;
end;


