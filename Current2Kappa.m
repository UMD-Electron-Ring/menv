function kappa = Current2Kappa( I, ele )
q = 1.602e-19;
M = 9.109e-31;
c = 2.9979e8;

E = 10.0e3;
T = q*E;
gamma = 1+T/(M*c^2);
beta = sqrt(1-gamma^-2);

% -- quad model parameters
g0 = 3.608e-2; % quad peak grad per amp
%q_sbfact = 0.8354; % SB HE factor
q_sbfact = 0.7208; % SB HE factor updated (2006)

% -- solenoid model parameters
B0 = 17.597; % peak field per amp
B_off = 1.35; % offset
B_sbfact = sqrt(0.6945); % SB HE factor

N = length(I);
kappa = zeros(1,N);

for i=1:N

if( nargin==1 || (nargin==2 && ele(i)=='Q') )
   kappa(i) = (I(i)*g0)/(gamma*M*beta*c/q)*q_sbfact;
elseif( nargin==2 && ele(i)=='S' )
	kappa(i) = ((I(i)*B0+B_off)*B_sbfact/(1e4*2*(gamma*M*beta*c/q)))^2;

else
   kappa(i) = 0;
   
end
end


