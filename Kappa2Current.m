function I = Kappa2Current( kappa, ele )
q = 1.602e-19;
M = 9.109e-31;
c = 2.9979e8;

E = 10.0e3;
T = q*E;
gamma = 1+T/(M*c^2);
beta = sqrt(1-gamma^-2);

N = length(kappa);
I = zeros(1,N);

% -- loop over all inputs
for i=1:N
  
    if( nargin==1 || (nargin==2 && ele(i)=='Q') )
        g0 = 3.608e-2;
        G = (gamma*M*beta*c/q)*kappa(i);
        I(i) = G/g0;
        %(gamma*M*beta*c/q)/g0; % not sure what this is here for...
    elseif( nargin==2 && ele(i)=='S' )
        B = 2*(gamma*M*beta*c/q)*sqrt(kappa(i)); %Tesla
        B = B*1e4; %Gauss
        B = B-1.35;
        B0 = 17.6; %Gauss/Amp
        I(i) = B/B0;
    end
    
end


