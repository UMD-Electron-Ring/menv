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

% -- quad model
g0 = 3.608e-2; % quad peak grad per amp
q_sbfact = 0.8354; % SB HE factor

% -- solenoid model
B0 = 17.597; % peak field per amp
B_off = 1.35; % offset
B_sbfact = sqrt(0.695); % SB HE factor

% -- loop over all inputs
for i=1:N
  
    if( nargin==1 || (nargin==2 && ele(i)=='Q') )
        G = (gamma*M*beta*c/q)*kappa(i);
        I(i) = G/g0/q_sbfact;
        %(gamma*M*beta*c/q)/g0; % not sure what this is here for...
    elseif( nargin==2 && ele(i)=='S' )
        B = 2*(gamma*M*beta*c/q)*sqrt(kappa(i))/B_sbfact; %Tesla
        B = B*1e4; %Gauss
        B = B-B_off;
        I(i) = B/B0;
    end
    
end


