% calculate reduced gravity, phase speed and length of soliton
% based on two-layer KdV solution

% parameters you input
rho1=eqstate(20,23); % (temp, sali)
rho2=eqstate(20,24);
h1=100; % upper layer thickness in meters
h2=900; % lower layer thickness in meters
a=20.0; % amplitude of soliton

% reduced gravity
g_reduced=9.81*abs(rho1-rho2)/rho1

% phase speed from linear solution
c0=sqrt(g_reduced*h1*h2/(h1+h2))

% phase speed from nonlinear solution
c=c0*(1.0+0.5*a*abs(h1-h2)/(h1*h2)) % h1>h2 is a positive soliton

% wave length
l=2.0*h1*h2/sqrt(3.0*a*abs(h1-h2))

