
function Rho = eqstate(T,S)
%  Fenygan Shi
%  2018/12/15
%
%  This script calculates the density of seawater using the Knudsen-Ekman
%  equation of state. It is consistent with the formula used in NHWAVE
%
%  usage: Rho = eqstate(T,S)
%  T is temperature DEG
%  S is salinity PSU
%  Rho is output of density
%
%  Reference: Fofonoff, N., The Sea: Vol 1, (ed. M. Hill). Interscience,
%  New York, 1962, pp 3-30.
%
      TF = T;
      SF = S;
      Rho0 = 1000.0;
      RHOF = SF*SF*SF*6.76786136*10^(-6)-SF*SF*4.8249614*10^(-4)+ ...
           SF*8.14876577&10^(-1)-0.22584586*10^(0);
      RHOF = RHOF*(TF*TF*TF*1.667*10^(-8)-TF*TF*8.164*10^(-7)+ ...
           TF*1.803*10^(-5));
      RHOF = RHOF+1.-TF*TF*TF*1.0843*10^(-6)+TF*TF*9.8185*10^(-5)-TF*4.786*10^(-3);
      RHOF = RHOF*(SF*SF*SF*6.76786136*10^(-6)-SF*SF*4.8249614*10^(-4)+ ...
           SF*8.14876577*10^(-1)+3.895414*10^(-2));
      RHOF = RHOF-(TF-3.98)^2*(TF+283.0)/(503.57*(TF+67.26));
      Rho=Rho0+RHOF;