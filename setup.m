
%  define constants, generate spectrum
%  one open channel assumed throughout

clear all
units
   
%  properties of the molecules

%  here we have NaRb
%ReducedMass = (110/2)/xm0;
%C6 = 1.525e6;        % atomic units; value for NaRb from Dajun's paper

%  here we have NaK
ReducedMass = (63/2)/xm0;
C6 = 561070;     % au, via Goulven

%  here we have Rb + KRb
%ReducedMass = 87*(87+40)/(87+87+40)/xm0;
%C6 = 8000;   % Koto, New J. Phys. 12, 073041 (2010).


% what are the natural units, in real units?
beta = (2*ReducedMass*C6)^(1/4);  % length scale, au
Ebeta = 1/2/ReducedMass/beta^2;   %  energy, au
taubeta = 2*pi/Ebeta;             %  time, au

Ebeta*t0;                        % energy, K
taubeta*tau0;                    % time, sec
abar = (pi/2^(3/2)/gamma(5/4)/gamma(1/2))^2;  %  GB scatt length

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  properties of the resonant spectrum 
MeanSpacing = 10.0;    %  units of E_beta
Meanx = 0.1;           % dimensionless

EnergySpectrumlo = -200*MeanSpacing;
EnergySpectrumhi = +200*MeanSpacing;

%  details of numerics
%rng('shuffle')
rng(1)

ET = 0;
Temp = 0.5;
numen = 200000;

EGrid = linspace(0,10*Temp,numen);
dGrid = EGrid(2)-EGrid(1);

numET = 1;

numSpect = 200;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MeanNu2 = Meanx*MeanSpacing/pi^2;
NumOpenChannels = 1;


  

   