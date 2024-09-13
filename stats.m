%  get a bunch of time delays and do statistics

clear ETmat QTmat

ETmat = (rand(1,numET)-0.5)*(EnergySpectrumhi-EnergySpectrumlo)/2;
numpos = 0;
for iSpect = 1: numSpect
    % go get a spectrum
    clear NumResonances EGOEmat Wmat EResmat GammaResmat ximat
    [ NumResonances, EGOEmat, Wmat, EResmat, GammaResmat, ximat ] = ...
    GenerateSpectrum( EnergySpectrumlo, EnergySpectrumhi, MeanSpacing, ...
                               MeanNu2, NumOpenChannels );
    for iET = 1 : length(ETmat)
        ET = ETmat(iET);
        QT = Thermal_Q(Temp, ET, ...
                      abar, EResmat, Wmat, ...
                      numen);
        if QT > 0
            numpos = numpos + 1;
            QTmat(numpos) = QT;
        end
    end
end
%  important step - normalize to RRKM lifetime
QTmat = (MeanSpacing/2/pi)*QTmat;
FractionPositive = numpos/(numET*numSpect)
MeanQ = mean(QTmat);
STDQ = std(QTmat);

QT10mat = log10(QTmat);
MeanQ10 = mean(QT10mat)
STDQ10 = std(QT10mat)

%figure(2)
%histogram(QTmat)
%xlabel('QT/TRRKM')
%ylabel('counts')
%hold off

figure(3)
Nbins=8;
histogram(QT10mat,Nbins)
xlabel("$\log_{10} (\tau / \tau_{RRKM})$",'Interpreter','latex','FontSize',20)
ylabel('Counts','Fontsize',20)
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  functions needed
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ NumResonances, EGOEmat, Wmat, EResmat, GammaResmat, ximat ] = ...
             GenerateSpectrum( Elo, Ehi, MeanSpacing, ...
                                MeanNu2, NumOpenChannels )

%  spectrum of resonances 
%  given: 
%         Elo = lower limit of spectrum
%         Ehi = upper limit of spectrum
%         MeanSpacing = mean spacing 
%         MeanW2 = mean value of coupling strength, ...
%                             <W_{mu,i}W_{nu,j}> = MeanW2
%         NumOpenChannels = what you think it is
%  on output:
%         NumResonances = number of resonances in this iteration
%         EGOEmat(1:NumResonances) = set of energies, ...
%                                   distributed according to GOE
%         Wmat(1:NumResonances,1:NumOpenChannels) = coupling matrix
%         EResmat(1:NumResonances)
%         GammaResmat(1:NumResonances),
%              where wigenvalue of the lambda-th resonance is
%              EResmat(lambda) - (1/2) GammaResmat
%         ximat(1:NumOpenChannels) = x_i parameter in channel i
%            Note: x is probably given a target value,  but the output is
%            based on the actual mean coupling in the numerical spectrum
%         


%Wigner-Dyson-distributed random numbers
%   based on a mapping from a uniform distribution

NumResonances = ceil((Ehi-Elo)/MeanSpacing);
x = rand(1,NumResonances);
% rescale the uniform distribution thusly:
spacing = (-(4/pi)*log(1-x)).^(1/2);
spacing = spacing*MeanSpacing;

% spacings determine the spectrum! these will be listed in increasing order
EGOEmat = zeros(NumResonances,NumResonances);
EGOEmat(1,1)=Elo;
for mu = 2: NumResonances
    EGOEmat(mu,mu) = EGOEmat(mu-1,mu-1)+spacing(mu);
end

%  generate coupling matrix elements
for mu = 1: NumResonances
    for i = 1: NumOpenChannels
        Wmat(mu,i) = normrnd(0,sqrt(MeanNu2));   % b/c MeanNu2 is the variance
    end
end

%    next get eigenenergies of the GOE Hamiltonina
Heff = EGOEmat - 1i*pi*Wmat*transpose(Wmat);
[V, D] = eig(Heff);
[d, ind] = sort(real(diag(D)));
Ds = D(ind,ind);
Vs = V(:,ind);
transpose(Vs)*Vs;
EResmat = diag(real(Ds),0);         %  complex energy is
GammaResmat = -2*diag(imag(Ds),0);  %  ERes - (i/2)GammaRes


%  dimensionless coupling constants, "var" is variance of the columns
ximat = (pi^2/MeanSpacing)*var(Wmat);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [QT] = Thermal_Q(Temp, ET, ...
                          abar, EResmat, Wmat, ...
                          numen)

%  Thermally averaged time delay 

%  Temp = temperature
%  ET = location in spectrum
% numen = number of energies in termperature integral

EGrid = linspace(0,10*Temp,numen);
dGrid = EGrid(2)-EGrid(1);
prod = 0.0;
sumGrid = 0;
for iGrid = 2: length(EGrid)
    Energy = EGrid(iGrid);
    prod = RelativeMB(Temp, Energy) ...
           * Time_Delay(Energy, ET, abar, ...
                        EResmat, Wmat);
    sumGrid = sumGrid + prod * dGrid;
end
QT = sumGrid;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Q = Time_Delay(Energy, ET, abar, ...
                                EResmat, Wmat)

% time delay 
%   Energy = collision energy
%   ET = shift of spectrum

%   Y = y-matrix,  Yp = its energy derivative

k = sqrt(Energy);
A = abar*k;
Ap = abar/2/k;
G = (1/3-abar^2)*k^2;
Gp = (1/3-abar^2);
eta = -abar*k;
etap = -abar/2/k;

Y = 0;
Yp = 0;
for mu = 1 : length(EResmat)
    Y = Y + Wmat(mu)^2/((Energy-ET)-EResmat(mu));
    Yp = Yp + Wmat(mu)^2/((Energy-ET)-EResmat(mu))^2;
end

Y = -pi*Y;
Yp = pi*Yp;

K = A*Y/(1+G*Y);
Kp = (Ap*Y + A*Yp)/(1 + G*Y) - A*Y*(Gp*Y + G*Yp)/(1+G*Y)^2;

Q = 2*( etap + Kp/(1+K^2) );
%Q = 2*Yp/(1+Y^2);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ fr ] = RelativeMB( T, E )
%  Maxwell-Boltzmann distribution for the relative motion

beta = 1/T;  % assumes that kT is given in the energy units of E
fr = ( 2*beta^(3/2)/pi^(1/2) )*E.^(1/2).*exp(-beta*E);

end


