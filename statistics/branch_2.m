%  first try: probability into and out of complex
%  given by overlap integrals squared; J conserved

clear all %  let's see python do that!


% fixed q. numbers (doubled!)
iR1 = 3; % rubidium 1 nuclear spin
iR2 = 3; % rubidium 1 nuclear spin
iK = 8; % potassium nuclear spin
s = 1; % potassium electronic spin

%  initial state q.  numbers  (doubled!)
Np = 0; % molecule rotation
MNp = 0;
Lp = 0; % partial wave
MLp = 0;
fp = 2; % total atom spin
mfp = 0;
mR2p = -3; % molecule rubidium nuclear spin projection
mKp = -8; % potassium nuclear spin projection

MF = MNp + MLp + mfp + mR2p + mKp;

%  generate states of the complex,
%  and the probabilities to reach them

F_complex = 1:2:15;  % hard-wired for a given initial
p_form = zeros(1, length(F_complex));  
icomplex = 0;
%for iR = 0: 2 : iR1 + iR2
for iR = 0: 2 : iR1 + iR2   % only even iR - fermions!
    for I = abs(iK-iR) : 2 : iK + iR
        for P = abs(I-s) : 2 : I+s
            for J = abs(Np-Lp) : 2 : Np+Lp                
                for F = abs(P-J) : 2 : P+J
                    if F >= abs(MF)
                        % iR = 0, mR1 = -1, I = 8, P = 7, J = 0, F = 7
                        % ms = 1, mI = 0, mP = -7, mJ = 0, mF = -7, f = 1
                        % [0 8 7 0 7]
                        icomplex = icomplex + 1;
                        iF = (F+1)/2;
                        % disp([iR, I, P, J, F])
                        xx = ...
                            overlap(iR1, iR2, iK, s, ...
                                     Np, MNp, Lp, MLp, fp, mfp, mR2p, mKp, ...
                                     J, iR, I, P, F, MF );
                        % disp(xx)
                        complex(icomplex,1) = icomplex;
                        complex(icomplex,2) = iR;
                        complex(icomplex,3) = I;
                        complex(icomplex,4) = P;
                        complex(icomplex,5) = J;
                        complex(icomplex,6) = F;
                        %p_complex_form(icomplex) = xx^2;
                        %p_form = append(p_form,F_complex,F,xx^2);
                        p_form(iF) = p_form(iF) + xx^2;
                    end
                end
            end
        end
    end
end
n_complex = icomplex   

%% 

%  tabulate possible final states
%  not all will actually project onto complex,
%  but never mind
Nmax = 0;  % restricted in exit channel by energy
fexit = 2; % final hyperfine is f=1
mfexit = 2;  % this is known actually
iexit = 0;
for N = 0 : 2 : Nmax
    for MN = -N : 2 : N
         if N==0 
             Lmin = 0;
             Lmax = 30;
             %Lmax = 60;
         end
         for L = Lmin : 4 : Lmax % parity restricts L to even/odd
             for ML = -L : 2 : L
                 for mf = mfexit : mfexit
                 %for mf = -fexit : 2 : fexit
                     for mR2 = -iR2 : 2 : iR2
                         for mK = -iK : 2 : iK 
                             if MN + ML + mf + mR2 + mK == MF
                                iexit = iexit + 1;
                                exit(iexit,1) = iexit;     
                                exit(iexit,2) = N;
                                exit(iexit,3) = MN;
                                exit(iexit,4) = L;
                                exit(iexit,5) = ML;
                                exit(iexit,6) = fexit;
                                exit(iexit,7) = mf;
                                exit(iexit,8) = mR2;
                                exit(iexit,9) = mK;
                             end
                         end
                     end
                 end
             end
         end
    end
end
numexit = iexit  % there are lots!  lots!

p_exit = zeros(numexit, length(F_complex) );
% for each of the exit channels, get overlap to complex states
for iexit = 1 : numexit
    if mod(iexit,100) == 0
        disp(iexit)
    end
    N = exit(iexit,2);
    MN = exit(iexit,3);
    L = exit(iexit,4);
    ML = exit(iexit,5);
    f = exit(iexit,6);
    mf = exit(iexit,7);
    mR2 = exit(iexit,8);
    mK = exit(iexit,9);
    %for iR = 0: 2 : iR1 + iR2
    for iR = 0: 4 : iR1 + iR2   % only even iR - fermions!
        for I = abs(iK-iR) : 2 : iK + iR
            for P = abs(I-s) : 2 : I+s
                for J = abs(N-L) : 2 : N+L                
                    for F = abs(P-J) : 2 : min((P+J),15)  %  F no higher than this
                        if F >= abs(MF)
                            iF = (F+1)/2;
                            xx = ...
                            overlap(iR1, iR2, iK, s, ...
                                     N, MN, L, ML, f, mf, mR2, mK, ...
                                     J, iR, I, P, F, MF );
                            %p_exit = append_exit(p_exit,F_complex,iexit,F,xx^2);
                            p_exit(iexit,iF) = p_exit(iexit,iF) + xx^2;
                        end
                    end
                end
            end
        end
    end      
end


