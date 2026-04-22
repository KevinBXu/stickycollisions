%  brancing ratios to the different N's
%  late summer, 2023


% probability to get into F, summed over all

sumN0 = 0;
sumN1 = 0;
sumN2 = 0;   % strangely, not doubled here

for iexit = 1 : numexit
    N = exit(iexit,2);
    sumF = 0;
    for iF = 1 : length(F_complex)
        sumF = sumF + p_form(iF)*p_exit(iexit,iF);
    end
    if N == 0 
        sumN0 = sumN0 + sumF;
    end
    if N == 2 
        sumN1 = sumN1 + sumF;
    end
    if N == 4 
        sumN2 = sumN2 + sumF;
    end
end

sumall = sumN1 + sumN1 + sumN2;
disp(sumN0/sumall)
disp(sumN1/sumall)
disp(sumN2/sumall)

%  final branching ratios


%disp('        mf       N        ratio')
%disp('===============================')
%idecay = 0;
%for mfd = -fd : 2 : fd
%    for Nd = 0 : 2 : Nmax
%        idecay = idecay + 1;
%        disp([ mfd/2 Nd/2 BR(idecay) ])
%    end
%end


%figure(1)
%scatter(F_complex/2,p_form,'filled')
%xlabel('F')
%ylabel('probability to  form complex')
%hold off

