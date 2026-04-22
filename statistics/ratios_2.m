%  brancing ratios to the different N's
%  late summer, 2023


% probability to get into F, summed over all

% output states in terms of [mF mK mRb mN mL] (doubled!)



state1 = [2 -8 -1 0 0];
state2 = [2 -6 -3 0 0];
state3 = [2 -8 -3 0 2];
state4 = [2 -8 1 0 -2];
state5 = [2 -8 3 0 -4];
fstates = [state1; state2; state3; state4; state5];

state1 = [2 -8 -3 0 -2];
state2 = [2 -8 -1 0 -4];
fstates = [state1; state2];

sums = zeros(size(fstates,1),1);

for iexit = 1 : numexit
    sumF = 0;
    for iF = 1 : length(F_complex)
        sumF = sumF + p_form(iF)*p_exit(iexit,iF);
    end
    
    % exit state in terms of [mF mK mRb mN mL] (doubled!)
    q_nums = [exit(iexit,7) exit(iexit,9) exit(iexit,8) exit(iexit,3) exit(iexit,5)];
    
    for istate = 1 : size(fstates,1)
        fstate = fstates(istate,:);
        matched = true;
        for i = 1 : 5
            if not(fstate(i) == q_nums(i)) 
                matched = false;
            end
        end
        if matched
            disp([istate sumF])
            sums(istate) = sums(istate) + sumF;
        end
    end

end

sumall = sum(sums);
disp("Ratio")
for istate = 1 : size(fstates,1)
    disp(sums(istate) / sumall)
end