function f = Likelihood_180821(Theta,output_type)
global beta lambda maxn t State Entry Exit Invest pr_down Merge T T2025 S Pi Kappa_c Kappa_e iterMLE LL ll...
    statespace newstate_e newstate_x newstate_i newstate_d newstate_m m_level b_level...
     V EV Policy W Lambda Weight Vpe EVpe Policy_pe runType;

    % run type
    %   0 - normal run using accelerated model
    %   1 - normal run using original model
    %   2 - normal run using full mex model
    %   3 - generate reference data using original model
    %   4 - compare accelerated model with reference data
    %   5 - compare full mex model with reference data
    runType = 2

    % run accelerated model
    if runType == 0 || runType == 4
        fLikelihood_7type(Theta, output_type, runType==4);
    end
    
    % run original model
    if runType == 1 || runType == 3
        oLikelihood(Theta, output_type, runType==3);
    end
    
    % run superfast model
    if runType == 2 || runType == 5
        gLikelihood_180821(Theta, output_type, runType==5);
    end


    %% Likelihood of observing the actual choices in the data

    LL = zeros(5,T-1);
    for t = 1:(T-1)

        % Numbers of firms that took particular actions
        fron = State(1, t);     % frontier level
        N = State(2:5, t);      % [4 x 1]: Number of firms by type
        E = Entry(1, t);        % [1 x 1]: Indicator of entry
        X = Exit(:, t);         % [4 x 1]: Indicator of exit (by exiter type)
        I = Invest(:, t);       % [4 x 1]: Indicator of innovation/investment (by innovator type)
        M1 = Merge(1:4, t);     % [4 x 1]: Indicator of merger, with level-1 as Acquirer (by Target type)
        M2 = Merge(5:8, t);     % [4 x 1]: Indicator of merger, with level-2 as Acquirer (by Target type)
        M3 = Merge(9:12, t);    % [4 x 1]: Indicator of merger, with level-3 as Acquirer (by Target type)
        M4 = Merge(13:16, t);   % [4 x 1]: Indicator of merger, with level-4 as Acquirer (by Target type)
        
        % Choice probabilities in the actual state at t (in data)
        astate = (statespace(:,1) == fron) .* (statespace(:,2) == N(1)) ...
              .* (statespace(:,3) == N(2)) .* (statespace(:,4) == N(3)) ...
              .* (statespace(:,5) == N(4)) == 1;    % Index of the Actual state at t in Data
        % [1 x 4]: Level-1, 2, 3, & 4's CCPs {exit, idle, invest, merge}, respectively
        p1 = Policy(1, 1:7, t, astate);   
        p2 = Policy(2, 1:7, t, astate);   
        p3 = Policy(3, 1:7, t, astate);   
        p4 = Policy(4, 1:7, t, astate);   
        p5 = Policy_pe(:, t, astate);   % [1 x 2]: Pot. ent.'s CCP {stay-out, enter}
        fprintf('Period %2.0f\n',t)
        fprintf('   Lv-1 CCP = (%0.2f,%0.2f,%0.2f,%0.4f,%0.4f,%0.4f,%0.4f)\n',p1(1),p1(2),p1(3),p1(4),p1(5),p1(6),p1(7))
        fprintf('   Lv-2 CCP = (%0.2f,%0.2f,%0.2f,%0.4f,%0.4f,%0.4f,%0.4f)\n',p2(1),p2(2),p2(3),p2(4),p2(5),p2(6),p2(7))
        fprintf('   Lv-3 CCP = (%0.2f,%0.2f,%0.2f,%0.4f,%0.4f,%0.4f,%0.4f)\n',p3(1),p3(2),p3(3),p3(4),p3(5),p3(6),p3(7))
        fprintf('   Lv-4 CCP = (%0.2f,%0.2f,%0.2f,%0.4f,%0.4f,%0.4f,%0.4f)\n',p4(1),p4(2),p4(3),p4(4),p4(5),p4(6),p4(7))
        fprintf('   Potential entrant CCP: Pr(out, enter) = (%0.4f, %0.4f)\n', p5(1), p5(2));

        % Fool-proof: Avoid "log(0) = -Inf" by replacing zero-CCPs with very small numbers
        p1(p1 < 1e-8) = 1e-8;
        p2(p2 < 1e-8) = 1e-8;
        p3(p3 < 1e-8) = 1e-8;
        p4(p4 < 1e-8) = 1e-8;
        p5(p5 < 1e-8) = 1e-8;

        % Log likelihood
        LL(1,t) = Weight(1,t) * (X(1) * log(p1(1)) + I(1) * log(p1(3)) + M1(1) * log(p1(4)) ...
                                + M1(2) * log(p1(5)) + M1(3) * log(p1(6)) + M1(4) * log(p1(7)) ...
                                + (1 - X(1) - I(1) - sum(M1)) * log(p1(2)));
        LL(2,t) = Weight(2,t) * (X(2) * log(p2(1)) + I(2) * log(p2(3)) + M2(1) * log(p2(4)) ...
                                + M2(2) * log(p2(5)) + M2(3) * log(p2(6)) + M2(4) * log(p2(7)) ...
                                + (1 - X(2) - I(2) - sum(M2)) * log(p2(2)));
        LL(3,t) = Weight(3,t) * (X(3) * log(p3(1)) + I(3) * log(p3(3)) + M3(1) * log(p3(4)) ...
                                + M3(2) * log(p3(5)) + M3(3) * log(p3(6)) + M3(4) * log(p3(7)) ...
                                + (1 - X(3) - I(3) - sum(M3)) * log(p3(2)));
        LL(4,t) = Weight(4,t) * (X(4) * log(p4(1)) + I(4) * log(p4(3)) + M4(1) * log(p4(4)) ...                                
                                + M4(2) * log(p4(5)) + M4(3) * log(p4(6)) + M4(4) * log(p4(7)) ...
                                + (1 - X(4) - I(4) - sum(M4)) * log(p4(2)));
        LL(5,t) = Weight(5,t) * (E * log(p5(2)) + (1 - E) * log(p5(1)));
        fprintf('    Log likelihood = %10.4f \n',sum(LL(:,t)))
    end

    % Prepare output of this function file
    ll = sum(sum(LL));
    if output_type == 1; 
        f = -ll;    % For the estimation (to maximize likelihood)
    elseif output_type == 2; 
        f = ll;     % For calculating standard errors
    end

    % Count the number of iterations
    iterMLE = iterMLE + 1;

end
