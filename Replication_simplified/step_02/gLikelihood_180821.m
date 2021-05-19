function gLikelihood_180821(Theta,output_type,enableTest)

% Calculate optimal choice probabilities & joint log likelihood
tic
global beta lambda maxn t State Entry Exit Invest pr_down Merge T T2025 S Pi Kappa_c Kappa_e iterMLE LL ll...
    statespace newstate_e newstate_x newstate_i newstate_d newstate_m m_level b_level...
     V EV Policy W Lambda Weight Vpe EVpe Policy_pe sigma;

% Display candidate parameter values
kappa0_x = 0;           % Exit cost := 0
kappa0_c = Theta(1);    % Fixed cost of operation (common part)
kappa0_i = Theta(2);    % Innovation cost to be estimated
kappa0_i4 = Theta(3);   % Level-4 (frontier) firm's innovation cost [same with other, for now]
kappa0_m = Theta(4);    % Merger cost to be estimated
kappa0_e = Theta(5);    % Entry cost to be estimated
delta = 0;
fprintf('\nMLE iter #%4.0f: trying (kappa0_x, kappa0_c, kappa0_i, kappa0_i4, kappa0_m, kappa0_e, sigma, delta, lambda0)\n',iterMLE)
fprintf('                     = (  %4.4f,   %4.4f,   %4.4f,    %4.4f,   %4.4f,   %4.4f,  %4.2f,  %4.2f,    %4.2f)\n',...
                                kappa0_x, kappa0_c, kappa0_i, kappa0_i4, kappa0_m, kappa0_e, sigma, delta, lambda)

% Poisson distribution of stochastic synergy draws
synergy = zeros(4,1);                       % Prob. mass function with support {0,1,2,3}
synergy(1) = exp(-lambda);                  % k = 0 (no synergy)
synergy(2) = lambda * exp(-lambda);         % k = 1 (mild synergy)
synergy(3) = (lambda^2) * exp(-lambda) / 2; % k = 2 (strong synergy)
synergy(4) = 1 - sum(synergy(1:3, 1));      % k = 3 (very strong synergy)
% "(lambda^3)*exp(-lambda)/6" is correct but we aggregate all k >= 3 cases

%% Loops (over time & state) begin here

LambdaT = zeros(5,1,S);

% Solve the game from t = T2025 - 1: 
for t = (T2025 - 1):-1:1

    tic
    fprintf('Solving period %2.0f...    ',t)

    % Parameters to be estimated
    kappa_x = kappa0_x;                                 % Cost of exit
    kappa_c = Kappa_c + kappa0_c * ones(size(Kappa_c)); % Fixed cost of operation
    kappa_i = kappa0_i;                                 % Cost of innovation
    kappa_i4 = kappa0_i4;                               % Cost of innovation for LV-4 firms
    kappa_m = kappa0_m + delta * t;                     % Redundant but needed here for run
    kappa_e = kappa0_e * Kappa_e(t);                    % Cost of entry
    
    likelihood_180821(statespace...
        ,maxn, [], beta...
        ,kappa_x, kappa_c, kappa_i, kappa_e, kappa_m...
        ,t, delta, m_level, synergy...
        ,newstate_i, newstate_m, newstate_e, newstate_x...
        ,Pi, S, T2025, Lambda, V, W, EV, Policy, Policy_pe, EVpe, Vpe, sigma, b_level, kappa_i4);


    LambdaT(:,:,:) = Lambda(:,t,:);

    %% Exogenous Stochastic Depreciation
    stochastic_depreciation(statespace, pr_down, LambdaT, newstate_d, Lambda, S, t);

    if enableTest
        filename = sprintf('refdata/iter%dt%d', iterMLE, t);

        uut_output = { Lambda(:,t,:), V(:,t,:), W(:,t,:), EV(:,:,t+1,:), EVpe(:,t+1,:), Policy(:,:,t,:), Vpe(t,:), Policy_pe(:,t,:) };

        fprintf('comparing against refdata %s\n', filename)

        load(filename, 'output');

        compare_outputs(output, uut_output)
    end

    toc
end    % End of: for t = (T2025 - 1):-1:1

toc

end