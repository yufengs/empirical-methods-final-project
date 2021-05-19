function f = foc(q)
% Firm's first-order condition for profit maximization, with log-linear demand function

% Parameters
global alpha0 alpha1 alpha2 alpha3 detrend_P0 detrend_P1 detrend_Q0 detrend_Q1 t x thai mc Qj;

% Price
Q = max(Qj + q, 0.001);     % max(., .) operation to avoid negative Q while searching for optimal q*
logQ_detrended = log(Q) - (detrend_Q0 + detrend_Q1 * 42);    % detrend Q
logP_detrended = (1 / alpha1) * logQ_detrended - (alpha2 / alpha1) * log(x) - (alpha0 / alpha1) - (alpha3 / alpha1) * thai;
logP = logP_detrended + (detrend_P0 + detrend_P1 * 42);      % retrend P
P = exp(logP);

% Slope of demand
dPdQ = (alpha1 * Q / P) ^ (-1);

% FOC: LHS - RHS = 0
f = P + dPdQ * q - mc;

end

