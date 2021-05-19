% MLE of a discrete-time model of entry, exit, innovation, merger, & "innovate-and-merge"
% With 32 frontier levels, 4 relative levels of firms, stochastic synergy & depreciation
% For twice-extended data (1996 Q1 - 2016 Q3) & monthly time interval

global alpha0 alpha1 alpha2 x mc Qj p...
    beta maxn lambda State Entry Exit Invest pr_down Merge T T2025 S Pi Kappa_c Kappa_e iterMLE LL ll...
    statespace newstate_e newstate_x newstate_i newstate_d newstate_m m_level b_level...
    V EV Policy W Lambda Weight Vpe EVpe Policy_pe EnableTest sigma;

diary diary_180822_sigma055.txt

%% Data Preparation
EnableTest = false;

% Load static parameters & basic data (including period profits)
load MakePi_EndoFrontier_180712.mat;   % Contains Pi (4 types x 83 quarters x 2240 states)
maxn = 5;
beta = .9 ^ (1/12); % Monthly discout factor
scaling = 1000;     % Scale-down kappa_c
sigma = 0.55;        % Logit scaling parameter
lambda = 1;         % Parameter of synergy distrib'n (Pareto)
Pi_qtr = Pi; 

% Fixed cost of operation (kappa_c) in data: SGA + CAPEX
% New estimate with 4 (own level) x 32 (frontier level):
Kappa_c_est = [0.00 	0.00 	0.00 	0.00 	0.00 	0.00 	0.00 	0.00 	0.00 	0.00 	0.00 	0.00 	1.70 	5.22 	8.74 	12.26 	15.78 	19.30 	22.82 	26.34 	29.86 	33.38 	36.90 	40.42 	43.94 	47.46 	50.98 	54.50 	58.03 	61.55 	65.07 	68.59; ...
               0.00 	0.00 	0.00 	0.00 	0.00 	0.00 	2.03 	5.55 	9.07 	12.59 	16.11 	19.63 	23.15 	26.67 	30.19 	33.71 	37.23 	40.75 	44.27 	47.79 	51.31 	54.83 	58.35 	61.87 	65.39 	68.91 	72.43 	75.95 	79.47 	82.99 	86.51 	90.03; ...
               2.35 	5.87 	9.39 	12.91 	16.43 	19.95 	23.47 	26.99 	30.51 	34.04 	37.56 	41.08 	44.60 	48.12 	51.64 	55.16 	58.68 	62.20 	65.72 	69.24 	72.76 	76.28 	79.80 	83.32 	86.84 	90.36 	93.88 	97.40 	100.92 	104.44 	107.96 	111.48; ...
               23.80 	27.32 	30.84 	34.36 	37.88 	41.40 	44.92 	48.44 	51.96 	55.48 	59.00 	62.52 	66.04 	69.56 	73.08 	76.60 	80.12 	83.64 	87.16 	90.69 	94.21 	97.73 	101.25 	104.77 	108.29 	111.81 	115.33 	118.85 	122.37 	125.89 	129.41 	132.93];
Kappa_c = (Kappa_c_est / 3) / scaling;  % Slice quarterly fixed cost into monthly numbers (and rescale)

% Convert quarterly data to monthly, by simple slicing
Pi = zeros(4,3*T,S);            % New, monthly profit matrix (Note: T is still quarterly at this point)
X = zeros(1,3*T);               % New, monthly PC shipments
P = zeros(3*T,S);               % New, monthly equilibrium prices
Q = zeros(3*T,S);               % New, monthly equilibrium outputs

Pi_mon = Pi_qtr ./ 3;           % Chopping quarterly profits into equal monthly slices
X_mon = X_data' ./ 3;           % Chopping quarterly PC shipments into equal monthly slices
P_mon = P1;                     % Don't chop prices
Q_mon = Q1 ./ 3;                % Chopping quarterly HDD shipments into equal monthly slices

for qtr = 1:T
    Pi(:, 3*(qtr-1)+1:3*(qtr-1)+3, :) = repmat(Pi_mon(:, qtr, :),[1 3 1]);
    X(1, 3*(qtr-1)+1:3*(qtr-1)+3) = repmat(X_mon(1, qtr),[1 3]);
    P(3*(qtr-1)+1:3*(qtr-1)+3, :) = repmat(P_mon(qtr, :),[3 1]);
    Q(3*(qtr-1)+1:3*(qtr-1)+3, :) = repmat(Q_mon(qtr, :),[3 1]);
end
T = 3 * T;                      % T was quarterly; now it's monthly

% Add a time trend to entry cost
Kappa_e = zeros(360,1);             % Scaling factor

% Number of active firms by type: (row 1 = frontier 1~32; row 2-5 = #firms of lv-1~4; column = time)
State = [1	1	1	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	3	3	3	3	3	3	4	4	4	4	4	4	5	5	5	5	5	5	5	5	5	6	6	6	7	7	7	7	7	7	8	8	8	8	8	8	9	9	9	9	9	9	9	9	9	9	9	9	9	9	9	10	10	10	10	10	10	10	10	10	11	11	11	11	11	11	11	11	11	11	11	11	12	12	12	12	12	12	12	12	12	12	12	12	13	13	13	13	13	13	13	13	13	13	13	13	13	13	13	13	13	13	13	13	13	13	13	14	14	14	14	14	14	14	14	14	14	14	14	14	14	14	14	15	15	15	15	15	15	15	15	15	15	15	15	16	16	16	16	16	16	16	16	16	16	16	17	17	17	17	17	17	17	17	17	17	17	17	17	17	17	17	17	17	17	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	18	19	19	19	19	19	19	19	19	19	19	19	19	19	19	19	19	19	19; ...
        6	5	4	4	4	4	3	3	3	2	2	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	1	1	1	2	3	4	3	2	2	2	1	1	1	1	1	1	1	1	1	1	1	1	0	0	0	0	0	0	0	0	0	0	0	0	1	1	1	1	1	1	1	1	1	0	0	0	1	1	1	1	1	1	2	2	3	4	4	4	4	4	4	4	4	4	4	4	4	4	4	3	2	2	2	2	2	2	1	1	1	1	1	2	3	4	3	2	2	2	2	2	2	1	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	1	1	2	3	4	3	3	3	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	0	0	0	0	0	0; ...
        4	5	6	6	6	6	6	6	5	6	6	7	8	8	7	6	6	6	5	5	4	3	3	3	3	3	3	3	3	3	2	2	2	1	1	1	1	2	2	3	3	3	3	4	4	4	3	3	3	3	3	3	3	3	4	5	5	4	4	4	3	3	3	4	4	4	5	5	5	5	5	5	4	4	4	4	5	4	3	1	2	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	4	4	4	4	4	4	4	4	4	4	4	4	3	3	3	4	4	4	3	3	3	4	4	5	4	4	4	4	4	5	3	3	2	1	1	1	1	1	1	1	1	1	1	1	1	1	1	2	3	3	3	3	3	2	3	3	3	3	4	3	2	0	1	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	3	3	2	2	2	3	2	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	1	1	1	1; ...
        1	1	1	1	1	1	1	0	1	1	1	1	1	1	2	3	3	3	3	3	4	5	5	5	5	5	5	5	4	3	3	3	2	3	3	3	3	3	3	2	3	4	5	4	3	2	3	3	4	4	4	4	5	5	4	3	4	5	5	5	6	6	6	3	3	3	2	2	2	2	2	2	2	2	2	2	1	1	1	2	2	2	2	2	2	2	2	2	2	2	3	3	3	3	3	2	2	2	2	3	3	3	3	3	3	3	3	3	3	3	3	2	2	2	3	3	3	3	3	2	2	2	2	2	2	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	0	0	0	0	0	0	0	1	1	1	1	1	0	0	0	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	0	0	0	1	2	1	1	1	1	1	1	1	1	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	1; ...
        2	2	2	2	2	2	2	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	3	4	5	5	5	6	6	6	6	6	6	6	6	5	4	3	3	4	5	5	5	4	4	4	4	3	3	3	3	2	1	1	1	1	1	1	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	1	1	1	1	1	2	2	2	2	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	1	0	0	0	0	0	0	0	0	1	1	1	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	1	1	1	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	1	1	1];

% We delieted some data on exit and merge to accomodate the simplified case
% Number of entry
Entry = zeros(1,T-1);
Entry(1, 36) = 1;   % Conner Technology (later ExcelStor) entry in Dec-98
    
% Number of exits by liquidation: (row = low & high; column = time period)
Exit = zeros(4,T-1);
Exit(2, 18) = 1;    % Micropolis in Jun-97
Exit(2, 57) = 1;    % NEC in Sep-00

% Number of innovations
Invest = [1	1	0	0	0	0	0	0	1	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	0	0	0	0	0	1	0	0	0	0	0	0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0; ...
        0	0	0	0	0	0	0	1	0	0	0	0	0	1	1	0	0	0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0; ...
        0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	0	0	1	0	0	0	0	0	0	0	0	0	0	0	1	1	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0; ...
        0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0];

% Probability of stochastic depreciation: 
% Levels went down 53 out of (1366 - 5 - 6) times-at-risk (# firm-years minus exits, excluding level-1 firms)
pr_down = zeros(4,1);
pr_down(1) = 0; 
pr_down(2) = 16 / (614 - 3 - 1);    % = 2.62% [= #down / (#active - #exit - #target)] 
pr_down(3) = 13 / (411 - 0 - 2);    % = 3.18%
pr_down(4) = 14 / (509 - 0 - 0);    % = 2.75%

% Number of mergers: (row 1 = level-1 buys level-1, row 2 = level-1 buys level-2, ..., row 16 = level-4 buys level-4; column = time period)
Merge = zeros(4*4,T-1);
Merge(4 * (3 - 1) + 3, 63) = 1; % Maxtor (level 3) buys Quantum (3) in Mar-01 
                                % ...and became level 4 (mild synergy from high-high)
Merge(4 * (4 - 1) + 2, 126) = 1; % Seagate (level 4) buys Maxtor (2) in May-06 
                                % ...and remained level 4 (no synergy from frontier-mid)
Merge(4 * (1 - 1) + 1, 195) = 1; % WDC (level 1) buys Hitachi (1) in Mar-12 
                                % ...and became level 4 (very^2 strong synergy from low-low)

% Prepare weighting for likelihood calculation
ActionFlag = zeros(5,T-1);                                          % Flag for types & time w/ observed actions
for type = 1:4                                                      % Flag for lv-1~4's actions
    ActionFlag(type,:) = Exit(type,:) + Invest(type,:) + sum(Merge(4*(type-1)+1:4*type,:));
end
ActionFlag(5,:) = Entry(1,:);                                       % Potential entrants as "type-5" player
Weight = [State(2:5,1:T-1); ones(1,T-1)] ./ repmat(maxn, 5, T-1);   % Weight by # of firms (incl. "nobody")
Weight(:,sum(ActionFlag) == 1) = ActionFlag(:,sum(ActionFlag) == 1);% Replace weight by 0 & 1 if action observed

% Prepare state-index matrices for transitions due to exit & innoavtion
load newstate_x.mat;        % Index of new state = f(current state, exiter's level)
load newstate_i.mat;        % Index of new state = f(current state, innovator's level)
load newstate_d.mat;        % Index of new state = f(current state, depreciator's level)
load newstate_e.mat;        % Index of new state = f(current state)

% Prepare state-index matrix for transitions due to merger
load newstate_m.mat;        % Index of new state = f(current state, levels of acquirer, target, & merged entity)

% Merged entity's level(s) as a function of the acquiring & the target firms
load m_level.mat;           % Merged entity's level can take (up to) 4 values depending on synergy realizations
                            % [a_level x t_level x 4 synergy cases]
load b_level.mat;

if false
    % This crops the matrices down to a size that can be run without lots
    % of memory, however the results are not going to be correct in anyway.
    S=2240;
    T=100;

    newstate_x = newstate_x(1:S,:);
    newstate_i = newstate_i(1:S,:);
    newstate_d = newstate_d(1:S,:);
    newstate_e = newstate_e(1:S,:);
    newstate_m = newstate_m(1:S,:,:,:);
    Pi = Pi(:,1:T,1:S);
    statespace = statespace(1:S,:);
    newstate_x = floor(mod(newstate_x,S))+1;
    newstate_i = floor(mod(newstate_i,S))+1;
    newstate_d = floor(mod(newstate_d,S))+1;
    newstate_e = floor(mod(newstate_e,S))+1;
    newstate_m = floor(mod(newstate_m,S))+1;
end

% Longer horizon: 249 (Jan1996-Sep2016) + 3 (Oct-Dec2016) + 12 months x 9 years (2017-2025) = 360 months
T2025 = T + 111;            

% Extend Pi until Dec-2025
Pi2016 = Pi;
Pi2025 = zeros(4, T2025, S);        % Jan1996 - Dec2025 (360 months)
Pi2025(:, 1:T, :) = Pi2016;         % Jan1996 - Sep2016 (249 months)

LinearDecline = zeros(111,1);
for t = 1:111
    LinearDecline(t) = 1 - (1/111) * t;
    Pi2025(:, (T + t), :) = Pi2025(:, T, :) .* LinearDecline(t);
    %Kappa_c2025(:, (T + t)) = Kappa_c2025(:, T) .* LinearDecline(t);
end
Pi = Pi2025;
clear Pi2025;% Kappa_c2025;

% Initialize solution variables
V = zeros(4,T2025,S);           % Ex-ante value functions
V(:,T2025,:) = zeros(4,1,S);    % True end of the HDD world
EV = zeros(4,11,T2025,S);       % Alternative-specific value functions 
Policy = zeros(4,11,T2025-1,S); % Optimal choice probabilities
W = zeros(5,T2025-1,S);         % Value function when it's not the firm's turn
Lambda = zeros(5,T2025,S);      % Umbrella value function

% Potential entrant's value & policy
Vpe = zeros(T2025,S);
EVpe = zeros(2,T2025,S);        % EV of {stay-out & enter}
Policy_pe = zeros(2,T2025-1,S); % CCP of {stay-out, enter}

%% Estimation

% Initialize MLE
iterMLE = 1;            % Iteration count
x0 = [0.0111553904807242, 0.503647209700581, 0.938950879578958, 1.36755860414072, 0.151411824768284];  % 18820ci4 estimates (sigma = 0.6)
options = optimset('Display','iter','TolFun',1e-7,'TolX',1e-6,'MaxIter',1000,'MaxFunEvals',2500);

% MLE
Theta = fminsearch(@Likelihood_180821,x0,options,1);

% Display estimates
fprintf('\n -------------------------------------------------------------------------------');
fprintf('\n Log likelihood:         %4.4f',ll);
fprintf('\n Parameter:       (kappa_x   kappa_c)  kappa_i   kappa_i4   kappa_m   kappa_e   (sigma   delta   lambda)');
fprintf('\n Coefficient:      0.0000    %4.4f    %4.4f    %4.4f    %4.4f    %4.4f    %4.4f    see iter    1.0', Theta, sigma);
fprintf('\n -------------------------------------------------------------------------------\n');

save Theta_180822_sigma055.mat Theta;
diary off;