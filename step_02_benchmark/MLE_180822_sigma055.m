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
load MakePi_EndoFrontier_180712.mat;   % Contains Pi (4 types x 83 quarters x 76160 states)
maxn = 14;
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

% Number of entry
Entry = zeros(1,T-1);
Entry(1, 36) = 1;   % Conner Technology (later ExcelStor) entry in Dec-98
    
% Number of exits by liquidation: (row = low & high; column = time period)
Exit = zeros(4,T-1);
Exit(1, 6) = 1;     % HP in Jun-96
Exit(2, 18) = 1;    % Micropolis in Jun-97
Exit(2, 30) = 1;    % JTS in Jun-98
Exit(2, 57) = 1;    % NEC in Sep-00
Exit(1, 168) = 1;   % ExcelStor in Dec-09

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
Merge(4 * (3 - 1) + 3, 63) = 1; % Merger #1: Maxtor (level 3) buys Quantum (3) in Mar-01 
                                % ...and became level 4 (mild synergy from high-high)
Merge(4 * (1 - 1) + 3, 84) = 1; % Merger #2: Hitachi (level 1) buys IBM (3) in Dec-02 
                                % ...and became level 3 (no synergy from low-high)
Merge(4 * (4 - 1) + 2, 126) = 1; % Merger #3: Seagate (level 4) buys Maxtor (2) in May-06 
                                % ...and remained level 4 (no synergy from frontier-mid)
Merge(4 * (2 - 1) + 1, 164) = 1; % Merger #4: Toshiba (level 2) buys Fujitsu (1) in Sep-09 
                                % ...and remained level 2 (no synergy from low-low)
Merge(4 * (3 - 1) + 1, 192) = 1; % Merger #5: Seagate (level 3) buys Samsung (1) in Dec-11 
                                % ...and remained level 3 (no synergy from high-low)
Merge(4 * (1 - 1) + 1, 195) = 1; % Merger #6: WDC (level 1) buys Hitachi (1) in Mar-12 
                                % ...and became level 4 (very^2 strong synergy from low-low)
                                
% "Coarse" versions Merger data, for better fit
Merge2 = zeros(4, T-1);         % By acquirer level (no distinction of target levels)
Merge2(3, 63) = 1;
Merge2(1, 84) = 1;
Merge2(4, 126) = 1;
Merge2(2, 164) = 1;
Merge2(3, 192) = 1;
Merge2(1, 195) = 1;

Merge3 = zeros(1, T-1);         % No distinction of acquirer- or target-levels
Merge3(63) = 1;
Merge3(84) = 1;
Merge3(126) = 1;
Merge3(164) = 1;
Merge3(192) = 1;
Merge3(195) = 1;

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
    S=5000;
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
%Likelihood2(x0,0);     % evaluate LL only once (for no-estimation run)
%Theta = x0;            % evaluate LL only once (for no-estimation run)
Theta = fminsearch(@Likelihood_180821,x0,options,1);

% Display estimates
fprintf('\n -------------------------------------------------------------------------------');
fprintf('\n Log likelihood:         %4.4f',ll);
fprintf('\n Parameter:       (kappa_x   kappa_c)  kappa_i   kappa_i4   kappa_m   kappa_e   (sigma   delta   lambda)');
fprintf('\n Coefficient:      0.0000    %4.4f    %4.4f    %4.4f    %4.4f    %4.4f    %4.4f    see iter    1.0', Theta, sigma);
fprintf('\n -------------------------------------------------------------------------------\n');

save Theta_180822_sigma055.mat Theta;
diary off;

%% Post-Estimation Exercises (x3):
% 1. Show the values and policies along the actual history in Data
% 2. Simulate 10,000 histories based on the ML estimates
% 3. Calculate welfare outcomes 

%% Task 1. Show the values and policies along the actual history in Data

V_actual = zeros(4,T);
Policy_actual = zeros(4,11,T-1);
V_actual_pe = zeros(T,1);
Policy_actual_pe = zeros(2,T-1);

for t = 1:T
    Fron = State(1,t);          % Frontier level
    N = State(2:5,t);           % [4 x 1]: Number of firms by type
    astate = find((statespace(:,1) == Fron) .* (statespace(:,2) == N(1)) .* (statespace(:,3) == N(2)) ...
               .* (statespace(:,4) == N(3)) .* (statespace(:,5) == N(4)) == 1); % State index
           
    V_actual(:,t) = V(:,t,astate);
    V_actual_pe(t) = Vpe(t,astate);
        
    if t < T
        Policy_actual(:,:,t) = Policy(:,:,t,astate);
        Policy_actual_pe(:,t) = Policy_pe(:,t,astate);
    end
end

V_actual = V_actual';

%% Task 2. Simulate 10,000 histories, based on the ML estimates

% Initialize
NS = 10000;                     % Number of simulations
State_simu = single(zeros(5,T2025,NS));     % Simulated histories
Exit_simu = single(zeros(4,T2025-1,NS));    % Simulated number of exits by type
Invest_simu = single(zeros(4,T2025-1,NS));  % Simulated number of innovations by type
Merge_simu = single(zeros(4*4,T2025-1,NS)); % Simulated number of mergers by type-buys-type
Entry_simu = single(zeros(T2025-1,NS));     % Simulated number of entry

seed = 645456;                  % Seed for random draws
rng(seed);                      % Set seed
U = single(rand(4,T2025-1,NS));             % Draws from U(0,1) to simulate Nature's choice, Player's choice, Synergy, & Depreciation

turn = single(zeros(T2025-1,NS));           % Record Nature's choices of (the types of) players: (t,ns)
target = single(zeros(T2025-1,NS));         % Record merger Target type
synergy = single(zeros(T2025-1,NS));        % Record Synergy realizations
depreciation = single(zeros(T2025-1,NS));   % Record Depreciation shock (which type got hit by it)

fprintf('\n Simulating 10,000 histories...');
for ns = 1:NS        % Loop 1: for each simulated history...
    
    fprintf('\n Simulating history #%5.0f',ns);
    State_simu(:,1,ns) = State(:,1);                            % State in period 1 is given by data
    
    for t = 1:(T2025-1)  % Loop 2: for each period...
    
        % This period's state has already been determined
        fron = State_simu(1, t, ns);
        n1 = State_simu(2, t, ns);
        n2 = State_simu(3, t, ns);
        n3 = State_simu(4, t, ns);
        n4 = State_simu(5, t, ns);
        n = sum(State_simu(2:5,t,ns),1);    % Number of all active firms
        n5 = 1;                             % One potential entrant always exists
        state = find((statespace(:,1) == fron) .* (statespace(:,2) == n1) .* (statespace(:,3) == n2) ...
                  .* (statespace(:,4) == n3) .* (statespace(:,5) == n4) == 1);
              
        % Simulation draws
        dn = U(1,t,ns);                                         % Random draw for Nature's choice
        dp = U(2,t,ns);                                         % Random draw for Player's choice
        ds = U(3,t,ns);                                         % Random draw for Synergy realization
        dd = U(4,t,ns);                                         % Random draw for Depreciation target
            
        % Given (n1,n2,...,n4,n5), Nature chooses a (type of) player, then the chosen type makes a choice
        
        % Nature decides which (type of) player moves
        turn(t, ns) = (dn < (n1 / maxn)) * 1 ...
                    + ((dn >= (n1 / maxn)) && (dn < ((n1 + n2) / maxn))) * 2 ...
                    + ((dn >= ((n1 + n2) / maxn)) && (dn < ((n1 + n2 + n3) / maxn))) * 3 ...
                    + ((dn >= ((n1 + n2 + n3) / maxn)) && (dn < ((n1 + n2 + n3 + n4) / maxn))) * 4 ...
                    + ((dn >= ((n1 + n2 + n3 + n4) / maxn)) ...
                        && (dn < ((n1 + n2 + n3 + n4 + n5) / maxn))) * 5;
        type = turn(t, ns);                                 % Decision-maker's type
        if type > 0 && type < 5                             % If incumbent's turn
            px = Policy(type, 1, t, state);                 % Exit probability
            pc = Policy(type, 2, t, state);                 % Stay-alone probability
            pi = Policy(type, 3, t, state);                 % Innovation probability
            pm1 = Policy(type, 4, t, state);                % Merger (with level-1 rival) probability
            pm2 = Policy(type, 5, t, state);                % Merger (with level-2 rival) probability
            pm3 = Policy(type, 6, t, state);                % Merger (with level-3 rival) probability
            pm4 = Policy(type, 7, t, state);                % Merger (with level-4 rival) probability
            pother = (px + pc + pi + pm1 + pm2 + pm3 + pm4);% Shorthand for the first 7 ccp's
            pmi1 = Policy(type, 8, t, state);               % Innovate & merge with lv-1
            pmi2 = Policy(type, 9, t, state);               % Innovate & merge with lv-2
            pmi3 = Policy(type, 10, t, state);              % Innovate & merge with lv-3
            pmi4 = Policy(type, 11, t, state);              % Innovate & merge with lv-4
            ps1 = exp(-lambda);                             % Probability that there is no synergy
            ps2 = lambda * exp(-lambda);                    % Probability that there is mild synergy
            ps3 = (lambda ^ 2) * exp(-lambda) / 2;          % Probability that there is strong synergy

            % Simulated choices (exit & innovation)
            Exit_simu(type, t, ns) = dp < px;                                   % = 1 if exit (0 o/w)
            Invest_simu(type, t, ns) = (dp > (px + pc) && dp < (px + pc + pi)) ...
                                        || (dp > pother);                       % = 1 if innov (0 o/w)

            % Simulated choices (merger): Target type = 1 thru 4 (remains 0 if no-merger)
            target(t, ns) = (dp > (px + pc + pi) && dp < (px + pc + pi + pm1)) * 1 ...
                          + (dp > (px + pc + pi + pm1) && dp < (px + pc + pi + pm1 + pm2)) * 2 ...
                          + (dp > (px + pc + pi + pm1 + pm2) ...
                            && dp < (px + pc + pi + pm1 + pm2 + pm3)) * 3 ...
                          + (dp > (px + pc + pi + pm1 + pm2 + pm3) ...
                            && dp < (px + pc + pi + pm1 + pm2 + pm3 + pm4)) * 4 ...
                          + (dp > pother && dp < (pother + pmi1)) * 1 ...
                          + (dp > (pother + pmi1) && dp < (pother + pmi1 + pmi2)) * 2 ...
                          + (dp > (pother + pmi1 + pmi2) ...
                            && dp < (pother + pmi1 + pmi2 + pmi3)) *3 ...
                          + (dp > (pother + pmi1 + pmi2 + pmi3) ...
                            && dp < (pother + pmi1 + pmi2 + pmi3 + pmi4)) * 4;

            % Synergy realization = {1, 2, 3, 4} (rationalization, mild/strong/very strong synergy)
            synergy(t, ns) = (ds < (ps1)) * 1 ...
                           + (ds > (ps1) && ds < (ps1 + ps2)) * 2 ...
                           + (ds > (ps1 + ps2) && ds < (ps1 + ps2 + ps3)) * 3 ...
                           + (ds > (ps1 + ps2 + ps3)) * 4;     
        elseif type == 5                                    % If potential entrant's turn
            pe = Policy_pe(2,t,state);                      % Entry probability
            Entry_simu(t, ns) = dp < pe;                    % = 1 if enter (0 otherwise)
        end

        % State transition  
        if type > 0 && type < 5
            
            if target(t,ns) > 0                                 % If merger                
                a_level = type;                                 % Acquirer's level (just re-labeling)
                t_level = target(t, ns);                        % Target's level (just re-labeling)
                s_case = synergy(t, ns);                        % Merger synergy outcome (4 cases)
                Merge_simu(4 * (a_level - 1) + t_level, t, ns) = 1;   % Record this merger
                % Find the resulting state due to merger
                if dp < pother
                    state_prime = newstate_m(state, a_level, t_level, m_level(a_level, t_level, s_case));
                else
                    state_prime = newstate_m(state, a_level, t_level, b_level(a_level, t_level, s_case));
                end
            elseif Exit_simu(type, t, ns) > 0                   % If exit                
                state_prime = newstate_x(state, type);          % Find the resulting state due to exit
            elseif Invest_simu(type, t, ns) > 0                 % If invest
                state_prime = newstate_i(state, type);           
            else                                                % If stay-alone or stay-out
                state_prime = state;                            % State does not change
            end
            
        elseif type == 5

            if Entry_simu(t, ns) > 0                            % If enter
                state_prime = newstate_e(state);
            else                                                % If stay-alone or stay-out
                state_prime = state;                            % State does not change
            end
            
        elseif type == 0
            
            state_prime = state;                                % When nobody moves, nothing happens
            
        end

        % Exogenous stochastic depreciation (record which type got hit)
        depreciation(t, ns) = (dd < n1 * pr_down(1)) * 1 ...
            + (dd > n1 * pr_down(1) && dd < (n1 * pr_down(1) + n2 * pr_down(2))) * 2 ...
            + (dd > (n1 * pr_down(1) + n2 * pr_down(2)) && ...
               dd <  (n1 * pr_down(1) + n2 * pr_down(2) + n3 * pr_down(3))) * 3 ...
            + (dd > (n1 * pr_down(1) + n2 * pr_down(2) + n3 * pr_down(3)) && ...
               dd < (n1 * pr_down(1) + n2 * pr_down(2) + n3 * pr_down(3) + n4 * pr_down(4))) * 4 ...
            + (dd > (n1 * pr_down(1) + n2 * pr_down(2) + n3 * pr_down(3) + n4 * pr_down(4))) * 0;
        if depreciation(t, ns) > 0                          % If someone got hit
            dtype = depreciation(t, ns);
            state_prime = newstate_d(state_prime, dtype);   % State changes
        end
        
        % Exogenous advance of technological frontier in certain periods
        if t == 18 || t == 24 || t == 30 || t == 39 || t == 42 || t == 48 || t == 69 || t == 78 || t == 102 || t == 141 || t == 183 || t == 231
            state_prime = newstate_i(state_prime, 4);
        end

        State_simu(:, t+1, ns) = statespace(state_prime, :)';   % Record next period's state
        
    end        % End of loop: for t = 1:(T2025-1)
    
end        % End of loop: for ns = 1:NS

fprintf('\n ...done. \n');

% Summarize 10,000 simulations into average state evolution and # of actions
State_avg = sum(State_simu, 3) / NS;                            % Average over NS simulated histories
State_avg_rd = round(State_avg);                                % Rounding frontier & firm counts
State_avg_prime = State_avg';                                   % For easy graph-making
State_avg_sum = sum(State_avg_prime(:, 2:5), 2);                % Total number of firms
State_ = State';                                        
Fron_fit = [State_(:,1), State_avg_prime(1:249, 1)];            % Frontier level: data vs. model
Exit_avg = sum(Exit_simu, 3) / NS;
Exit_avg_rd = round(Exit_avg);
Exit_avg_prime = Exit_avg';
Invest_avg = sum(Invest_simu, 3) / NS;
Invest_avg_rd = round(Invest_avg);
Invest_avg_prime = Invest_avg';
Merge_avg = sum(Merge_simu, 3) / NS;
Merge_avg_rd = round(Merge_avg);
Merge_avg_prime = Merge_avg';
Entry_avg = sum(Entry_simu, 2) / NS;
Entry_avg_rd = round(Entry_avg);

% Summarize some key "moments"
Exit_avg_sum_t248 = sum(sum(Exit_avg_prime(1:248,:)));
Invest_avg_sum_t248 = sum(sum(Invest_avg_prime(1:248,:)));
Merge_avg_sum_t248 = sum(sum(Merge_avg_prime(1:248,:)));
Entry_avg_sum_t248 = sum(Entry_avg(1:248,:));
State_avg_t249 = State_avg_prime(249,:);
State_avg_sum_t249 = State_avg_sum(249);

Exit_avg_sum_t359 = sum(sum(Exit_avg_prime));
Invest_avg_sum_t359 = sum(sum(Invest_avg_prime));
Merge_avg_sum_t359 = sum(sum(Merge_avg_prime));
Entry_avg_sum_t359 = sum(Entry_avg);
State_avg_t360 = State_avg_prime(360,:);
State_avg_sum_t360 = State_avg_sum(360);

%% Task 3. Calculate social-welfare outcomes

% Load "master" CS & PS matrices [T2025 x S]
load MakeCSPS_T2025_180713.mat;                  
clear CS;
load CS_cutoff500.mat;
load CS_cutoff1000.mat;

% Pick up CS & PS along the simulated history
State_simu_id = zeros(T2025,NS);        % State-ID# version of 10,000 simulated histories
P_simu = zeros(T,NS);
CS_cutoff500_simu = zeros(T2025, NS);
CS_cutoff1000_simu = zeros(T2025, NS);
PS_simu = zeros(T2025, NS);

for ns = 1:NS
    
    fprintf('\n Calculating (CS, PS, FC) in simulation #%5.0f (of #%5.0f)', ns, NS)
    
    for t = 1:T2025

        Fron = State_simu(1, t, ns);    % Frontier level
        N = State_simu(2:5,t,ns);       % [4 x 1]: # firms by type
        n = sum(N);                     % # of all firms
        
        % index of state (Frontier, n1,...,n4)
        state = find((statespace(:,1) == Fron) .* (statespace(:,2) == N(1)) .* (statespace(:,3) == N(2)) ...
                  .* (statespace(:,4) == N(3)) .* (statespace(:,5) == N(4)) == 1);
              
        State_simu_id(t,ns) = state;
        if t <= T
            P_simu(t, ns) = P(t, state);    % I didn't calculate P for t > 234
        end
        
        CS_cutoff500_simu(t, ns) = CS_cutoff500(t, state);
        CS_cutoff1000_simu(t, ns) = CS_cutoff1000(t, state);
        %CS_simu(t, ns) = CS(t, state);
        PS_simu(t, ns) = PS(t, state);

    end
end

% Average across 10,000 simulations
P_avg = sum(P_simu, 2) ./ NS;
CS_cutoff500_avg = sum(CS_cutoff500_simu, 2) / NS;
CS_cutoff1000_avg = sum(CS_cutoff1000_simu, 2) / NS;
PS_avg = sum(PS_simu, 2) / NS;
SW_cutoff500 = [CS_cutoff500_avg, PS_avg, CS_cutoff500_avg + PS_avg];
SW_cutoff1000 = [CS_cutoff1000_avg, PS_avg, CS_cutoff1000_avg + PS_avg];

% Prepare discount factor
beta = .9 ^ (1/12);     % Monthly
Beta = zeros(360,1);    % Vector for all periods
for t = 1:360
    Beta(t,1) = beta ^ (t-1);
end

% Using CS with '$500/GB-cutoff'

% Welfare PDVs for all periods (1996-2025)
SWpv_cutoff500 = Beta' * SW_cutoff500;

% Welfare PDVs for 1996-2010
Beta_1996_2010 = Beta(1:180, 1);
SW_cutoff500_1996_2010 = SW_cutoff500(1:180, :);
SWpv_cutoff500_1996_2010 = Beta_1996_2010' * SW_cutoff500_1996_2010;

% Welfare PDVs for 2011-2025
Beta_2011_2025 = Beta(181:360, 1);
SW_cutoff500_2011_2025 = SW_cutoff500(181:360, :);
SWpv_cutoff500_2011_2025 = Beta_2011_2025' * SW_cutoff500_2011_2025;

% Using CS with '$1000/GB-cutoff'

% Welfare PDVs for all periods (1996-2025)
SWpv_cutoff1000 = Beta' * SW_cutoff1000;

% Welfare PDVs for 1996-2010
SW_cutoff1000_1996_2010 = SW_cutoff1000(1:180, :);
SWpv_cutoff1000_1996_2010 = Beta_1996_2010' * SW_cutoff1000_1996_2010;

% Welfare PDVs for 2011-2025
SW_cutoff1000_2011_2025 = SW_cutoff1000(181:360, :);
SWpv_cutoff1000_2011_2025 = Beta_2011_2025' * SW_cutoff1000_2011_2025;

% Rename & save

% Save the "baseline" (N = 3) results for comparison with CF policies (later)
base_SW_cutoff500 = SW_cutoff500;
base_SWpv_cutoff500 = SWpv_cutoff500;
base_SWpv_cutoff500_1996_2010 = SWpv_cutoff500_1996_2010;
base_SWpv_cutoff500_2011_2025 = SWpv_cutoff500_2011_2025;

% Save the "baseline" (N = 3) results for comparison with CF policies (later)
base_SW_cutoff1000 = SW_cutoff1000;
base_SWpv_cutoff1000 = SWpv_cutoff1000;
base_SWpv_cutoff1000_1996_2010 = SWpv_cutoff1000_1996_2010;
base_SWpv_cutoff1000_2011_2025 = SWpv_cutoff1000_2011_2025;

% convert to single-precision to save storage space
EVpe = single(EVpe);
Lambda = single(Lambda);
Policy_actual = single(Policy_actual);
Policy_actual_pe = single(Policy_actual_pe);
Policy_pe = single(Policy_pe);
V = single(V);
W = single(W);
EV = single(EV);
Policy = single(Policy);
CS_cutoff500_simu = single(CS_cutoff500_simu);
CS_cutoff1000_simu = single(CS_cutoff1000_simu);
PS_simu = single(PS_simu);
SW_cutoff500 = single(SW_cutoff500);
SW_cutoff1000 = single(SW_cutoff1000);
State_simu = single(State_simu);
State_simu_id = single(State_simu_id);
Entry_simu = single(Entry_simu);
Exit_simu = single(Exit_simu);
Invest_simu = single(Invest_simu);
Merge_simu = single(Merge_simu);
depreciation = single(depreciation);
synergy = single(synergy);
target = single(target);
turn = single(turn);

diary off

save result_KeyResult_180822_sigma055.mat ...
     Theta V_actual V_actual_pe Policy_actual Policy_actual_pe State_avg_prime ...
     State_avg_sum Entry_avg Exit_avg_prime Invest_avg_prime Merge_avg_prime Fron_fit T NS ll ...
     Exit_avg_sum_t248 Invest_avg_sum_t248 Merge_avg_sum_t248 Entry_avg_sum_t248 State_avg_t249 State_avg_sum_t249 ...
     Exit_avg_sum_t359 Invest_avg_sum_t359 Merge_avg_sum_t359 Entry_avg_sum_t359 State_avg_t360 State_avg_sum_t360;
delete result_Theta_sigma055.mat;
save result_more_180822_sigma055.mat EVpe Lambda Policy_actual Policy_actual_pe Policy_pe Theta V V_actual V_actual_pe W;
save result_EV_180822_sigma055.mat -v7.3 EV;
save result_Policy_180822_sigma055.mat -v7.3 Policy;
save result_Simulation_180822_sigma055.mat CS_cutoff500_simu CS_cutoff1000_simu PS_simu SW_cutoff500 SWpv_cutoff500 SW_cutoff1000 SWpv_cutoff1000 State_simu Entry_simu Exit_simu Invest_simu Merge_simu depreciation synergy target turn;
save result_P_avg_180822_sigma055.mat P_avg;
save compareSW_base_newCS_sigma055.mat ...
    base_SW_cutoff500 base_SWpv_cutoff500 base_SWpv_cutoff500_1996_2010 base_SWpv_cutoff500_2011_2025 ...
    base_SW_cutoff1000 base_SWpv_cutoff1000 base_SWpv_cutoff1000_1996_2010 base_SWpv_cutoff1000_2011_2025;
