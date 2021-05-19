# Construct CS & PS "master" matrices

## Data Preparation
# Load static parameters & basic data [including period profits]
using MAT
vars = matread("MakePi_EndoFrontier_180712.mat")
Alpha0 = get(vars,"Alpha0",nothing)
Alpha1 = get(vars,"Alpha1",nothing)
Alpha2 = get(vars,"Alpha2",nothing)
Alpha3 = get(vars,"Alpha3",nothing)
detrend_P0 = get(vars,"detrend_P0",nothing)
detrend_P1 = get(vars,"detrend_P1",nothing)
detrend_Q0 = get(vars,"detrend_Q0",nothing)
detrend_Q1 = get(vars,"detrend_Q1",nothing)
P1 = get(vars,"P1",nothing)
Pi = get(vars,"Pi",nothing)
Q1 = get(vars,"Q1",nothing)
S = get(vars,"S",nothing)
statespace = get(vars,"statespace",nothing)
T = get(vars,"T",nothing)
Thai_data = get(vars,"Thai_data",nothing)
X_data = get(vars,"X_data",nothing)
beta = .9 ^ (1/12)  # Monthly discout factor()ac
scaling = 1000      # Re-scale profits to avoid numerical explosion due to exp(pi) & exp(EV)
lambda = 1          # Parameter of synergy distrib'n [Pareto]

# Keep the quarterly variables under new names
Pi_qtr = Pi;      
Alpha0_qtr = Alpha0
Alpha1_qtr = Alpha1
Alpha2_qtr = Alpha2
Alpha3_qtr = Alpha3

# Fixed cost of operation [kappa_c] in data: SGA + CAPEX
Kappa_c_qtr = zeros(4, T);      # Later prepare new version in "Data\Financials\Igami_151012_STX.xlsx' in 'Quarterize" tab
Kappa_c_qtr = Kappa_c_qtr / scaling

# Convert quarterly data to monthly; by simple slicing
Pi = zeros(4, 3*T, S)           # New, monthly profit matrix [Note: T is still quarterly at this point]
Kappa_c = zeros(4, 3*T)         # New, monthly fixed-cost 
X = zeros(1, 3*T)               # New, monthly PC shipments
Thai = zeros(1, 3*T)            # New, monthly Thai-flood dummy
P = zeros(3*T, S)               # New, monthly equilibrium prices
Q = zeros(3*T, S)               # New, monthly equilibrium outputs
Alpha0 = zeros(1, 3*T)          # New, monthly demand parameter estimates
Alpha1 = zeros(1, 3*T)
Alpha2 = zeros(1, 3*T)
Alpha3 = zeros(1, 3*T)

Pi_mon = Pi_qtr ./ 3            # Chopping quarterly profits into equal monthly slices
Kappa_c_mon = Kappa_c_qtr ./ 3  # Chopping quarterly profits into equal monthly slices
X_mon = X_data' ./ 3            # Chopping quarterly PC shipments into equal monthly slices
P_mon = P1                      # Don't chop prices
Q_mon = Q1 ./ 3                 # Chopping quarterly HDD shipments into equal monthly slices

for qtr = 1:T
    Pi[:, 3*(qtr-1)+1:3*(qtr-1)+3, :] = permutedims(repeat(Pi_mon[:, qtr, :], outer = [1 1 3]), [1,3,2])
    Kappa_c[:, 3*(qtr-1)+1:3*(qtr-1)+3] = repeat(Kappa_c_mon[:, qtr], outer = [1 3])
    X[1, 3*(qtr-1)+1:3*(qtr-1)+3] = transpose(fill(X_mon[1, qtr], 3))
    Thai[1, 3*(qtr-1)+1:3*(qtr-1)+3] = transpose(fill(Thai_data[qtr], 3))
    P[3*(qtr-1)+1:3*(qtr-1)+3, :] = permutedims(repeat(P_mon[qtr, :], outer = [1 3]), [2,1])
    Q[3*(qtr-1)+1:3*(qtr-1)+3, :] = permutedims(repeat(Q_mon[qtr, :], outer = [1 3]), [2,1])
    Alpha0[1, 3*(qtr-1)+1:3*(qtr-1)+3] = transpose(fill(Alpha0_qtr[qtr], 3))
    Alpha1[1, 3*(qtr-1)+1:3*(qtr-1)+3] = transpose(fill(Alpha1_qtr[qtr], 3))
    Alpha2[1, 3*(qtr-1)+1:3*(qtr-1)+3] = transpose(fill(Alpha2_qtr[qtr], 3))
    Alpha3[1, 3*(qtr-1)+1:3*(qtr-1)+3] = transpose(fill(Alpha3_qtr[qtr], 3))
end
T = 3 * T                       # T was quarterly; now it's monthly

## Task 3'. Calculate social-welfare outcomes

# Construct CS & PS [for all time & state]
CS = zeros(T,S)                 # Consumer surplus [t,state]
PS = zeros(T,S)                 # Producer surplus [t,state]
using Printf
for t = 1:T
    
    @printf("\n Constructing CS & PS [in all states] for period %3.0f of %3.0f", t, T)
    
    for state = 1:S
        
        f = statespace[state, 1]    # Frontier level
        N = statespace[state, 2:5]  # [1 x 4] vector of # firms: (n1, ..., n4)
        n = sum(N, dims = 2)        # # all firms
        
        # Consumer surplus [CS]
        x = X[t]
        thai = Thai[t]
        alpha0 = Alpha0[t]
        alpha1 = Alpha1[t]
        alpha2 = Alpha2[t]
        alpha3 = Alpha3[t]
        price = P[t,state]          # Equilibrium price
        qty = Q[t,state]            # Equilibrium [aggregate] output
        
        p_cut = 1000                # The very first 3.5'' HDDs in 1982 had the average price of $759.44/unit & could hold only 8.612MB, which translates into "price/GB" of [$759.44 / 8.612MB] x 1000 = 88.1839 x 1000 = 88,183.9 [$/GB]. Their 1982 shipment was 0.0015 million [= 1,500] units, which is practically zero.
        A = detrend_Q0 + detrend_Q1 * 42    # For de-/re-meaning of Q
        B = detrend_P0 + detrend_P1 * 42    # For de-/re-meaning of P
        logP_detrended = log(p_cut) - B     # De-meaning P
        logQ_detrended = alpha0 + alpha1 * logP_detrended + alpha2 * log(x) + alpha3 * thai
        logQ = logQ_detrended + A           # Re-meaning Q
        q_cut = exp(logQ)
        omega = alpha0 + A - alpha1 * B     # Adjusting "alpha0" for the re-meaning of P & Q
        
        CS[t,state] = (alpha1 / (1 + alpha1)) * (x ^ (-alpha2 / alpha1)) * exp(-omega / alpha1) * exp(-alpha3 / alpha1 * thai) * ((qty ^ ((1 + alpha1) / alpha1)) - (q_cut ^ ((1 + alpha1) / alpha1))) - price * (qty - q_cut)
        
        # Producer surplus [PS]
        pi = Pi[:, t, state]        # [4 x 1] vector of profit by type()
        PS[t,state] = N' * pi       # [1 x 4] x [4 x 1] = [1 x 1]: Industry-wide profit
        
    end
end

# Fool-proof: replace negative values with zeros()
CS[CS .< 0] .= 0
PS[PS .< 0] .= 0

# More fool-proofing
CS[:,1] = zeros(T,1)                # Correct "CS = -Inf" when N = 0 [in state #1]

## Linear extrapolation thru Dec-2025

#T2025 = T + 111
T2025 = 360

CS2015 = CS
CS2025 = zeros(T2025, S)
CS2025[1:T, :] = CS2015

PS2015 = PS
PS2025 = zeros(T2025, S)
PS2025[1:T, :] = PS2015

LinearDecline = zeros(111,1)
for t = 1:111
    LinearDecline[t] = 1 - (1/111) * t
    CS2025[(T+t), :] = CS2025[T, :] .* LinearDecline[t]
    PS2025[(T+t), :] = PS2025[T, :] .* LinearDecline[t]
end

CS_cutoff1000 = CS2025

matwrite("CS_cutoff1000.mat", Dict(
    "CS_cutoff1000" => CS_cutoff1000
); compress = true)