# Make Pi [Period Profits]: Endogenous MC-frontier with 4 types of firms & log-linear demand

# Demand-related Data
P_data = [178.617 	145.795 	139.330 	116.438 	121.133 	102.600 	68.707 	57.098 	42.892 	33.633 	28.683 	23.242 	20.526 	15.173 	11.425 	8.975 	6.7894 	5.2892 	4.3132 	4.4868 	4.1653 	3.6522 	3.1868 	2.4713 	1.9218 	1.7831 	1.6154 	1.3589 	1.2174 	1.1000 	0.9789 	0.9351 	0.7460 	0.6813 	0.6062 	0.5747 	0.5483 	0.5234 	0.4859 	0.4507 	0.4696 	0.4269 	0.3810 	0.3665 	0.3131 	0.2868 	0.2703 	0.2395 	0.2244 	0.2003 	0.1670 	0.1307 	0.1277 	0.1097 	0.1006 	0.0882 	0.0850 	0.0787 	0.0727 	0.0607 	0.0590 	0.0546 	0.0514 	0.0671 	0.0702 	0.0660 	0.0587 	0.0552 	0.0538 	0.0522 	0.0532 	0.0485 	0.0487 	0.0478 	0.0429 	0.0426 	0.0407 	0.0401 	0.0392 	0.0366 	0.0363 	0.0351 	0.0346]' # Price per gigabytes [$/GB]
Q_data = [0.021 	0.023 	0.028 	0.037 	0.042 	0.048 	0.071 	0.082 	0.100 	0.112 	0.143 	0.179 	0.184 	0.217 	0.299 	0.403 	0.477 	0.546 	0.704 	0.803 	0.820 	0.839 	0.992 	1.386 	1.532 	1.395 	1.740 	2.262 	2.285 	2.362 	3.131 	3.748 	3.853 	3.507 	4.694 	6.072 	5.919 	5.685 	6.677 	7.969 	7.439 	7.452 	8.920 	9.950 	10.716 	11.328 	15.207 	17.028 	17.836 	17.513 	21.622 	20.504 	20.508 	24.042 	28.236 	35.009 	38.531 	35.694 	41.213 	48.365 	48.640 	50.074 	53.196 	34.117 	33.635 	36.714 	37.029 	41.619 	40.888 	36.907 	37.878 	42.623 	41.150 	38.560 	45.670 	41.180 	38.85 	31.59 	33.96 	37.07 	32.65 	29.15 	32.28]' # Total exabytes shipped [EB]
X_data = [14.498	14.468	15.15	16.972	15.774	15.546	16.463	18.936	17.397	16.806	17.864	21.415	19.974	19.988	22.047	26.678	23.109	22.85	25.99	28.191	24.597	23.103	23.13	26.739	23.491	21.803	25.712	29.902	25.602	24.879	29.949	34.094	31.445	30.143	32.021	37.289	34.108	33.582	35.448	40.312	35.94	31.84	35.3	36.11	34.05	32.6	36.7	39.16	35.74	35.2	38.9	35.04	31.02	31.51	34.15	35.85	36.1	34.7	36.45	36.9	35.25	35.05	36.81	36.45	33.8	33.3	33.9	35.45	33.1	32.8	34.85	34.55	32.1	32.6	33.3	32.9	29.60 	27.80 	28.70 	28.60 	25.17 	24.80 	25.23]' # Demand shifter: Desktop PC shipment [million]
Thai_data = hcat(zeros(1,63),ones(1,20))'  # Thai-flood dummy [=1 from 2011Q4, | 64th quarter]
T = size(P_data,1)                         # Length of sample period, based on Price data

# Frontier of productivity
F = 32                                     # Number of frontier-MC levels
FrontierGrid = zeros(F,1)
for f = 1:F
    FrontierGrid[f] = 5 - 0.5 * (f - 1)    # Grid of log(MC) = [5.0, 4.5, ..., -10.5]'
end

# Mapping frontier level to its corresponding time period [quarter] in history, 
# for the purpose of initializing the best-response iterations with "close" quantities in data.
FTmap = [1; 2; 7; 9; 11; 14; 15; 17; 19; 24; 27; 31; 35; 43; 48; 52; 56; 62; 78; 63; 63; 63; 63; 63; 63; 63; 63; 63; 63; 63; 63; 63]   

# 2018/7/10 demand estimates with detrended [log_Q & log_P], original [log_X], & "Z & Thai_flood as IV"
Alpha0 = ones(T,1) * -.9202086  # constant term         
Alpha1 = ones(T,1) * -1.043296  # price coefficient [= elasticity]
Alpha2 = ones(T,1) *  .2757225  # demand-shifter coeff.
Alpha3 = ones(T,1) * 0          # Thai flood coeff. (not used in this ver.)
detrend_P0 = 3.900766           # constant term from:   reg log_P time
detrend_P1 = -0.1026352         # time coeff. from:     reg log_P time
detrend_Q0 = -2.378813          # constant term from:   reg log_Q time
detrend_Q1 = 0.0913376          # time coeff. from:     reg log_Q time

# Number of types
Type = 4

# Construct state space: frontier * market structure = (frontier, n1, n2, n3, n4)
maxn = 13                  # Maximum number of firms in the industry = 13 [data in 1996Q1]
using ElasticArrays
statespace = ElasticArray{Int}(undef,Type+1,0)  # Grid of states, each row of which = (frontier, n1, n2, n3, n4)
using Printf
@printf("\n Defining the state space...    \n")
i = 0                      # Index of state
using TimerOutputs
to = TimerOutput()
@timeit to "Step 1-1" begin
for f = 1:F
    for n1 = 0:maxn
        for n2 = 0:(maxn - n1)
            for n3 = 0:(maxn - n1 - n2)
                for n4 = 0:(maxn - n1 - n2 - n3)
                    global i = i + 1
                    append!(statespace, [f, n1, n2, n3, n4])
                end
            end
        end
    end
end
statespace = transpose(statespace)
# Size of the state space
S = size(statespace,1)
@printf("\n The size of the state space is: %5.0f \n",S)
    
# Marginal-cost estimates by firm type: [F x Type] with [frontier, type]
logMC = zeros(F, Type)                         # Grid of MC [by type] in log $
logMC[:, 4] = FrontierGrid                     # Top-level [level-4] MC is equivalent to the frontier
logMC[:, 3] = FrontierGrid + 0.1 * ones(F,1)   # Level-3 MC is 0.1 higher than frontier [in log$]
logMC[:, 2] = FrontierGrid + 0.2 * ones(F,1)   # Level-2 MC is 0.2 higher than frontier [in log$]
logMC[:, 1] = FrontierGrid + 0.3 * ones(F,1)   # Level-1 MC is 0.3 higher than frontier [in log$]
MCgrid = exp.(logMC)                           # Grid of MC [by frontier & type] in $

# Initialize period profit
Pi = zeros(Type,T,S)       # Period profit [type, t, state]
Q1 = zeros(T,S)            # Aggregate output [t, state]
P1 = zeros(T,S)            # Price [t, state]

    using Roots
# Construct period profit in all periods & all state = (f, n1, n2, n3, n4)
for t = 1:T 
    
    @printf("\n Calculating equilibrium profits in period... %5.0f \n",t)

    # Pick up parameters in period t
    alpha0 = Alpha0[t]
    alpha1 = Alpha1[t]
    alpha2 = Alpha2[t]
    alpha3 = Alpha3[t]
    
    # Pick up variables in period t
    x = X_data[t]
    thai = Thai_data[t]
    p = P_data[t]
    
    # For each state = frontier * market structure: (f, n1, n2, n3, n4)
    state = 1               # Initialize state #
 
    for f = 1:F
        equiv_t = FTmap[f]  # Mapping frontier level to corresponding time in history
        MCf = MCgrid[f, :]  # For given 'f', use relevant MC by firm level [x4]
        
        for n1 = 0:maxn
            for n2 = 0:(maxn - n1)
                for n3 = 0:(maxn - n1 - n2)
                    for n4 = 0:(maxn - n1 - n2 - n3)
                                
                        @printf("\n Calculating equilibrium profits in period %2.0f, state %5.0f [of %5.0f]...    ",t,state,S) 
                        n = n1 + n2 + n3 + n4       # Number of all (active) firms

                        if n .== 0
                            # Equilibrium output when no firm exists
                            Qsum = 0                               # Aggregate output
                            p1 = 0                                 # Price
                            q1 = 0                                 # Output of a level-1 firm
                            q2 = 0                                 # Output of a level-2 firm
                            q3 = 0                                 # Output of a level-3 firm
                            q4 = 0                                 # Output of a level-4 firm
                            
                        else()
                            # Computing equilibrium output profile when at least one firm is active

                            # Prepare marginal cost of each [active] firm
                            MC = zeros(1,n)                        # Initialize all (active) firms' MC

                            if n1 .> 0
                                MC[1,1:n1] = MCf[1] * ones(1,n1)   # Each level-1 firm has MC_1
                            end

                            if n2 .> 0
                                MC[1,(n1 + 1):(n1 + n2)] = MCf[2] * ones(1,n2)
                            end

                            if n3 .> 0
                                MC[1,(n1 + n2 + 1):(n1 + n2 + n3)] = MCf[3] * ones(1,n3)
                            end

                            if n4 .> 0
                                MC[1,(n1 + n2 + n3 + 1):(n1 + n2 + n3 + n4)] = MCf[4] * ones(1,n4)
                            end

                            # Preparing the iterative best-response loop
                            
                            Q0 = zeros(1,n)    # Initialize all firms' initial output
                            Q0[1,:] = (sum(Q_data[equiv_t,:]) / n) * ones(1,n)
                            # Initial output is average per-capita output in data 
                            # (dividing the actual aggregate output by N)
                                                        
                            Q1t = Q0
                            iter = 1           # Initialize the # of iterations for while-loop for Nash equilibrium
                            gap = 1000         # Initialize the Nash equilibrium convergence criterion

                            # Iterative best-response loop
                            while gap[1] > .0001

                                # For each firm, calculate optimal output [given rivals']
                                for i = 1:n
                                    mc = MC[i]
                                    q0 = Q1t[i]
                                    Qj = sum(Q1t) - q0
                                    function foc(q)
                                        # Firm's first-order condition for profit maximization; with log-linear demand function
                                        
                                        # Price
                                        Q = max(Qj + q, 0.001)     # max(., .) operation to avoid negative Q while searching for optimal q*
                                        logQ_detrended = log(Q) - (detrend_Q0 + detrend_Q1 * 42)     # detrend Q
                                        logP_detrended = (1 / alpha1) * logQ_detrended .- (alpha2 / alpha1) * log(x) .- (alpha0 / alpha1) .- (alpha3 / alpha1) * thai
                                        logP = logP_detrended .+ (detrend_P0 .+ detrend_P1 * 42)     # retrend P
                                        P = exp(logP)
                                        
                                        # Slope of demand
                                        dPdQ = (alpha1 * Q / P) ^ (-1)
                                        
                                        # FOC: LHS - RHS = 0
                                        P + dPdQ * q - mc
                                        
                                    end
                                    q = max(fzero(foc,q0),0)        # Best-response output of firm i
                                    Q1t[i] = q                      # Record i's best-response output
                                end

                                # Check convergence to Nash equilibrium
                                Gap = Q1t - Q0                      # Difference between previous
                                gap = Gap * Gap'
                                Q0 = Q1t
                                iter = iter + 1
                            end

                            # Equilibrium outcomes
                            Qsum = sum(Q1t)                         # Aggregate output
                            logQ_detrended = log(Qsum) - (detrend_Q0 + detrend_Q1 * 42)
                                                                    # Detrend Q                                                                    
                            logP_detrended = (1 / alpha1) * logQ_detrended - (alpha0 / alpha1) - (alpha2 / alpha1) * log(x) - (alpha3 / alpha1) * thai
                                                                    # Log-linear inverse demand
                            logP = logP_detrended + (detrend_P0 + detrend_P1 * 42)
                                                                    # Retrend P
                            p1 = exp(logP)                          # Price [in $]

                            # Equilibrium output of each type of firms [taking average, just in case]                                    
                            if n1 .> 0
                                q1 = sum(Q1t[1,1:n1]) / n1
                            else()
                                q1 = 0
                            end

                            if n2 .> 0
                                q2 = sum(Q1t[1,(n1 + 1):(n1 + n2)]) / n2
                            else()
                                q2 = 0
                            end

                            if n3 .> 0
                                q3 = sum(Q1t[1,(n1 + n2 + 1):(n1 + n2 + n3)]) / n3
                            else()
                                q3 = 0
                            end

                            if n4 .> 0
                                q4 = sum(Q1t[1,(n1 + n2 + n3 + 1):(n1 + n2 + n3 + n4)]) / n4
                            else()
                                q4 = 0
                            end

                        end
 
                        # Equilibrium profit of each type of firms
                        pi1 = (p1 - MCf[1]) * q1
                        pi2 = (p1 - MCf[2]) * q2
                        pi3 = (p1 - MCf[3]) * q3
                        pi4 = (p1 - MCf[4]) * q4
                        
                        # Record outcomes
                        Q1[t,state] = Qsum
                        P1[t,state] = p1
                        Pi[:,t,state] = [pi1 pi2 pi3 pi4]'
                        state = state + 1      # Next state index
                        
                    end
                end
            end
        end
    end
end
end
@printf("\n ") 

# Fool-proof: replace negative [Q, P, Pi] with zero
Q1[Q1 .< 0] .= 0
P1[P1 .< 0] .= 0
Pi[Pi .< 0] .= 0

using MAT
matwrite("MakePi_EndoFrontier_180712.mat", Dict(
    "T" => T,
    "S" => S,
	"Q1" => Q1,
	"P1" => P1,
    "Pi" => Pi,
    "Alpha0" => Alpha0,
    "Alpha1" => Alpha1,
    "Alpha2" => Alpha2,
    "Alpha3" => Alpha3,
    "X_data" => convert(Array{Float64},X_data),
    "Thai_data" => convert(Array{Float64},Thai_data),
    "detrend_P0" => detrend_P0,
    "detrend_P1" => detrend_P1,
    "detrend_Q0" => detrend_Q0,
    "detrend_Q1" => detrend_Q1,
    "statespace" => convert(Array{Float64},statespace)
); compress = true)

using Logging
io = open("log.txt", "w+")
logger = SimpleLogger(io)
with_logger(logger) do
    @info(to)
end
flush(io)
close(io)