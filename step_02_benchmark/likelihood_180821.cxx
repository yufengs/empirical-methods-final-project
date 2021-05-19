#include "mex.h"

#define WITH_BLAS 1
#include "array.hpp"

#define USE_THREADS 1

#if _WIN32
__declspec(thread) ArrayPool * arrayPool = nullptr;
#else
__thread ArrayPool * arrayPool = nullptr;
#endif


#if USE_THREADS
#include "workers.hpp"

ThreadPool * threadPool = nullptr;
#endif

void copysub2d(ArrayND &out, int t, int s, const ArrayND &src)
{
    mexAssert(out.ndims > 2);

    if(out.ndims == 3)
    {
        for(int i = 1; i <= src.dims[0]; i++)
            out(i,t,s) = src(i);
    }

    if(out.ndims == 4)
    {
        for(int j = 1; j <= src.dims[1]; j++)
        {
            for(int i = 1; i <= src.dims[0]; i++)
                out(i,j,t,s) = src(i,j);
        }
    }
}


void runIt(
         int          nlhs,
         mxArray      *plhs[],
         int          nrhs,
         const mxArray *prhs[],
         int startS, int endS
         )  
{
    if (nrhs != 31) {
        mexErrMsgIdAndTxt("MATLAB:likelihood7:nargin", 
            "MEXCPP requires 31 input arguments.");
    } else if (nlhs != 0) {
        mexErrMsgIdAndTxt("MATLAB:likelihood7:nargout",
            "MEXCPP requires 0 output arguments.");
    }

    ArrayND statespace(prhs[0]);
    double maxn = mxGetScalar(prhs[1]);
    ArrayND foo(prhs[2]);
    double beta = mxGetScalar(prhs[3]);
    double kappa_x = mxGetScalar(prhs[4]);
    ArrayND kappa_c_est(prhs[5]);
    double kappa_i = mxGetScalar(prhs[6]);
    double kappa_e = mxGetScalar(prhs[7]);
    double kappa0_m = mxGetScalar(prhs[8]);
    int t = mxGetScalar(prhs[9]);
    double delta = mxGetScalar(prhs[10]);
    ArrayND m_level(prhs[11]);
    ArrayND synergy(prhs[12]);
    ArrayND newstate_i(prhs[13]);
    ArrayND newstate_m(prhs[14]);
    ArrayND newstate_e(prhs[15]);
    ArrayND newstate_x(prhs[16]);
    ArrayND Pi(prhs[17]);
    //int S = mxGetScalar(prhs[18]);
    int T = mxGetScalar(prhs[19]);
    ArrayND Lambda(prhs[20]);
    ArrayND V(prhs[21]);
    ArrayND W(prhs[22]);
    ArrayND EV(prhs[23]);
    ArrayND Policy(prhs[24]);
    ArrayND Policy_pe(prhs[25]);
    ArrayND EVpe(prhs[26]);
    ArrayND Vpe(prhs[27]);
    double sigma = mxGetScalar(prhs[28]);
    ArrayND b_level(prhs[29]);
    double kappa_i4 = mxGetScalar(prhs[30]);

    //mexPrintf("state:%d -> %d\n", startS, endS);

    if(arrayPool == nullptr)
        arrayPool = new ArrayPool();

    for(int state=startS; state <= endS; state++)
    {
        arrayPool->resize(100*1024*1024);

        int f_old = statespace(state, 1);       // Frontier level in current state
        ArrayND kappa_c(4, 1);                  // Fixed cost (by own level) given the current frontier
        ArrayND n_old(4, 1);                    // Number of firms (by type) in current state: [4 x 1]

        for(int i = 1; i <= 4; i++)
        {
            kappa_c(i) = kappa_c_est(i, f_old); 
            n_old(i) = statespace(state, i+1);  
        }

        double n = sum(n_old);                  // Number of all firms
        
        // Alternative "state_prime" to manage exogenous shifts of frontier (Corrected in 8/21/2018, to reflect the latest "Monthly Data" as of 7/10/2018)
        int state2 = state;
        if(t == 18 || t == 24 || t == 30 || t == 39 || t == 42 || t == 48 || t == 69 || t == 78 || t == 102 || t == 141 || t == 183 || t == 231)
        {
         state2 = newstate_i(state, 4);
        }

      // "Active" value functions (V) & optimal choice probabilities (Policy)
        
        // Alternative-specific expected values of exit, stay-alone, & innovation (by firm type)
        ArrayND ev_x = zeros(4,1);                          // Zero value once you exit (absorbing state)
        ArrayND ev_c = Lambda.slice(Range(1,4),t+1,state2);  // If idle, the same number of firms tomorrow
        ArrayND ev_i = zeros(4,1);
        ArrayND Kappa_x = zeros(4,1);                       // Distinguish kappa_x by level
        for(int i = 1; i <= 4; i++)                         // Innovator's level
        {
            Kappa_x(i) = kappa_x;                           
        }
        // if (t < 144)                                        // Higher types never exit before great recession
        // {
        //    Kappa_x(3) = 999;                               
        //    Kappa_x(4) = 999;
        // }
        ArrayND Kappa_i = zeros(4,1);                       // Distinguish kappa_i by level
        for(int i = 1; i <= 3; i++)                         // Innovator's level
        {
            ev_i(i) = Lambda(i+1, t+1, newstate_i(state2, i));   // If innovate, your level goes up by 1
            Kappa_i(i) = kappa_i;                               // Lesser firms' kappa_i
        }
        // No relative change (but frontier moves up) if your level is already 4
        ev_i(4) = Lambda(4, t+1, newstate_i(state2, 4));     
        Kappa_i(4) = kappa_i4;                              // Frontier firm's kappa_i
            
        // Optimal choice probabilities of exit, stay-alone, & innovation (numerators to calculate such prob's)
        ArrayND numer_x = exp((-kappa_c - Kappa_x + beta * ev_x) / sigma);   // Numerator for pr(exit)
        ArrayND numer_c = exp((-kappa_c + beta * ev_c) / sigma);                         // Numerator for pr(idle)
        ArrayND numer_i = exp(((-kappa_c - Kappa_i) + beta * ev_i) / sigma);             // Numerator for pr(inno)
        
        // Gov't policy toward merger alters effective kappa_m
        double kappa_m;

        if (n <= 3)
            kappa_m = 999;          // Gov't bans mergers if N = 3 or less
        else
            kappa_m = kappa0_m + delta * t;      // Cost of merger proposal 

        ArrayND ev_m(4,4);
        ArrayND ev_b(4,4);

        for(int a_level=1; a_level <=4; a_level++)
        {
            for(int t_level=1; t_level <=4; t_level++)
            {  
                // If acquirer and/or target types don't exist, ev_m = 0
                if( n_old(a_level) == 0 || n_old(t_level) == 0 || (n_old(t_level) == 1 && t_level == a_level) )
                {
                    ev_m(a_level,t_level) = 0;    
                    ev_b(a_level,t_level) = 0;
                }
                else
                {
                    ev_m(a_level,t_level) = 
                        synergy(1) * Lambda(m_level(a_level, t_level, 1), t+1,
                                        newstate_m(state2, a_level, t_level, m_level(a_level, t_level, 1)))
                      + synergy(2) * Lambda(m_level(a_level, t_level, 2), t+1,
                                        newstate_m(state2, a_level, t_level, m_level(a_level, t_level, 2)))
                      + synergy(3) * Lambda(m_level(a_level, t_level, 3), t+1,
                                        newstate_m(state2, a_level, t_level, m_level(a_level, t_level, 3)))
                      + synergy(4) * Lambda(m_level(a_level, t_level, 4), t+1,
                                        newstate_m(state2, a_level, t_level, m_level(a_level, t_level, 4)));
                    ev_b(a_level,t_level) = 
                        synergy(1) * Lambda(b_level(a_level, t_level, 1), t+1,
                                        newstate_m(state2, a_level, t_level, b_level(a_level, t_level, 1)))
                      + synergy(2) * Lambda(b_level(a_level, t_level, 2), t+1,
                                        newstate_m(state2, a_level, t_level, b_level(a_level, t_level, 2)))
                      + synergy(3) * Lambda(b_level(a_level, t_level, 3), t+1,
                                        newstate_m(state2, a_level, t_level, b_level(a_level, t_level, 3)))
                      + synergy(4) * Lambda(b_level(a_level, t_level, 4), t+1,
                                        newstate_m(state2, a_level, t_level, b_level(a_level, t_level, 4)));
                }
            }
        }

        // Merger proposals under the alternative specification of Nash Bargaining with 50-50 surplus split
        ArrayND offer_base = beta * ev_c;           // A merger proposal offers the target's outside option
        ArrayND offer_plus = 0.0 * beta * (ev_m - ev_c * ones(1,4));    // 0% of Acquirer's gains (= TIOLI)
        // ArrayND offer_plus = 0.5 * beta * (ev_m - ev_c * ones(1,4));    // 50% of Acquirer's gains

        // Optimal choice probabilities of mergers with levels 1 thru 4 (numerators to calculate such prob's)
        ArrayND numer_m = exp(((-kappa_c - kappa_m * ones(4,1)) * ones(1,4) - ones(4,1) * transpose(offer_base) - offer_plus + beta * ev_m) / sigma); // [4 x 4]
        ArrayND numer_b = exp(((-kappa_c - kappa_i * ones(4,1) - kappa_m * ones(4,1)) * ones(1,4) - ones(4,1) * transpose(offer_base) - offer_plus + beta * ev_b) / sigma); // [4 x 4]
        ArrayND n_old_4x4 = repmat(transpose(n_old), 4, 1);   // [4 x 4] indicates irrelevant elements (targets don't exist)
        n_old_4x4 = n_old_4x4 - eye(4);     // Subtract 1 from # of target firms when a_level == t_level
                                            // (i.e., on the diagonal)
        n_old_4x4(n_old_4x4 < 0) = 0;
        numer_m(n_old_4x4 == 0) = 0;        // Eliminating CCP of mergers with non-existing target types
        numer_b(n_old_4x4 == 0) = 0;
        
        // Choice probabilities of exit, stay-alone, innovate, & mergers
        ArrayND denom = numer_x + numer_c + numer_i + sum(numer_m, 2) + sum(numer_b, 2);            // [4 x 1]
        ArrayND ccp = horzcat(numer_x, numer_c, numer_i, numer_m, numer_b) / (denom * ones(1,11));  // [4 x 11]

        // Fool-proof: Cleaning ccp (= 0 when the decision-maker's type does not exist)
        ArrayND n_old_4x11 = repmat(n_old, 1, 11);      // [4 x 11] indicator of irrelevant elements
        ccp(n_old_4x11 == 0) = 0;
        
        // Record alternative-specific values and policies (=0 if n==0)
        ArrayND EV_val = horzcat(ev_x, ev_c, ev_i, ev_m, ev_b);     // [4 x 11 (x 1 x 1)]
        ArrayND Policy_val = ccp;                                   // [4 x 11 (x 1 x 1)]

        copysub2d(EV, t+1, state, EV_val);
        copysub2d(Policy, t, state, Policy_val);



        // Record "active" (ex-ante) value functions 
        ArrayND VT = Pi.slice(Range(1,Pi.dims[0]), t, state) + sigma * (.5772 * ones(4,1) + log(denom));// [4 x 1]

        // Fool-proof: Cleaning V (= 0 if the type doesn't exist in the current state)
        for(int type = 1; type <= 4; type++)
        {
            if(n_old(type) == 0)
                VT(type) = 0;
        }
        copysub2d(V, t, state, VT);

        // Potential entrants' (EV, Policy, V)
        if (n >= 13)
            kappa_e = 999;  // No entry allowed when industry is "full"
        else
            kappa_e = kappa_e;

        double ev_o = Lambda(5, t+1, state2);               // If stay-out, same state tomorrow
        double ev_e = Lambda(1, t+1, newstate_e(state2));   // If enter, become level-1
        double numer_o = exp((beta * ev_o) / sigma);                 // Numerator for prob of stay-out
        double numer_e = exp((- kappa_e + beta * ev_e) / sigma);     // Numerator for prob of entry
        double denom_pe = numer_o + numer_e;               // denominator of CCP (scalar)
        ArrayND ccp_pe = make2x1(numer_o / denom_pe, numer_e / denom_pe);

        ArrayND EVpe_val = vertcat(ev_o, ev_e);         // [2 (x1 x1)]
        copysub2d(EVpe, t+1, state, EVpe_val);

        copysub2d(Policy_pe, t, state, ccp_pe);

        double Vpe_val = sigma * (.5772 + log(denom_pe)); // Value of pot. ent. (scalar)        

        Vpe(t,state) = Vpe_val;      

        //// "Passive" value functions (W)
        // Rival turn-to-move probabilities (conditional on not my turn, how likely each of the rival types moves)
        ArrayND pr_rivalturn = (ones(5,1) * transpose(vertcat(n_old, 1)) - eye(5)) / ((maxn - 1) * ones(5,5));
            // (my type, rival type): [5 x 5]
            // The numerator = the # firms of each type (minus 1 for my own type)
            // The denominator = the # firms of all types (minus 1, to exclude myself)
                
        // Probabiblities that I become the merger target of the rival
        ArrayND pr_targetme4 = ones(4,4) / (n_old * ones(1,4) - eye(4));   // (my type, rival type): [4 x 4]
            // The numerator is one (= me).
            // The denominator is the # firms of my type (minus 1 if the acquirer is also of my type)

        pr_targetme4(n_old * ones(1,4) - eye(4) < 1) = 0;           // Fool-proof: Avoid division by 0 or -1
                                                                    // (i.e., replace Inf & -1 with 0)

        ArrayND pr_targetme = vertcat(horzcat(pr_targetme4, zeros(4,1)), zeros(1,5));
                                                                    // Pot. ent. can't be an acquirer
                                                                    // Pot. ent. can't be a target
        
        // W_ante: i's values, cond'l on rival types' turn-to-move (but *before* j's choice)
        // W_post: i's values, cond'l on j's choice (passive version of "alternative-specific Vs")

        ArrayND W_ante = ArrayND(5, 5);
        ArrayND W_post = ArrayND(13, 1);

        for(int i = 1; i <= 5; i++)
        {
            for(int j = 1; j <= 4; j++)
            {           
                // Prepare W_post
                W_post(1) = beta * Lambda(i, t+1, newstate_x(state2, j));    // If rival j exits
                W_post(2) = beta * Lambda(i, t+1, state2);                   // If rival j stays alone
                W_post(3) = beta * Lambda(i, t+1, newstate_i(state2, j));    // If rival j innovates

                W_post(4) = beta * (synergy(1) * Lambda(i, t+1, newstate_m(state2, j, 1, m_level(j, 1, 1)))
                                  + synergy(2) * Lambda(i, t+1, newstate_m(state2, j, 1, m_level(j, 1, 2)))
                                  + synergy(3) * Lambda(i, t+1, newstate_m(state2, j, 1, m_level(j, 1, 3)))
                                  + synergy(4) * Lambda(i, t+1, newstate_m(state2, j, 1, m_level(j, 1, 4))));
                                                                   // If j merges with level-1

                W_post(5) = beta * (synergy(1) * Lambda(i, t+1, newstate_m(state2, j, 2, m_level(j, 2, 1)))
                                  + synergy(2) * Lambda(i, t+1, newstate_m(state2, j, 2, m_level(j, 2, 2)))
                                  + synergy(3) * Lambda(i, t+1, newstate_m(state2, j, 2, m_level(j, 2, 3)))
                                  + synergy(4) * Lambda(i, t+1, newstate_m(state2, j, 2, m_level(j, 2, 4))));
                                                                   // If j merges with level-2

                W_post(6) = beta * (synergy(1) * Lambda(i, t+1, newstate_m(state2, j, 3, m_level(j, 3, 1)))
                                  + synergy(2) * Lambda(i, t+1, newstate_m(state2, j, 3, m_level(j, 3, 2)))
                                  + synergy(3) * Lambda(i, t+1, newstate_m(state2, j, 3, m_level(j, 3, 3)))
                                  + synergy(4) * Lambda(i, t+1, newstate_m(state2, j, 3, m_level(j, 3, 4))));
                                                                   // If j merges with level-3

                W_post(7) = beta * (synergy(1) * Lambda(i, t+1, newstate_m(state2, j, 4, m_level(j, 4, 1)))
                                  + synergy(2) * Lambda(i, t+1, newstate_m(state2, j, 4, m_level(j, 4, 2)))
                                  + synergy(3) * Lambda(i, t+1, newstate_m(state2, j, 4, m_level(j, 4, 3)))
                                  + synergy(4) * Lambda(i, t+1, newstate_m(state2, j, 4, m_level(j, 4, 4))));
                                                                   // If j merges with level-4

                W_post(8) = beta * (synergy(1) * Lambda(i, t+1, newstate_m(state2, j, 1, b_level(j, 1, 1)))
                                  + synergy(2) * Lambda(i, t+1, newstate_m(state2, j, 1, b_level(j, 1, 2)))
                                  + synergy(3) * Lambda(i, t+1, newstate_m(state2, j, 1, b_level(j, 1, 3)))
                                  + synergy(4) * Lambda(i, t+1, newstate_m(state2, j, 1, b_level(j, 1, 4))));
                                                                   // If j merges with level-1

                W_post(9) = beta * (synergy(1) * Lambda(i, t+1, newstate_m(state2, j, 2, b_level(j, 2, 1)))
                                  + synergy(2) * Lambda(i, t+1, newstate_m(state2, j, 2, b_level(j, 2, 2)))
                                  + synergy(3) * Lambda(i, t+1, newstate_m(state2, j, 2, b_level(j, 2, 3)))
                                  + synergy(4) * Lambda(i, t+1, newstate_m(state2, j, 2, b_level(j, 2, 4))));
                                                                   // If j merges with level-2

                W_post(10) = beta * (synergy(1) * Lambda(i, t+1, newstate_m(state2, j, 3, b_level(j, 3, 1)))
                                  + synergy(2) * Lambda(i, t+1, newstate_m(state2, j, 3, b_level(j, 3, 2)))
                                  + synergy(3) * Lambda(i, t+1, newstate_m(state2, j, 3, b_level(j, 3, 3)))
                                  + synergy(4) * Lambda(i, t+1, newstate_m(state2, j, 3, b_level(j, 3, 4))));
                                                                   // If j merges with level-3

                W_post(11) = beta * (synergy(1) * Lambda(i, t+1, newstate_m(state2, j, 4, b_level(j, 4, 1)))
                                  + synergy(2) * Lambda(i, t+1, newstate_m(state2, j, 4, b_level(j, 4, 2)))
                                  + synergy(3) * Lambda(i, t+1, newstate_m(state2, j, 4, b_level(j, 4, 3)))
                                  + synergy(4) * Lambda(i, t+1, newstate_m(state2, j, 4, b_level(j, 4, 4))));
                                                                   // If j merges with level-4
                                                                            
                // Adjustments when my (firm i's) type becomes firm j's merger target
                if(i < 5)
                {
                    double ones[11]  = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
                    double zeros[11] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

                    Array2D adjust01 = Array2D(11, 1, ones);
                    Array2D adjust02 = Array2D(11, 1, zeros);

                    adjust01(i + 3) = adjust01(i + 3) - pr_targetme(i, j);  // When someone else becomes target
                    adjust01(i + 7) = adjust01(i + 7) - pr_targetme(i, j);  // When someone else becomes target
                    adjust02(i + 3) = pr_targetme(i, j) * (offer_base(i) + offer_plus(j,i)); 
                    adjust02(i + 7) = pr_targetme(i, j) * (offer_base(i) + offer_plus(j,i)); 
                                                                            //  When I (firm i) become target

                    for(int k = 1; k <= 11; k++)
                      W_post(k,1) = W_post(k,1) * adjust01(k) + adjust02(k);

                    //W_post(1:11, 1) = W_post(1:11, 1) .* adjust01 + adjust02;
                }

                // Calculate W_ante
                //W_ante(i, j) = [ccp(j,:), 0, 0] * W_post; // [1 x 1] = [1 x 13] x [13 x 1] 

                double dp = 0.0;

                for(int k = 1; k <= 11; k++)
                  dp += ccp(j,k) * W_post(k); // last 2 elems in W_post are ignored

                W_ante(i, j) = dp;

            } // for j
            
            if(i < 5)
            {
                int j = 5;  // Note: Because only one potential entrant exists, "i = j = 5" cannot happen.
                
                // Prepare W_post
                W_post(12) = beta * Lambda(i, t+1, state2);                  // If pot. ent. j stays out
                W_post(13) = beta * Lambda(i, t+1, newstate_e(state2));      // If pot. ent. j enters
                
                // Calculate W_ante 
                W_ante(i,j) = ccp_pe(1)*W_post(12) + ccp_pe(2)*W_post(13);
                //W_ante(i, j) = [zeros(1, 11), Policy_pe(:, t, state)'] * W_post;    
            }
          
        } // for i

        // Probability that Nature chooses nobody
        double pr_nobody = 1 - ((n + 1) / maxn);

        // Adjustment for such cases of "nobody's turn"
        ArrayND W_zero = ArrayND(5,1);
        
        for(int i = 1; i <= 5; i++)
        {
            W_zero(i) = pr_nobody * beta * Lambda(i, t+1, state2);
        }

        // Calculate the "passive" value function(s)
        ArrayND WT = vertcat(Pi.slice(Range(1,Pi.dims[0]), t, state), 0) - vertcat(kappa_c, 0) + sum(smul(pr_rivalturn, W_ante), 2) + W_zero;
        
        // Fool-proof (4): Cleaning W (= 0 if the type doesn't exist in the current state)
        for(int type = 1; type <= 4; type++)
        {
            if(n_old(type) == 0)
                WT(type) = 0;
        }
        copysub2d(W, t, state, WT);

        //// Umbrella value functions (Lambdas) = weighted sums of active & passive values (V & W)
        
        // My turn-to-move probability
        ArrayND pr_myturn4 = zeros(4,1);            // my type: [4 x 1] (later becomes [5 x 1]; see below)
        if(n > 0)
        {
            pr_myturn4 = ones(4,1) / maxn;               
            pr_myturn4(n_old == 0) = 0;             // Fool-proof: = 0 if no firm (of that type) exists
        }

        ArrayND pr_myturn = vertcat(pr_myturn4, 1/maxn);    // Pot. ent. ("type 5") always exists


        // Umbrella value functions (Lambdas) are weighted sums of active & inactive values (V & W)
        ArrayND LambdaT = smul(pr_myturn, vertcat(VT, Vpe_val)) + smul(ones(5,1) - pr_myturn, WT);

        // Fool-proof: Cleaning Lambda (= 0 if the type doesn't exist in the current state)
        for(int type = 1; type <= 4; type++)
        {
            if(n_old(type) == 0)
                LambdaT(type) = 0;
        }

        copysub2d(Lambda, t, state, LambdaT);
    }
}


#if USE_THREADS
void mexFunction(
                 int          nlhs,
                 mxArray      *plhs[],
                 int          nrhs,
                 const mxArray *prhs[]
                 )
{
    int S = mxGetScalar(prhs[18]);
    
    int dS = S / 4;

	if(threadPool == nullptr)
		threadPool = new ThreadPool(4);
    
    auto t1 = threadPool->enqueue([&]() {
        runIt(nlhs, plhs, nrhs, prhs, (dS*0)+1, (dS*1));
    });
    
    auto t2 = threadPool->enqueue([&]() {
        runIt(nlhs, plhs, nrhs, prhs, (dS*1)+1, (dS*2));
    });
    
    auto t3 = threadPool->enqueue([&]() {
        runIt(nlhs, plhs, nrhs, prhs, (dS*2)+1, (dS*3));
    });
    
    auto t4 = threadPool->enqueue([&]() {
        runIt(nlhs, plhs, nrhs, prhs, (dS*3)+1, S);
    });
    
    t1.get();
    t2.get();
    t3.get();
    t4.get();
}
#else
void mexFunction(
                 int          nlhs,
                 mxArray      *plhs[],
                 int          nrhs,
                 const mxArray *prhs[]
                 )
{
    int S = mxGetScalar(prhs[18]);
    
    runIt(nlhs, plhs, nrhs, prhs, 1, S);
}
#endif

void foo()
{
    ArrayND a(2,2);
    a.dump();  // force linkage of dump function
}
