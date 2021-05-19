#include "mex.h"
#include "array.hpp"

/*

        n = sum(n_old);                     // Number of all firms
        
        % Probabilities of depreciation shock
        pr_hitus = n_old * pr_down;     % Prob's that "we" (particular levels of firms) get hit
        pr_hitme = ones(4,1) ./ n_old;  % Cond'l on hitting "us", prob. that it actually hits "me"
        pr_hitme(n_old == 0) = 0;       % Fool-proof: = 0 if no firm (of my type) exists

        tmp = zeros(1,5);

        % Umbrella values that account for exogenous state change (due to stochastic depreciation)
        tmp(1) = (1 - n * pr_down) * Lambda0(1, 1, state)...
            + pr_hitus(1) * (pr_hitme(1) * Lambda0(1, 1, newstate_d(state, 1))...
                        + (1 - pr_hitme(1)) * Lambda0(1, 1, newstate_d(state, 1)))...
            + pr_hitus(2) * Lambda0(1, 1, newstate_d(state, 2))...
            + pr_hitus(3) * Lambda0(1, 1, newstate_d(state, 3))...
            + pr_hitus(4) * Lambda0(1, 1, newstate_d(state, 4));

        tmp(2) = (1 - n * pr_down) * Lambda0(2, 1, state)...
            + pr_hitus(1) * Lambda0(2, 1, newstate_d(state, 1))...
            + pr_hitus(2) * (pr_hitme(2) * Lambda0(1, 1, newstate_d(state, 2))...
                        + (1 - pr_hitme(2)) * Lambda0(2, 1, newstate_d(state, 2)))...
            + pr_hitus(3) * Lambda0(2, 1, newstate_d(state, 3))...
            + pr_hitus(4) * Lambda0(2, 1, newstate_d(state, 4));

        tmp(3) = (1 - n * pr_down) * Lambda0(3, 1, state)...
            + pr_hitus(1) * Lambda0(3, 1, newstate_d(state, 1))...
            + pr_hitus(2) * Lambda0(3, 1, newstate_d(state, 2))...
            + pr_hitus(3) * (pr_hitme(3) * Lambda0(2, 1, newstate_d(state, 3))...
                        + (1 - pr_hitme(3)) * Lambda0(3, 1, newstate_d(state, 3)))...
            + pr_hitus(4) * Lambda0(3, 1, newstate_d(state, 4));

        tmp(4) = (1 - n * pr_down) * Lambda0(4, 1, state)...
            + pr_hitus(1) * Lambda0(4, 1, newstate_d(state, 1))...
            + pr_hitus(2) * Lambda0(4, 1, newstate_d(state, 2))...
            + pr_hitus(3) * Lambda0(4, 1, newstate_d(state, 3))...
            + pr_hitus(4) * (pr_hitme(4) * Lambda0(3, 1, newstate_d(state, 4))...
                        + (1 - pr_hitme(4)) * Lambda0(4, 1, newstate_d(state, 4)));

        tmp(5) = Lambda0(5,1,state);

        LambdaT(:,state) = tmp;
*/

// stochastic_depreciation(n_old, pr_down, LambdaT, newstate_d, Lambda, S, t);
//
// This function writes the results directly into the Lambda array.
//
void mexFunction(
         int          nlhs,
         mxArray      *plhs[],
         int          nrhs,
         const mxArray *prhs[]
         )  
{
    if (nrhs != 7) {
        mexErrMsgIdAndTxt("MATLAB:stochastic_depreciation:nargin", 
            "MEXCPP requires 7 input arguments.");
    } else if (nlhs != 0) {
        mexErrMsgIdAndTxt("MATLAB:stochastic_depreciation:nargout",
            "MEXCPP requires 0 output argument.");
    }

    ArrayND statespace(prhs[0]);
    // double pr_down = mxGetScalar(prhs[1]);
    ArrayND pr_down(prhs[1]);
    ArrayND Lambda0(prhs[2]);
    ArrayND newstate_d(prhs[3]);
    ArrayND Lambda(prhs[4]);
    int S = mxGetScalar(prhs[5]);
    int t = mxGetScalar(prhs[6]);

    ArrayND LambdaT = ArrayND::create2D(5, 1, plhs[0]);
        
    // Probabilities of depreciation shock
    double _n_old[4];
    double _pr_hitus[4];
    double _pr_hitme[4];

    ArrayND n_old(4,1, _n_old);
    ArrayND pr_hitus(4,1, _pr_hitus);
    ArrayND pr_hitme(4,1, _pr_hitme);

    for(int state=1; state <= S; state++)
    {
        // Number of firms (by type) in current state: [4 x 1]
        for(int i = 1; i <= 4; i++)
            n_old(i) = statespace(state, i+1);  

        double n = sum(n_old);                  // Number of all firms
        double pr_hitanyone = n_old(1) * pr_down(1) + n_old(2) * pr_down(2) 
                            + n_old(3) * pr_down(3) + n_old(4) * pr_down(4);

        for(int i=1; i <= 4; i++)
        {
            pr_hitus(i) = n_old(i) * pr_down(i);
            pr_hitme(i) = 1.0 / n_old(i);

            if(n_old(i) == 0)
                pr_hitme(i) = 0.0;
        }

        // Umbrella values that account for exogenous state change (due to stochastic depreciation)
        LambdaT(1) = (1 - pr_hitanyone) * Lambda0(1, state)
            + pr_hitus(1) * (pr_hitme(1) * Lambda0(1, newstate_d(state, 1))
                        + (1 - pr_hitme(1)) * Lambda0(1, newstate_d(state, 1)))
            + pr_hitus(2) * Lambda0(1, newstate_d(state, 2))
            + pr_hitus(3) * Lambda0(1, newstate_d(state, 3))
            + pr_hitus(4) * Lambda0(1, newstate_d(state, 4));

        LambdaT(2) = (1 - pr_hitanyone) * Lambda0(2, state)
            + pr_hitus(1) * Lambda0(2, newstate_d(state, 1))
            + pr_hitus(2) * (pr_hitme(2) * Lambda0(1, newstate_d(state, 2))
                        + (1 - pr_hitme(2)) * Lambda0(2, newstate_d(state, 2)))
            + pr_hitus(3) * Lambda0(2, newstate_d(state, 3))
            + pr_hitus(4) * Lambda0(2, newstate_d(state, 4));

        LambdaT(3) = (1 - pr_hitanyone) * Lambda0(3, state)
            + pr_hitus(1) * Lambda0(3, newstate_d(state, 1))
            + pr_hitus(2) * Lambda0(3,  newstate_d(state, 2))
            + pr_hitus(3) * (pr_hitme(3) * Lambda0(2, newstate_d(state, 3))
                        + (1 - pr_hitme(3)) * Lambda0(3, newstate_d(state, 3)))
            + pr_hitus(4) * Lambda0(3, newstate_d(state, 4));

        LambdaT(4) = (1 - pr_hitanyone) * Lambda0(4, state)
            + pr_hitus(1) * Lambda0(4, newstate_d(state, 1))
            + pr_hitus(2) * Lambda0(4, newstate_d(state, 2))
            + pr_hitus(3) * Lambda0(4, newstate_d(state, 3))
            + pr_hitus(4) * (pr_hitme(4) * Lambda0(3, newstate_d(state, 4))
                        + (1 - pr_hitme(4)) * Lambda0(4, newstate_d(state, 4)));

        LambdaT(5) = Lambda0(5, state);

        Lambda(1,t,state) = LambdaT(1);
        Lambda(2,t,state) = LambdaT(2);
        Lambda(3,t,state) = LambdaT(3);
        Lambda(4,t,state) = LambdaT(4);
        Lambda(5,t,state) = LambdaT(5);
    }
}