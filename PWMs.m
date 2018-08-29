function [ alpha, zeta ] = PWMs( Q , flag )
%  Probability-weighted moments (PWMs) estimate procedure.
%
%     Inputs:
%                  Q - input data
%               flag - use / not use PWMs method, mainly for test, True for OPEN
%     Outputs:
%       [alpha zeta] - two parameters of PWMs methods
%     
%     NOTE: Regular moments estimate method is also used in this method

%
%  Created by Liren
%  May 2018
%
%  References:
%
%  Hosking, J. & Wallis, J. Regional Frequency Analysis: An Approach Based on
%               L-moments (Cambridge Univ. Press, 1997).
%
%  Hirabayashi, Y., Kanae, S., Emori, S., Oki, T. and Kimoto, M.: Global 
%               projections of changing risks of floods and droughts in a changing climate, 
%               Hydrol. Sci. J., 53(4), 754–772, doi:10.1623/hysj.53.4.754, 2008.
%
%  Rasmussen, P. F. and Gautam, N.: Alternative PWM-estimators of the gumbel distribution, 
%               J. Hydrol., 280(1–4), 265–271, doi:10.1016/S0022-1694(03)00241-5, 2003.

    Q      = sort(Q, 3);
    N      = size(Q, 3);
    M100   = mean(Q, 3);

    n      = [0: 1: N-1]./(N-1);
    nx     = size(Q, 1);
    ny     = size(Q, 2);
    for ix = 1 : nx
        for iy = 1 : ny
            M110(ix, iy) = 0.;
            for i = 1:N
                M110(ix, iy) = M110(ix, iy) + (Q(ix, iy, i) * n(i));
            end
        end
    end
    M110 = M110 ./ N;

    if flag == true
        L1    = M100;
        L2    = 2 .* M110 - M100;
        alpha = L2 ./ (log(2));
        zeta  = L1 - (alpha .* 0.57721);
    else
        sigmaM = std(Q, 0, 3);
        muM    = mean(Q, 3);
        A      = pi/sqrt(6);
        
        alpha  = A * (1 ./ sigmaM);
        zeta   = muM - (0.57722 .* A) .* sigmaM;
    end
end
