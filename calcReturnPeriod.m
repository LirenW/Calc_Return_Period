function [ X ] = calcReturnPeriod( Q, YEAR )
%  Calculate the return period using PWMs method
%
%     Inputs:
%                  Q - input data, should be [lon lat time]
%               YEAR - the return level
%     Outputs:
%                  X - the return period at specific (YEAR) level
%     
%  Created by Liren
%  May 2018

N             = size(Q, 3);
Q             = sort(Q, 3);

YEAR
[alpha, zeta] = PWMs(Q, true);
para          = 1 - 1./YEAR;
LNF           = -log(para);
X             = zeta - alpha .* (log(LNF));

end
