function [datesCDS, survProbs, intensities] = bootstrapCDS_JT(datesCDS, spreadsCDS, recovery)
% Computes the survival probability and intensity from CDS with JT formula

% INPUT:
%
% datesDF:      dates for which we have  adiscount factor
% discounts:    discount factors which are given
% datesCDS:     dates in which we have a CDS quoted spread
% spreadsCDS:   CDS quoted spread 
% recovery:     recovery rate of the issuer
%
% OUTPUT:
% datesCDS:     dates in which we have a CDS quoted spread
% survProbs:    piecewise constant survival probability for the issuer
% intensities:  piecewise constant intensities for the issuer

% Setting the daycount convention
daycount = 3; % act/365

% Computing intensities and probabilities
intensities = 1e4*spreadsCDS/(1-recovery);
survProbs = exp(cumsum(-1e-4*intensities.*yearfrac(datesCDS(1:end-1),datesCDS(2:end),daycount)));

end