function [datesCDS, survProbs, intensities] = bootstrapCDS(datesDF, discounts, datesCDS, spreadsCDS, flag, recovery)
% Computes the survival probability and intensity from CDS

% INPUT:
%
% datesDF:      dates for which we have  adiscount factor
% discounts:    discount factors which are given
% datesCDS:     dates in which we have a CDS quoted spread
% spreadsCDS:   CDS quoted spread 
% flag:         1 if we want to compute the approximated version 
%               2 if we want to compute the exact version 
%               3 if we want to compute the JT version
% recovery:     recovery rate of the issuer
%
% OUTPUT:
% datesCDS:     dates in which we have a CDS quoted spread
% survProbs:    piecewise constant survival probability for the issuer
% intensities:  piecewise constant intensities for the issuer

% Compute the needed discounts
discounts_CDS = [Disc_interp(discounts,datesDF,datenum(2024,02,02)); discounts(13:18)];

% Compute survival probability and intensities for the needed version
if flag==1 || flag==2
    [datesCDS, survProbs, intensities] = bootstrapCDS_exact_approx(discounts_CDS, datesCDS, spreadsCDS, flag, recovery);
elseif flag==3
    [datesCDS, survProbs, intensities] = bootstrapCDS_JT(datesCDS, spreadsCDS, recovery);
end

end