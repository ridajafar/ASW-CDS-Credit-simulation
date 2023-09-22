function [datesCDS, survProbs, intensities] = bootstrapCDS_exact_approx(discounts_CDS, datesCDS, spreadsCDS, flag, recovery)
% Computes the survival probability and intensity from CDS

% INPUT:
%
% discounts_CDS:discount factors for the dates in datesCDS
% datesCDS:     dates in which we have a CDS quoted spread
% spreadsCDS:   CDS quoted spread 
% flag:         1 if we want to compute the approximated version 
%               2 if we want to compute the exact version
% recovery:     recovery rate of the issuer
%
% OUTPUT:
% datesCDS:     dates in which we have a CDS quoted spread
% survProbs:    piecewise constant survival probability for the issuer
% intensities:  piecewise constant intensities for the issuer

% Setting daycount convention
probability_daycount = 6; % 30/360
intensities_daycount = 3; % act/365

% Compute lambda piecewise constant
P_survival = [1; zeros(length(spreadsCDS),1)];
intensities = zeros(length(spreadsCDS),1);

P_survival(2) = (1 - recovery - (flag==2)*yearfrac(datesCDS(1),datesCDS(2),probability_daycount)/2*spreadsCDS(1))/(1 - recovery + spreadsCDS(1)*yearfrac(datesCDS(1),datesCDS(2),probability_daycount)*(1-(flag==2)/2) ); %%  Check yearfrac mode
intensities(1) = -log(P_survival(2))/yearfrac(datesCDS(1),datesCDS(2),intensities_daycount)/1e-4;

for i = 2:length(intensities)

    num = ( (1-recovery)*( sum(discounts_CDS(1:i).*P_survival(1:i)) ...
        - sum((discounts_CDS(1:i-1).*P_survival(2:i))) ) ...
        - spreadsCDS(i)*sum(yearfrac(datesCDS(1:i-1),datesCDS(2:i),probability_daycount).*discounts_CDS(1:i-1).*(P_survival(2:i)...
                            +(flag==2)/2*(P_survival(1:i-1)-P_survival(2:i)))) ...
        - (flag==2)*spreadsCDS(i)*yearfrac(datesCDS(i),datesCDS(i+1),probability_daycount)/2*discounts_CDS(i)*P_survival(i) );

    den = ( spreadsCDS(i)*yearfrac(datesCDS(i),datesCDS(i+1),probability_daycount)*discounts_CDS(i)*(1-(flag==2)/2) ...
        + (1-recovery)*discounts_CDS(i) );

    P_survival(i+1) = num/den;

    intensities(i) = -log(P_survival(i+1)/P_survival(i))/yearfrac(datesCDS(i),datesCDS(i+1),intensities_daycount)/1e-4;

end

survProbs = P_survival(2:end);


end