function [asw]= assetswap(discounts,c,dirtyPriceBond,times)
% Computes the spread in asset swap

% INPUT:
%
% discounts:        discount factors on coupon dates
% c:                coupon rate
% dirtyPriceBond:   dirty price of the bond
% times:            settlement + times in which we receive the coupons
%
% OUTPUT:
% asw:              spread in asset swap

% Setting daycount convention
asw_daycount = 6; % 30/360

% Computing ptice of the bond
BPV = sum(yearfrac(times(1:end-1),times(2:end),asw_daycount).*discounts); % which day convention
PriceBond = c*BPV+discounts(end);

% Computing spread in asw
asw = (PriceBond-dirtyPriceBond)/BPV;

end

