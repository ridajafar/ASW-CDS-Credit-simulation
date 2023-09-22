% runAssignment3
% group 7, AY2022-2023
% 

clear;
close all;
clc;
format long;

%% Market data

if ispc()   % Windows version
    formatData='dd/mm/yyyy'; 
    [datesSet, ratesSet] = readExcelData('MktData_CurveBootstrap_AY22-23', formatData);
else        % MacOS version
    datesSet = load("datesSet.mat");
    datesSet = datesSet.datesSet;
    ratesSet = load("ratesSet.mat") ;
    ratesSet = ratesSet.ratesSet;
end

% Bootstrap
[dates, discounts] = BootStrap(datesSet, ratesSet);

%% EXERCISE 1

% Define the parameters
face_value = 1;
dirtyPriceBond = 1.01*face_value;
date1y = datenum(2024,02,02);
times = [dates(1); date1y; dates(13:14)];
c = 3.9/100;
discounts_asw = [Disc_interp(discounts,dates,date1y); discounts(13:14)];

% Compute the spread in asset swap
s_asw = assetswap(discounts_asw,c,dirtyPriceBond,times);


%% EXERCISE 2

% Parameters
recovery = 0.4;
dates_CDS = [dates(1); date1y; dates(13:18)];

% Complete set of CDS spreads
CDS_spreads = [30; 34; 37; 39; 40; 40]*1e-4;
CDS_spreads = interp1(dates_CDS([2:6,8]),CDS_spreads,dates_CDS(2:8),'spline');

% Without the accrual
[ ~, survProbs_approx, intensities_approx] = bootstrapCDS(dates, discounts, dates_CDS, CDS_spreads, 1, recovery);

% With the accrual
[ ~, survProbs_exact, intensities_exact] = bootstrapCDS(dates, discounts, dates_CDS, CDS_spreads, 2, recovery);

% With JT
[ ~, survProbs_JT, intensities_JT] = bootstrapCDS(dates, discounts, dates_CDS, CDS_spreads, 3, recovery);

% Showing that the accrual term is really negligible
disp('The difference in intensities between the approximated version and the exact one is:')
disp(intensities_exact-intensities_approx)

% Showing that the JT intensities are very similar to the progressive mean
mean_int = cumsum(intensities_approx)./(1:7)';
disp('The difference between the intensities computed with JT formula and progressive mean of the approximated intensities is:')
disp(abs(intensities_JT-mean_int))

% Plotting the intensities
plot_piecewise_intensities(intensities_approx,intensities_exact,intensities_JT)


%% EXERCISE 3

%% Part a
% Simulate the default times with fixed seed
seed=42;
rng(seed);
bp = 1e-4;
T_final = 30;
lambda1 = 4 * bp;
lambda2 = 10 * bp;
theta = 5;
M = 1e5;

u = rand(M,1);
tau=(-log(u)./(lambda1)).*(-log(u)./(lambda1)<=theta)+((-log(u)-lambda1.*theta)./(lambda2)+theta).*(-log(u)./(lambda1)>theta);
tau=tau';
u=u';
A=[tau ; u];
tau = sort(tau); %we will need also a sorted vector of tau

%% Part b
% Maximum likelihood estimation and Fisher Matrix confidence
% interval
n1 = length(find(tau<theta));
n2 = length(find(tau>=theta));

% Estimator lambda 1
lambda1_MLE=n1/(sum(tau(find(tau<theta)))+n2*theta);

% Estimator lambda 2
lambda2_MLE=n2/(sum(tau(find(tau>=theta))-theta));

% Confidence Intervals (95%)
I=fisher_matrix([lambda1_MLE,lambda2_MLE],theta);
alpha=0.95;
C_MLE_lambda1=[lambda1_MLE - norminv((1+alpha)/2)/(sqrt(M*I(1,1)))  ,  lambda1_MLE + norminv((1+alpha)/2)/(sqrt(M*I(1,1)))];
C_MLE_lambda2=[lambda2_MLE - norminv((1+alpha)/2)/(sqrt(M*I(2,2)))  ,  lambda2_MLE + norminv((1+alpha)/2)/(sqrt(M*I(2,2)))];

%% Other methods for estimation and CI 

%% 1) Another method for estimation
% We could also compute the estimation in an empirical way:
% We have the tau from the simulation and they have a
% relation with our parameters so we can invert this relation and impose
% equality between empirical cumulative probability and the exact one
index = find(tau<=theta); %remember tau has already been sorted
cumulative_prob = 1-index/M;
lambda1_vec = -log(cumulative_prob)./tau(index);
lambda1_EMP = mean(lambda1_vec);

cumulative_prob = 1-(index(end)+1:(length(tau)-1))/M;
lambda2_vec = (-log(cumulative_prob)-lambda1_EMP*theta)./(tau(index(end)+1:end-1)-theta);
lambda2_EMP = mean(lambda2_vec);

% Confidence Intervals (95%)
I=fisher_matrix([lambda1_EMP,lambda2_EMP],theta);
alpha=0.95;
C_EMP_lambda1=[lambda1_EMP - norminv((1+alpha)/2)/(sqrt(M*I(1,1)))  ,  lambda1_EMP + norminv((1+alpha)/2)/(sqrt(M*I(1,1)))];
C_EMP_lambda2=[lambda2_EMP - norminv((1+alpha)/2)/(sqrt(M*I(2,2)))  ,  lambda2_EMP + norminv((1+alpha)/2)/(sqrt(M*I(2,2)))];


%% 2) Another method for estimation and CI
% We use linear regression.

% Find index which correspond to tau less than theta
index_1 = find(A(1,:)<theta);
A_1 = A(:,index_1);
tau1= A_1(1,:); %default time
u1= A_1(2,:); %survival probability

% Find index which correspond to tau less than 30y (time horizon we are
% interested in)
index_2 = find(A(1,:)>=theta & A(1,:)<30);
A_2 = A(:,index_2);
tau2= A_2(1,:);
u2 = A_2(2,:);

% When creating a linear fit of the survival probability vs default time it
% can be seen that tau vs log(survival probability) is the most linear, in
% compliance with the relation between the two, indeed the quadratic
% coefficient is much more negligible if we use the log of u.
attempt_log1 = fit(tau1', log(u1)', 'poly2');
attempt_lin1 = fit(tau1', (u1)', 'poly2');
attempt_log2 = fit(tau2', log(u2)', 'poly2');
attempt_lin2 = fit(tau2', u2', 'poly2');

% Fit in the two parts to find an estimate of lambda1 and lambda2
fit_tau1 = fitlm(tau1, log(u1));
fit_tau2 = fitlm(tau2, log(u2));

% Converting to array to get the estimation
lambda1_estimate_Slope = - table2array(fit_tau1.Coefficients(2,1));
lambda2_estimate_Slope = - table2array(fit_tau2.Coefficients(2,1));

std_error1 = table2array(fit_tau1.Coefficients(2,2));
std_error2 = table2array(fit_tau2.Coefficients(2,2));

% Computing the 95% Confidence interval
Confidence_interval_1 = [lambda1_estimate_Slope-norminv(0.975)*std_error1, lambda1_estimate_Slope+norminv(0.975)*std_error1];
Confidence_interval_2 = [lambda2_estimate_Slope-norminv(0.975)*std_error2, lambda2_estimate_Slope+norminv(0.975)*std_error2];

% Plot
figure();
plot(fit_tau1);
hold on;
title('Lambda visualization');
xlabel('Years');
ylabel('Default Probability');
plot(fit_tau2);
hold off;
