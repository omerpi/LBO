%%
clc, clear, close all;

cf = [-5201.8, 4550.7];
irr(cf);

%%
Rate         = 0.12/12;   % 12 percent APR = 1 percent per month
NumPeriods   = 30*12;     % 30 years = 360 months
PresentValue = 10e5;

[Principal, Interest, Balance, Payment] = amortize(Rate, ...
NumPeriods, PresentValue);

format bank

Payment

%%
plot(Balance,'b'), hold('on')
plot(cumsum(Principal),'--k')
plot(cumsum(Interest),':r')

xlabel('Payment Month')
ylabel('Dollars')
grid('on')
title('Outstanding Balance, Cumulative Principal & Interest')
legend('Outstanding Balance', 'Cumulative Principal', ... 
'Cumulative Interest')


d_a = 100;
