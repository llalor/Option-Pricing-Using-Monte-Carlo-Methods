function [Price, CI] = BermudanPutCond(S0,K,r,T1,T2,sigma,NRepl)
randn('state',0)
muT1 = (r-sigma^2/2)*T1;
siT1 = sigma*sqrt(T1);
Samples = randn(NRepl,1);
PriceT1 = S0*exp(muT1 + siT1*Samples); %getting all the possible prices at node T1
PayoffT1 = exp(-r*T1)*mean(max(K-PriceT1,0)); 

% Pricing put at T1
[calls, puts] = blsprice(PriceT1,K,r,T2-T1,sigma);

% Selecting the best payoff
Values = exp(-r*T1)*max(puts,PayoffT1);
[Price, dummy, CI] = normfit(Values);