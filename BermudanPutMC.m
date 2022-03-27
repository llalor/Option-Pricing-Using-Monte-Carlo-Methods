function [Price, CI] = BermudanPutMC(S0,K,r,T1,T2,sigma,NRepl1,NRepl2)
% compute auxiliary quantities outside the loop
randn('state',0)
DeltaT = T2-T1;
muT1 = (r-sigma^2/2)*T1;
muT2 = (r-sigma^2/2)*(T2-T1);
siT1 = sigma*sqrt(T1);
siT2 = sigma*sqrt(T2-T1);

% vector to contain payoffs
DiscountedPayoffs = zeros(NRepl1*NRepl2, 1);
% sample at time T1
Samples1 = randn(NRepl1,1); %repititions at T1
PriceT1 = S0*exp(muT1 + siT1*Samples1); %Price for all the samples at T1

for k=1:NRepl1 %for k, 1 to all repititions in 1
    Samples2 = randn(NRepl2,1); %Repititions in T2
    PriceT2 = PriceT1(k)*exp(muT2 + siT2*Samples2); %all the prices here
    ValuePut = exp(-r*DeltaT)*mean(max(K-PriceT2,0));
    PayoffT1 = exp(-r*T1)*mean(max(K-PriceT1,0));
    if ValuePut > PayoffT1
        DiscountedPayoffs(1+(k-1)*NRepl2:k*NRepl2) = ...
            exp(-r*T2)*max(K-PriceT2, 0);
    else 
        DiscountedPayoffs(1+(k-1)*NRepl2:k*NRepl2) = ...
        exp(-r*T1)*PayoffT1;
    end
end 

[Price, dummy, CI] = normfit(DiscountedPayoffs);

