%% Project 
clc
clear
%% Inputs
S0 = 100;
K = 105;
r = 0.05;
sigma = 0.2;
T1 = 0.5;
T2 = 1;
NRepl1 = 100;
NRepl2 = 100;
NRepl = NRepl1*NRepl2;
Smax = 200;
T=1;
M=100;
N=102;
NSteps = 52; %2 for twice exercisable
omega=1; %Overrelaxation parameter which can be adjusted between 0 and 2. Algorithm likely to go bust otherwise. 
tol = 0.001;% Specified Tolerance for the Gauss Siedel iteration procedure
fhandles = {@(x)ones(length(x),1), @(x)x, @(x)x.^2};
ExerciseTimes = 52; %2 for twice exercisable
Ex = ExerciseTimes;
%% Run function
MC = BermudanPutMC(S0,K,r,T1,T2,sigma,NRepl1,NRepl2); 
MC2Cond = BermudanPutCond(S0,K,r,T1,T2,sigma,NRepl);
table(MC,MC2Cond)
LongstaffSchwartzMC = GenericLS(S0,K,r,T,sigma,NSteps,NRepl,fhandles);
Implicit = BermudanPutImplicit(S0,K,r,T,sigma,Smax,M,N,omega,tol,Ex);
table(LongstaffSchwartzMC, Implicit)

%% Brute Force Monte Carlo

DeltaT = T2-T1;
muT1 = (r-sigma^2/2)*T1;
muT2 = (r-sigma^2/2)*(T2-T1);
siT1 = sigma*sqrt(T1);
siT2 = sigma*sqrt(T2-T1);

% vector to contain payoffs
DiscountedPayoffs = zeros(NRepl1*NRepl2, 1);
% sample at time T1
Samples1 = randn(NRepl1,1); %repetitions at T1
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


%% Monte Carlo Methods
[Call, Put] = blsprice(S0,K,r,T2,sigma);
randn('state',0);
[Price, CI] = BermudanPutMC(S0,K,r,T1,T2,sigma,NRepl1,NRepl2);
rand('state',0);
[PriceCond, CICond] = BermudanPutCond(S0,K,r,T1,T2,sigma,NRepl1*NRepl2);
fprintf(1,'European Put = %f\n', Put);
fprintf(1,'MC ->      Price = %f   CI = (%f, %f) \n', ...
    Price, CI(1), CI(2));
fprintf(1,'           Variance = %6.4f%%\n', ...
    100*(CI(2)-CI(1))/Price);
fprintf(1,'MC+Cond -> Price = %f   CI = (%f, %f) \n', ...
    PriceCond, CICond(1), CICond(2));
fprintf(1,'           Variance = %6.4f%%\n', ...
    100*(CICond(2)-CICond(1))/PriceCond);

%% LongstaffSchwartzMC

dt=T/NSteps;
discountVet = exp(-r*dt*(1:NSteps)');
NBasis = length(fhandles); %number of basis functions
alpha = zeros(NBasis,1); %regression parameters
RegrMat = zeros(NRepl,NBasis);
%generate sample paths
SPaths=AssetPaths(S0,r,sigma,T,NSteps,NRepl); %Generates sample paths
SPaths(:,1)=[]; % get rid of starting policies or initial price
CashFlows = max(0, K - SPaths(:,NSteps));
ExerciseTime = NSteps*ones(NRepl,1);
for step = NSteps-1:-1:1
    InMoney = find(SPaths(:,step) < K);
    XData = SPaths(InMoney,step);
    RegrMat = zeros(length(XData), NBasis);
    for k=1:NBasis
        RegrMat(:,k) = feval(fhandles{k}, XData);
    end
    YData = CashFlows(InMoney).*discountVet(ExerciseTime(InMoney)-step);
    alpha = RegrMat \ YData;
    IntrinsicValue = K-XData;
    ContinuationValue = RegrMat*alpha;
    Index = find(IntrinsicValue > ContinuationValue);
    ExercisePaths = InMoney(Index);
    CashFlows(ExercisePaths) = IntrinsicValue(Index);
    ExerciseTime(ExercisePaths) = step;
end % for
price = max(K-S0, mean(CashFlows.*discountVet(ExerciseTime)))


%% Check LS
randn('state',0)
[Call, Put] = blsprice(S0,K,r,T2,sigma);
priceLS = GenericLS(S0,K,r,T,sigma,NSteps,NRepl,fhandles);
[LatS, LatPrice] = binprice(S0,K,r,T,T/NSteps,sigma,0);
priceBin = LatPrice(1,1);
table(priceBin, priceLS)


%% Implicit

dS = Smax/M;  %%getting the rounded stock steps
dt = T/N; %%getting rounded time step

oldval = zeros(M-1,1); %%setting a matrx of zeros one step short to account for matlab not starting at 0
newval = zeros(M-1,1); %%same as above
vetS = linspace(0,Smax,M+1)'; %%setting a linearly spaced matrix with specified values of 0 (min) smax (maximum) M+1 (stock steps+1)
veti = 0:M; %%Matrix from zero to the last stock step
vetj = 0:N; %%Matrix from zero to the last time step

payoff = max(K-vetS(2:M),0); %%payoff is equal to the max difference between strike and the stock price starting in the second row
pastval = payoff; %%past value is the payoff
boundval = K*exp(-r*dt*(N-vetj)); %%boundary value is discounted strike price minus Smin which is zero.

a = 0.5*(r*dt*veti-sigma^2*dt*(veti.^2)); %% 
b = sigma^2*dt*(veti.^2)+r*dt; %% 
c = 0.5*(r*dt*veti+sigma^2*dt*(veti.^2)); %% 

M2 = diag(a(3:M),-1)+diag(1+b(2:M))+diag(c(2:M-1),1); %%M2 the diagonal matrix of alpha beta gamma 

iternmax = 0; 
itern = 0;

aux = zeros(M-1,1); %%matrix of zeros, 1 column wide, M-1 rows long
for j=N:-1:1
    aux(1) = a(2)*(boundval(1,j)+boundval(1,j+1)); %%setting aux matrix as some combination of alpha(2nd input) we use second point because matlab recognises 1 as start not zero 
    
    rhs = pastval(:) + aux; %%RHS as shown in the formula 
    oldval = pastval; 
    error = realmax; 
    
    if mod(j,2)==0 %Used for weekly exercisable, set N = 102
    %if j == round(N/Ex) %use for twice exercisable
            
    while tol<error %%while tollerance is less that the error. i.e. we want the error to be less that tolerance!
        newval(1) = max(payoff(1),oldval(1)+omega/(1+b(2))*(rhs(1)-(1+b(2))*oldval(1)+c(2)*oldval(2))); %Calculating first newval
        for k = 2:M-2
            newval(k)= max(payoff(k),oldval(k)+omega/(1+b(k+1))*(rhs(k)-a(k+1)*newval(k-1)-(1+b(k+1))*oldval(k)+c(k+1)*oldval(k+1))); % Creating a loop to find all the new values and then we take the one where the tolerance criterion is met 
        end
        newval(M-1) = max(payoff(M-1),oldval(M-1)+omega/(1+b(M))*(rhs(M-1)-a(M)*newval(M-2)- (1+b(M))*oldval(M-1)));
        error = norm(newval-oldval); 
        oldval = newval;
        iternmax = iternmax +1;
    end
    pastval = newval;
    
    else
        while tol<error %%while tollerance is less than the error. 
        newval(1) = oldval(1)+omega/(1+b(2))*(rhs(1)-(1+b(2))*oldval(1)+c(2)*oldval(2)); %Calculating again the first newval
        for k = 2:M-2
            newval(k)= oldval(k)+omega/(1+b(k+1))*(rhs(k)-a(k+1)*newval(k-1)-(1+b(k+1))*oldval(k)+c(k+1)*oldval(k+1)); %%Creating a loop to find all the new values and then we take the one where the tolerance criterion is met 
        end
        newval(M-1) = oldval(M-1)+omega/(1+b(M))*(rhs(M-1)-a(M)*newval(M-2)-(1+b(M))*oldval(M-1));
        error = norm(newval-oldval); 
        oldval = newval;
        itern = itern+1;
    end
    pastval = newval; 
        
    end
end
    
    %%changing the pastval from payoff to newval, if the early exercise is better
newval = [boundval(1);newval;0]; %%newval is matrix of all the values of the put option at different nodes, taking into account the possibility of early exercise!

price = interp1(vetS,newval,S0); %%interpolating, when we know all the stock steps and derivative prices, what it will be at the initial price of S0


figure1=figure();
plot(vetS, newval(:,1),'x-');
hold on
title('Bermudan Put Implicit Method')
xlabel('Stock Price')
ylabel('Payoff')
hold off
grid on