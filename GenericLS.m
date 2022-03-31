function price = GenericLS(S0,K,r,T,sigma,NSteps,NRepl,fhandles)
randn('state',0)
dt=T/NSteps;
discountVet = exp(-r*dt*(1:NSteps)');
NBasis = length(fhandles); %number of basis functions
alpha = zeros(NBasis,1); %regression parameters
RegrMat = zeros(NRepl,NBasis);
%generate sample paths
SPaths=AssetPaths(S0,r,sigma,T,NSteps,NRepl); %Generates sample paths
SPaths(:,1)=[]; % get rid of starting policies or inital price
%
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
price = max(K-S0, mean(CashFlows.*discountVet(ExerciseTime)));



