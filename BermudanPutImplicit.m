function price = BermudanPutImplicit(S0,K,r,T,sigma,Smax,M,N,omega,tol,Ex)

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
