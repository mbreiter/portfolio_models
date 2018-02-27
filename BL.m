function  x_optimal = BL(Q, tau, P, q, Omega, lambda, mktNoShares, currentPrices)

    % Find the total number of assets
    n = size(Q,1); 

    % Market portfolio
    x_mkt = mktNoShares .* currentPrices ./ sum(mktNoShares .* currentPrices);
    
    % Market equilibrium returns
    pi = lambda * Q * x_mkt;
    
    % Calculate the updated expected returns from our views
    mu_view = inv(inv(tau*Q) + P'*inv(Omega)*P) * ( inv(tau*Q)*pi + P'*inv(Omega)*q);

    % We incorporate our updated expected returns into our objective function 
    f = -mu_view;
    
    % We do not have any inequality constraints. 
    A = [];
    b = [];
    
    % Aeq and beq will house our equality constraint, ensuring that the
    % proprotions in our portifolio sum to 1
    Aeq = ones(1,n);
    beq = 1;  
    
    % We allow shortselling and otherwise have no bounding constraints
    ub = [];
    lb = [];
   
    % We use quadprog to optimize
    x_optimal = quadprog(lambda*Q,f,A,b,Aeq,beq,lb,ub);
        
end