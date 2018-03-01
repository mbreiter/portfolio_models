function  x_optimal = MVO_card(mu, Q, targetRet, card)
    % cardinality constraint - need to hold 12
    % no shortselling permitted
    
    n = size(Q,1);
    
    % upper asset bound set arbitarily high
    upper_asset = 1e17;
    % lower asset bound set arbitarily low
    lower_asset = -1e17;
    
    % asset covariance matrix
    Q = [Q zeros(n); zeros(n,n*2)];
    
    % inequality constraints
    % lower bound constraints
    lower = lower_asset * ones(n,1);
    lower_A = [-eye(n) diag(lower)];
    lower_b = zeros(n,1);
    
    % upper bound constraints
    upper = upper_asset * ones(n,1);
    upper_A = [eye(n) -diag(upper)];
    upper_b = zeros(n,1);
    
    A = [lower_A; upper_A; [mu; zeros(n, 1)]'];
    b = targetRet;
    
    % equality constraints
    cardinality_Aeq = [zeros(1,n) ones(1,n)];
    cardinality_beq = card;
    budget_Aeq = [ones(1,n) zeros(1,n)];
    budget_beq = 1;
    Aeq = [cardinality_Aeq; budget_Aeq];
    beq = [cardinality_beq; budget_beq];
    
    % variable types and limits
    var_types = [repmat('C', n, 1); repmat('B', n, 1)];
    lb = zeros(2*n,1);
    ub = ones(2*n,1);
    
    % gurobi model
    clear model;
    model.Q = sparse(Q);
    model.obj = zeros(1,n*2);
    model.A = [sparse(A); sparse(Aeq)];
    model.rhs = full([lower_b;upper_b;b;beq]);
    model.sense = [repmat('<',2*n,1);'>';repmat('=',2,1)];
    model.vtype = var_types;
    model.ub = ub;
    model.lb = lb;
    
    % gurobi parameters
    clear params;
    params.TimeLimit = 100;
    params.OutputFlag = 0;

    results = gurobi(model,params);
    x_optimal = results.x(1:n); 
end