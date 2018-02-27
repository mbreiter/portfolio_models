function  x_optimal = MVO_card(mu, Q, targetRet, card)
    % cardinality constraint - need to hold 12
    % no shortselling permitted
    
    n = size(Q,1); 
    
    % asset covariance matrix
    Q = [Q zeros(n); zeros(n,n*2)];
    
    % inequality constraints
    % Upper & lower bound constraints
    A = [mu' zeros(1,n)];
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
    model.rhs = [b; beq];
    model.sense = ['>', '=', '='];
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