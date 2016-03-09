%% Calculate f(X_t) or f(X_s)

function v = f(X,J)
v = trace((J'*X)*(J'*X)');
end