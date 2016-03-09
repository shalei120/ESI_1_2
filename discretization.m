

function X_ast = discretization(X_tilde_ast,K,N)
%% Intialize X by computing R as:
R = zeros(K);
R(:,1) = X_tilde_ast(randi(N),:)';

c = zeros(N,1);
for k = 2:K
    c = c + abs(X_tilde_ast * R(:,k-1));
    cf = find(c == min(c));
    R(:,k) = X_tilde_ast(cf(1),:);
end

fprintf('Step 4 completed...\n')
%% Initialize convergence monitoring parameter phi
phi_bar_ast = 0;

%% Find the optimal discrete solution $X^\ast$ by:
step = 1;
eps = 1e-10;
while 1
    X_tilde = X_tilde_ast * R;
    
    X_ast = zeros(N,K);
    for i = 1:N
        X_ast(i,:) = X_tilde(i,:) == max(X_tilde(i,:));
    end
    
    %% Find the optimal orthonormal matrix $R^\ast$ by:
    [U,omega,U_tilde]=svd(X_ast'*X_tilde_ast);
    phi_bar = trace(omega);
    fprintf('Iteration %d | Difference = %f\n',step,phi_bar - phi_bar_ast)
    if(abs(phi_bar - phi_bar_ast) < eps)
        break;
    end
    phi_bar_ast = phi_bar;
    R = U_tilde * U';
    step = step + 1;
end

end