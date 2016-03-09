%% The main program of muticlass spectral clustering
function X_ast = SpectralClustering(W,D,K)

%% Basic parameters
N = size(W,1); % the number of templates

%% Compute the degree matrix
%D=diag(W*ones(N,1));
%W = W + 0.1;
fprintf('Step 1 completed...\n')
%% Find the optimal eigensolution $Z^\ast$ by:
kernel = D^(-0.5)*W*D^(-0.5);
%opts.tol = 1e-3; 
[V,S]=eigs(kernel);
Z = D^(-0.5)*V(:,1:K);
fprintf('Step 2 completed...\n')

%% Normalize $Z^\ast$ by:
X_tilde_ast = ((Z*Z').*eye(N))^(-0.5)*Z;
fprintf('Step 3 completed...\n')
X_ast = X_tilde_ast;
end