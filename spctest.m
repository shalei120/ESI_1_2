% I = imread('1.bmp');
% I=rgb2gray(I);
% I = reshape(I,93*93,1);
% W = ones(93*93);
% for i = 1:93*93
%     for j = i+1:93*93
%         W(i,j) = I(i) == I(j);
%         W(j,i) = W(i,j);
%     end
% end
W=WTC;
W = W + 0.1;

fprintf('Step 0 completed...\n')

N = size(W,1);
K = 5;    

D=diag(W*ones(N,1));

fprintf('Step 1 completed...\n')

kernel = D^(-0.5)*W*D^(-0.5);
[V,S]=eigs(kernel);
Z = D^(-0.5)*V(:,1:K);
fprintf('Step 2 completed...\n')

X_tilde_ast = ((Z*Z').*eye(N))^(-0.5)*Z;
fprintf('Step 3 completed...\n')

R = zeros(K);
R(:,1) = X_tilde_ast(randi(N),:)';

c = zeros(N,1);
for k = 2:K
    c = c + abs(X_tilde_ast * R(:,k-1));
    cf = find(c == min(c));
    R(:,k) = X_tilde_ast(cf(1),:);
end

fprintf('Step 4 completed...\n')

phi_bar_ast = 0;

step = 1;
eps = 1e-10;
while 1

    X_tilde = X_tilde_ast * R;

    X_ast = zeros(N,K);
    for i = 1:N
        X_ast(i,:) = X_tilde(i,:) == max(X_tilde(i,:));
    end

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
% I = X_ast * [0 64 128 192 256]';
% I = reshape(I,93,93);
% I = uint8(I);
% imshow(I)
