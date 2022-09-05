function [x_hat, obj] = ADMM(y, lambda, rho, ...
    Dforward, DAdj, A, A_adj, max_iter)
% Solves \min \| y - A(x) \|^2 + lambda \| D(x) \|_1 
% using ADMM algorithm
% We reformulate this problem as 
% \min f(x) + g(Gx), 
% where f(x)=0, g(Gx) = \| y - A(x) \|^2 + lambda \| D(x) \|_1.
% We use auxilary variable u1 = A(x), u2 = D(x), w1, w2
% \min \| y - u1 \|^2 + lambda \| u2 \|_1 + ...
% \rho / 2 * (\| A(x) - u1 + w_1 \|^2) + ...
% \rho / 2 * (\| D(x) - u2 + w_2 \|^2) + ...
% - \rho/2 * \| w_1 \|^2 + \rho/2 * \| w_1 \|^2.
%
% Input
%   y -- measurements
%   A -- function handle of A
%   A_adj -- function handle of A adjoint
%   Dforward -- function handle of differential in xy-direction
%   DAdj -- function handle of Dforward adjoint
%   lambda --- coefficients of TV norm
%   rho --- coefficient of scaled form ADMM
%   max_iter --- maximum iterations
%
% Output
%   x      -- solution to optimization problem
%   u1     -- augmented variable u1 = A(x)
%   u2     -- augmented variable u2 = Dforward(x)
%   w1     -- scaled lagrange variable of u1
%   w2     -- scaled lagrange variable of u2
%
% Algorithm at k-step:
%   1. update x_{k+1}:
%      x_{k+1} = argmin_x f(x) + |A(x)-u1^{k}+w1^{k}|^2 + |D(x)-u2^{k}+w2^{k}|^2
%      Take the derivative of it wrt. x and set to zero.
%      (A^T*A + D^T*D) x = A^T*(u1^{k}-w1^{k}) + D^T*(u2^{k}-w2^{k});
%
%   2. update u1_{k+1}
%      u1_{k+1} = argmin_u1 0.5*|y - u1|^2 + \rho/2*|A(x_{k+1}) - u1 + w1_{k}|^2
%      Take the derivative of it wrt. u1 and set to zero.
%      u1_{k+1} = 1/(1+\rho)*(y + \rho*(A(x_{k+1}) + w1_{k}));
%
%      update u2_{k+1}
%      u2_{k+1} = argmin_u2 \lambda |u2|_1 + \rho/2 * |D(x_{k+1}) - u2 + w2_{k}|^2
%      It has a close-form solution
%      u2_{k+1} = S_{\lambda/\rho} (D(x_{k+1}) + w2_{k});
%
%   3. update w1_{k+1} and w2_{k+1}
%      w1_{k+1} = w1_{k} + (A(x_{k+1}) - u1_{k+1});
%      w2_{k+1} = w2_{k} + (D(x_{k+1}) - u2_{k+1}).


x_hat = zeros(size(A_adj(y)));
u1_hat = zeros(size(A(x_hat)));
u2_hat = zeros(size(Dforward(x_hat)));
w1_hat = zeros(size(u1_hat));
w2_hat = zeros(size(u2_hat));

if ~exist('rho', 'var')
rho = 1; % todo
% rho = 1e-1;
end

% gradient = @(x, alpha, w) -A_adj(y-A(x)) + rho*(x-alpha-w);
% least_square = @(x) 0.5*sum(vec(y-A(x)).^2);

% t0 = 2; %%starting step size
% t_beta = 0.5; %%shrinkage per iteration

% figure, hold on
obj_old = 1e100;
for iter = 1:max_iter
%     imshow(reshape(x_hat, [1764, 2220]));
    % update x_hat (this operation is very slow)
    A_x = @(x) (A_adj(A(x))+DAdj(Dforward(x)));
    b_x = gpuArray(A_adj(u1_hat - w1_hat) + DAdj(u2_hat - w2_hat));
    x_hat = cgs(A_x, b_x, 1e-6, 100, [], [], vec(x_hat));
    
%     H = 1024; W = 2048; kH=256; kW=256;
%     figure(1), imshow(reshape(x_hat, [H+kH-1, W+kW-1]), []);
    
    % update auxilary variable u1 and u2
    u1_hat = 1/(1+rho) * (y + rho * (A(x_hat) + w1_hat));
    u2_hat = soft_thresholding(Dforward(x_hat) + w2_hat, lambda/rho);
    
    % update w_hat
    w1_hat = w1_hat + A(x_hat) - u1_hat;
    w2_hat = w2_hat + Dforward(x_hat) - u2_hat;
    
    % print status
    obj = 0.5 * norm(vec(y-A(x_hat)))^2 + ...
        lambda * norm(vec(Dforward(x_hat)),1);
    fprintf('iter = %d obj = %f\n', iter, obj);
    if (2*abs(obj_old - obj)/(1e-8+obj+obj_old)) < 1e-6
        break;
    end
    
    if iter >=20 && obj > 1000
        break;
    end
    obj_old = obj;
end
end

function y = vec(x)
y = x(:);
end

function y = soft_thresholding(x, beta)
y = max(0, x - beta) - max(0, -x-beta);
end
