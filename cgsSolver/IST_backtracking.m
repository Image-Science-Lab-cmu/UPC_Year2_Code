function [x_hat, obj] = IST_backtracking(y, beta, psi, psi_adj, A, A_adj, max_iter)
%Solves \min \| y - A(x) \|^2 +beta \| Psi(x) \|_1 
%
% Input
%   y -- measurements
%   A -- function handle of A
%   A_adj -- function handle of A adjoint
%   psi -- function handle of psi
%   psi_adj -- function handle op psi adjoint
%   beta --- scaling parametere
%   max_itere --- maximum iterations
%
% Output
%   x -- solution to optimization problem

x_hat = zeros(size(A_adj(y)));

gradient = @(x) -A_adj(y-A(x));
least_square = @(x) 0.5*sum(vec(y-A(x)).^2);

t0 = 2; %%starting step size
t_beta = 0.5; %%shrinkage per iteration

obj_old = 1e100;
for iter = 1:max_iter
    t = t0;
    g = gradient(x_hat);
    while true
        Gt = (x_hat - psi_adj(soft_thresholding(psi(x_hat-t*g), t*beta)))/t;
        Gt = (x_hat - (soft_thresholding((x_hat-t*g), t*beta)))/t;
        test_lhs = least_square(x_hat - t*Gt);
        test_rhs = least_square(x_hat) - t*vec(g)'*vec(Gt) + 0.5*t*sum(vec(Gt).^2);
        if test_lhs > test_rhs
            t = t_beta * t;
        else
            break;
        end
    end
    % update x_hat 
%     x_hat = psi_adj(soft_thresholding(psi(x_hat-t*g), t*beta));
    x_hat = (soft_thresholding((x_hat-t*g), t*beta));
%    imshow(x_hat); drawnow
    
    % print status
%     obj = 0.5 * norm(vec(y-A(x_hat)))^2 + beta * norm(vec(psi(x_hat)),1);
    obj = 0.5 * norm(vec(y-A(x_hat)))^2 + beta * norm(vec((x_hat)),1);
    fprintf('iter = %d obj = %f\n', iter, obj);
    if (2*abs(obj_old - obj)/(1e-8+obj+obj_old)) < 1e-6
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
