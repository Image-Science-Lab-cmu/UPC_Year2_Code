function [x_hat_curr, obj] = FISTA(y, beta, psi, psi_adj, A, A_adj, max_iter, init)
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

if exist('init', 'var')
    x_hat_pre = init;
else
    x_hat_pre = A_adj(y);      % x_{-1}
end
x_hat_curr = x_hat_pre;                 % x_{0}

gradient = @(x) -A_adj(y-A(x));
least_square = @(x) 0.5*sum(vec(y-A(x)).^2);

t = 5e-1; %%starting step size
% t_beta = 0.5; %%shrinkage per iteration

obj_old = 1e100;
for iter = 1:max_iter
    
%     x_hat_pre_pre = x_hat_pre;  % x_{k-2}
    x_hat_pre = x_hat_curr;     % x_{k-1}
    
    % update x_hat
%     v = x_hat_pre + (iter-2)/(iter+1) * (x_hat_pre - x_hat_pre_pre);
    v = x_hat_pre;
    g = gradient(v);  
%     x_hat_curr = psi_adj(soft_thresholding(psi(v-t*g), t*beta));
    x_hat_curr = psi_adj(psi(v-t*g));
    figure(1), imshow(reshape(x_hat_curr, 1024+255, 2048+255),[]);
%    imshow(x_hat); drawnow
    
    % print status
    obj = 0.5 * norm(vec(y-A(x_hat_curr)))^2 + beta * norm(vec(psi(x_hat_curr)),1);
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
