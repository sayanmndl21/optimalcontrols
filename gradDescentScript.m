clear; clc;
f = @(x) (x(1)-4).^2 + (x(1)-4).*(x(2)-2) + 3*(x(2)-2).^2 +3;
[xopt,fopt,niter,gnorm,dx]=grad_descent([3;3],f,@gradf,0.1,1000)

function g = gradf(x)
g = [2*(x(1)-4) + (x(2)-2)
    (x(1)-4) + 6*(x(2)-2)];
end