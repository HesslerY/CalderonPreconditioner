function [J, flag, relres, iter, resvec] = iterative_sol(method, A, b, tol, maxit, M)
% Wrapper for calling iterative solver
GMRES_restart = 15;

if strcmp(method,'gmres')
    if nargin == 6
       [J, flag, relres, iter, resvec] = feval(method, A, b, GMRES_restart, tol, maxit, M);
    else
        [J, flag, relres, iter, resvec] = feval(method, A, b, GMRES_restart, tol, maxit);
    end
else
    if nargin == 6
        [J, flag, relres, iter, resvec] = feval(method, A, b, tol, maxit, M);
    else
        [J, flag, relres, iter, resvec] = feval(method, A, b, tol, maxit);
    end
end

