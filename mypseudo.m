%% Pseudoinverse script

% IN: An m x n matrix A.
% OUT: An n x m matrix piA, the pseudoinverse of A.

%% Pseudoinverse function declaration and preliminaries

function[piA] = mypseudo(A, tol)
    
    % Size declaration 
    [m, n] = size(A);
    
    % Conjugate transpose of A
    ctA = A';
    
    % Declaration of the initial matrix according to Tychonoff's Regularization Theorem
    r = 0.0001;
    piA = inv(ctA*A + r*eye(n))*ctA;
    
    % Checking pseudoinverse conditions
    while (nw_comp_mat(A*piA*A, A) > tol) || (nw_comp_mat(piA*A*piA, piA) > tol) || (nw_comp_mat((A*piA)', A*piA) > tol) || (nw_comp_mat((piA*A)', piA*A) > tol)
        
        % Pseudoinverse actualization
        piA = 2*piA - piA*A*piA;
        
    end

end