% err = comp_mat(A,B,'st')
%
% ||A-B||/||B||
% Prints 'st: err'

function err = nw_comp_mat(A, B)

    err = sqrt(sum(sum(abs(A-B).^2)))/sqrt(sum(sum(abs(B).^2)));

end