% err = comp_mat(A,B,'st')
%
% ||A-B||/||B||
% Prints 'st: err'

function err = comp_mat(A, B, st)

err = sqrt(sum(sum(abs(A-B).^2)))/sqrt(sum(sum(abs(B).^2)));
fprintf('%s: %.2e\n', st, err);