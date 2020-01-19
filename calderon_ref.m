%% Y esto tambien

% Calderon precondiioner with refined mesh

disp(' ')

N = size(Ze,1);

% Refining the mesh
[obj2, trian_ch, trian_pr] = refine_mesh6(obj);
obj2 = get_edge(obj2);
obj2.N = length(obj2.ln);
N2 = obj2.N;
%user_plot_geom3d(obj2); 
%plot_obj_numbers(obj2);

% Obtaining RWG as lc of refined RWG
lc_mat = linear_comb6(obj, obj2, trian_ch);
R = lc_mat.';
piR = pinv(R);

G2 = nxD2_mat(obj, obj2, trian_pr);
piG2 = pinv(G2);

% Obtaining BC as lc of refined RWG
lc_mat_bc = linear_comb_bc(obj, obj2, trian_ch);
P = lc_mat_bc.';
piP = pinv(P);

EM_data.field = 1;
Ze2 = user_impedance(1:N2, 1:N2, obj2, EM_data); 

%T2jmr = P.'*Ze2 * P*inv(Gm)*R.'*Ze2*R; % Andriulli
%T2jmr =    Ze2*piG2*Ze;
%T2jmr = R'*Ze2*piG2*Ze;
%T2jmr = P'*Ze2*(G2\Ze);
T2jmr = P'*Ze2*piG2*Ze;
fprintf('cond T           = %.2e\n', cond(Ze));
fprintf('cond T2 JMR      = %.2e\n', cond(T2jmr));

Mjmr = P'*Ze2*piG2;
%Mjmr = Ze2*piG2;

%Jjmr = -T2jmr\(Mjmr*Ei);
Jjmr = -T2jmr\(P'*Ze2*(G2\Ei));


%% Metodo Andriulli
G = nxD_mat(obj2);
Gm = R.'*G*P;
Q = P*inv(Gm)*R.';
T2A = P.'*Ze2*Q*Ze2*R;
fprintf('cond T2 Andriulli = %.2e\n\n', cond(T2A));

%% Andriulli method via orthogonal projections

Gb = D_mat(obj2); % Gram matrix for refined RWG functions
nxGb = nxD_mat(obj2); % Gram matrix for refined nxRWG functions with refined RWG functions
Gbc = P.'*Gb*P; % Gram matrix for BC functions
G1 = D_mat(obj); % Gram matrix for original RWG (and nxRWG functions)

Tand = P.'*Ze2*P*inv(inv(G1)*Gm)*inv(G1)*R.'*Ze2*R; % Developed Andriulli Matrix
Tand2 = P.'*Ze2*P*inv(Gbc)*Gm*inv(G1)*R.'*Ze2*R; % Developed Pablo Matrix
Tand3 = P.'*Ze2*P*inv(Gbc)*Gm'*inv(G1)*R.'*Ze2*R; % Refined Pablo Matrix
fprintf('cond Gb = %.2e\n', cond(Gbc));
fprintf('cond Tand = %.2e\n', cond(Tand));
fprintf('cond Tand2 = %.2e\n', cond(Tand2));
fprintf('cond Tand3 = %.2e\n\n', cond(Tand3));

comp_mat(inv(inv(G1)*Gm), inv(Gbc)*Gm', 'Comparison between projection matrices');

%% Definition of operators

G3 = nxD_mat(obj);

TPrwg = inv(G1)*R.'*Ze2*R; % T operator + projection on nxRWG
TPbc = inv(Gbc)*P.'*Ze2*P; % T operator + projection on nxBC
PIbc = inv(Gbc)*Gm'; % Projection from nxRWG on BC
PInrwg = inv(G1)*Gm; % Projection from BC to nxRWG
PIrwg = inv(G1)*G3; % Projection from nxRWG to RWG
PIrwgb = inv(Gb)*G2'; % Projection from nxRWG to RWGb

fprintf('-----------------------\nCONDITION NUMBERS FOR DIFFERENT PROJECTION OPERATORS\n')

fprintf('Condition number for projection from nxRWG to RWG: %.5e\n', cond(PIrwg));
fprintf('Condition number for projection from nxRWG to BC: %.5e\n', cond(PIbc));
fprintf('Condition number for projection from nxRWG to RWGb: %.5e\n', cond(PIrwgb));
fprintf('Condition number for (inverse of) projection from BC to nxRWG: %.5e\n\n', cond(inv(PInrwg)));

%% Methods using operators and condition numbers

Tjmr1 = TPrwg*PIrwg*TPrwg;
Tjmr2 = TPbc*PIbc*TPrwg;
Tjmr3 = inv(G1)*R.'*Ze2*PIrwgb*TPrwg;
Tandr = TPbc*inv(PInrwg)*TPrwg;

fprintf('Condition number for JMR projection method 1: %.5e\n', cond(Tjmr1));
fprintf('Condition number for JMR projection method 2: %.5e\n', cond(Tjmr2));
fprintf('Condition number for JMR projection method 3: %.5e\n', cond(Tjmr3));
fprintf('Condition number for andriulli projection method: %.5e\n\n', cond(Tandr));

%% Checking pseudoinverse conditions

% mypiG2 = mypseudo(G2, 0.00000000001);
% 
% fprintf('Checking pseudoinverse conditions:\n\n');
% 
% comp_mat(G2*piG2*G2, G2, 'Error in G2*piG2*G2 = G2');
% comp_mat(piG2*G2*piG2, piG2, 'Error in piG2*G2*piG2 = piG2');
% comp_mat((G2*piG2)', G2*piG2, 'Error in (G2*piG2) = G2*piG2');
% comp_mat((piG2*G2)', piG2*G2, 'Error in (piG2*G2) = piG2*G2');
% 
% fprintf('\n');
% 
% fprintf('Checking pseudoinverse conditions for user defined pseudoinverse:\n\n');
% 
% comp_mat(G2*mypiG2*G2, G2, 'Error in G2*mypiG2*G2 = G2');
% comp_mat(mypiG2*G2*mypiG2, mypiG2, 'Error in mypiG2*G2*mypiG2 = mypiG2');
% comp_mat((G2*mypiG2)', G2*mypiG2, 'Error in (G2*mypiG2) = G2*mypiG2');
% comp_mat((mypiG2*G2)', mypiG2*G2, 'Error in (mypiG2*G2) = mypiG2*G2');
% 
% fprintf('\n');
% 
% fprintf('Checking pseudoinverse conditions for P*inv(Gm):\n\n');
% 
% comp_mat(G2*P*inv(Gm)*G2, G2, 'Error in G2*P*inv(Gm)*G2 = G2');
% comp_mat(P*inv(Gm)*G2*P*inv(Gm), P*inv(Gm), 'Error in P*inv(Gm)*G2*P*inv(Gm) = P*inv(Gm)');
% comp_mat((G2*P*inv(Gm))', G2*P*inv(Gm), 'Error in (G2*P*inv(Gm)) = G2*P*inv(Gm)');
% comp_mat((P*inv(Gm)*G2)', P*inv(Gm)*G2, 'Error in (P*inv(Gm)*G2) = P*inv(Gm)*G2');
% 
% fprintf('\n');
% 
% comp_mat(piG2, P*inv(Gm), 'Error in piG2 = P*inv(Gm)');
% comp_mat(piG2, mypiG2, 'Error in piG2 = mypiG2');

MA = P.'*Ze2*P*inv(Gm);

JA = -T2A\(MA*Ei);


%% Errores de cálculo directo
fprintf('\n');
comp_mat(Jjmr, Je, 'error in Jjmr (direct)     ');
comp_mat(JA, Je,   'error in JA   (direct)     ');

%% Iterative solution

% astr = {'bicg','bicgstab','cgs','gmres','lsqr','qmr','tfqmr'};

% Already defined in run_3d
%method = 'lsqr';
%maxit = ceil(0.5*N);
%tol = 1e-3;

fprintf('\nIterative solution JMR, N = %d, tol = %g, maxit = %d\n', N, tol, maxit);
tic; [Jjmrit, flag_jmr, relres_jmr, iter_jmr, resvec_jmr] = iterative_sol(method, T2jmr, -Mjmr*Ei, tol, maxit); toc
%tic; [Jjmrit, flag_jmr, relres_jmr, iter_jmr, resvec_jmr] = iterative_sol(method, T2jmr, -P'*Ze2*(G2\Ei), tol, maxit); toc
fprintf('Nit = %d, error:            %.3e\n\n', prod(iter_jmr), norm(Jjmrit-Je)/norm(Je) );

fprintf('\nIterative solution Andriulli, N = %d, tol = %g, maxit = %d\n', N, tol, maxit);
tic; [JAit, flag_A, relres_A, iter_A, resvec_A] = iterative_sol(method, T2A, -MA*Ei, tol, maxit); toc
fprintf('Nit = %d, error:            %.3e\n\n', prod(iter_A), norm(JAit-Je)/norm(Je) );

figure;
semilogy(resvec_jmr/resvec_jmr(1),'linewidth',1);
hold on
semilogy(resvec_A/resvec_A(1),'linewidth',1);
hold off
grid
title(sprintf('Iterative solution, N = %d, tol = %g, maxit = %d\n', N, tol, maxit));
legend('JMR','Andriulli')
ylabel('Relative error in residual')
xlabel('Iteration')



%% Comprobaciones

% err_jmr = comp_mat(T2jmr, T2A, 'error in T2jmr vs T2A      ');
% fprintf('\n');
% comp_mat(piG2, P*inv(Gm), 'error in G2+ == P*inv(Gm)  ');
% comp_mat(G2*P*inv(Gm), eye(N), 'error in G2*P*inv(Gm) == I ');
% comp_mat(G2*piG2,eye(N), 'error in G2*G2+ == I       ');
% 
% % Gm == nxD2*P, (irrelevant, works for any P)
% errDb = comp_mat(G2*P, Gm, 'error in Gm == G2com*P     ');
% 
% % nxD2 = R.'*G
% errnxD2 = comp_mat(R.'*G, G2, 'error in G2 == R.''*G       ');

%% Norm of the error vectors in JM and Andriulli methods: V1

fprintf('------------------------\nFIRST VERSION: EXTENDED BOTH BASES \n')

% Define the Gram matrix for the super vector space
Gram = [Gb, nxGb'; nxGb, Gb];

% Define the original vector and projections for every approximation
Vo = R*TPrwg*Je;
Vand = P*inv(PInrwg)*TPrwg*Je;
Vjm1 = P*PIbc*TPrwg*Je;
Vjm2 = R*PIrwg*TPrwg*Je;
Vjm3 = PIrwgb*TPrwg*Je;

% Vectors on direct sum superspace and error norm calculation

Eand = [- Vand; Vo];
Ejm1 = [- Vjm1; Vo];
Ejm2 = [- Vjm2; Vo];
Ejm3 = [- Vjm3; Vo];

fprintf('Error with Andriulli approximation: %.5e,\n', Eand'*Gram*Eand);
fprintf('Error with JMR approximation 1 on RWG: %.5e,\n', Ejm1'*Gram*Ejm1);
fprintf('Error with JMR approximation 2 on RWG: %.5e,\n', Ejm2'*Gram*Ejm2);
fprintf('Error with JMR approximation 3 on RWG: %.5e.\n\n', Ejm3'*Gram*Ejm3);

%% Norm of the error vectors in JM and Andriulli methods: V2

fprintf('------------------------\nFIRST VERSION: EXTENDED ONLY DOMAIN BASE \n')

% Define the Gram matrix for the super vector space
Gram = [Gb, G2'; G2, G1];

% Define the original vector and projections for every approximation
Vo = TPrwg*Je;
Vand = P*inv(PInrwg)*TPrwg*Je;
Vjm1 = P*PIbc*TPrwg*Je;
Vjm2 = R*PIrwg*TPrwg*Je;
Vjm3 = PIrwgb*TPrwg*Je;

% Vectors on direct sum superspace and error norm calculation

Eand = [- Vand; Vo];
Ejm1 = [- Vjm1; Vo];
Ejm2 = [- Vjm2; Vo];
Ejm3 = [- Vjm3; Vo];

fprintf('Error with Andriulli approximation: %.5e,\n', Eand'*Gram*Eand);
fprintf('Error with JMR approximation 1 on RWG: %.5e,\n', Ejm1'*Gram*Ejm1);
fprintf('Error with JMR approximation 2 on RWG: %.5e,\n', Ejm2'*Gram*Ejm2);
fprintf('Error with JMR approximation 3 on RWG: %.5e.\n', Ejm3'*Gram*Ejm3);
