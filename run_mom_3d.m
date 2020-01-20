% MoM 3D
%
% Juan M. Rius et al., AntennaLab, Universitat Politecnica de Catalunya (Spain), 
% v1        2007
% v1.1.2    September 2015
% Isolated from ACA solver in October 2018

%% User configurable parameters 

% Incident field
th_i = 180;
ph_i = 0;
rot_EH = 90;

% Electromagnetic
lambda = 1/200;
k = 2*pi/lambda; 
eta = 120*pi;
field = 1; % EFIE->1, MFIE->2, CFIE->3

%
% Rint_s = 0.2;       % MoM Integration radius (meters). Rint=0 is enough if basis functions are very small.
% Rint_f = Rint_s;
% Ranal_s = 0;
Rint_s = 1;       % MoM Integration radius (meters). Rint=0 is enough if basis functions are very small.
Rint_f = Rint_s;
Ranal_s = 1;
corr_solid = 0;
flag = 0;

EM_data = struct('lambda',lambda, 'k',k, 'eta',eta, 'field',field, 'Rint_s',Rint_s, 'Rint_f',Rint_f, 'Ranal_s',Ranal_s, 'corr_solid',corr_solid, 'flag',flag );

%% Object geometry

geom = 'sphere';
Ne = 192*4; 
radio = 2e-4; 
param = struct('R',radio, 'Ne',Ne);

% geom = 'cube';
% Nt = 192; 
% L = 2e-6; 
% param = struct('L',L, 'Nt',Nt);

% geom = 'tetrahedron';
% param = struct('R',1, 'Ne',192);

% geom = 'sq_plate';
% param = struct('Nx',10, 'Lx',1, 'Nz',10, 'Lz',1, 'x',0, 'y',0, 'z',0, 'cor',0);

% geom = 'sq_plate';
% param = struct('Nx',1, 'Lx',1, 'Nz',2, 'Lz',2, 'x',0.5, 'y',0, 'z',1, 'cor',0);

% Compute object geometry data
cd objects
obj = feval(geom,param);
cd ..

obj = get_edge(obj);

%obj = get_edge_jmr(obj);
obj.N = length(obj.ln); N = obj.N;
obj.name = geom;

% Display object geometry with feeding in red
h_obj3d = user_plot_geom3d(obj);
%plot_obj_numbers(obj); 
drawnow;

%% MoM matrices
[Ei, Hi] = test_fields(obj, th_i, ph_i, rot_EH, k, eta);
D = D_mat(obj);

EM_data.field = 1;
Ze = user_impedance(1:N, 1:N, obj, EM_data); 
disp('Direct solution EFIE:')
tic; Je = -Ze\Ei; toc;

EM_data.field = 2;
Zh = user_impedance(1:N, 1:N, obj, EM_data); 
disp('Direct solution MFIE:')
tic; Jh = -(Zh-D)\Hi; toc;  % With flag == 0
%tic; Jh = -Zh\Hi; toc;   % With flag == 1

%% Iterative solution

% astr = {'bicg','bicgstab','cgs','gmres','lsqr','qmr','tfqmr'};
%method = 'lsqr';
method = 'gmres';
maxit = N; %ceil(0.5*N);
tol = 1e-4;

fprintf('\nIterative solution EFIE, N = %d, tol = %g, maxit = %d\n', N, tol, maxit);
tic; [Jeit, flag_e, relres_e, iter_e, resvec_e] = iterative_sol(method, Ze, -Ei, tol, maxit); toc
fprintf('Nit = %d, error J it vs direct = %.3e\n\n', prod(iter_e), norm(Jeit-Je)/norm(Je) );

fprintf('\nIterative solution MFIE, N = %d, tol = %g, maxit = %d\n', N, tol, maxit);
tic; [Jhit, flag_h, relres_h, iter_h, resvec_h] = iterative_sol(method, Zh-D, -Hi, tol, maxit); toc
fprintf('Nit = %d, error J it vs direct = %.3e\n\n', prod(iter_h), norm(Jhit-Jh)/norm(Jh) );

figure;
semilogy(resvec_e/resvec_e(1),'linewidth',1);
hold on
semilogy(resvec_h/resvec_h(1),'linewidth',1);
hold off
grid
title(sprintf('Iterative solution, N = %d, tol = %g, maxit = %d\n', N, tol, maxit));
legend('EFIE','MFIE')
ylabel('Relative error in residual')
xlabel('Iteration')

%% RCS
ang_e1 = 1; ang_e2 = 179; ang_h1 = 1; ang_h2 = 179; M_e = 179; M_h = 179;
[RCSbe_Je, RCSbh_Je] = bist_rcs(Je, obj, k, eta, th_i, ph_i, rot_EH, ang_e1, ang_e2, ang_h1, ang_h2, M_e, M_h);
[RCSbe_Jh, RCSbh_Jh] = bist_rcs(Jh, obj, k, eta, th_i, ph_i, rot_EH, ang_e1, ang_e2, ang_h1, ang_h2, M_e, M_h);

%% Post-processing
user_plot_obj_current(Je,obj, EM_data);
title('Je computed original')
%user_plot_radpat3d(Je, obj, EM_data)

if strcmp(geom,'sphere')
    % Exact results for sphere
    run_sphe; 
    comp_mat(Je, Jexact, 'Error J EFIE');
    comp_mat(Jh, Jexact, 'Error J MFIE');    
else
    ang_e = linspace(ang_e1,ang_e2,M_e);
    ang_h = linspace(ang_h1,ang_h2,M_h);
    
    figure; 
    plot(ang_e, 10*log10([RCSbe_Je RCSbe_Jh]))
    title('Bistatic RCS E-plane')
    legend('EFIE','MFIE')
    
    figure; 
    plot(ang_h, 10*log10([RCSbh_Je RCSbh_Jh]))
    title('Bistatic RCS H-plane')
    legend('EFIE','MFIE')
end    


