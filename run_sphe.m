% Script for computation of exact results for sphere
%
% Version 3.1a: Incluye calculo de la J ex

if ~strcmp(geom,'sphere')
	error('Exact solution available only for sphere');
else
	if th_i==180 && ph_i==0 && rot_EH==90
		disp('Sphere: Computation of current');tic;
      	Jexact = j_ex(param.R,D,obj,k,eta);
    else
		disp('Warning:');
		disp('Exact current computation only available for incidence direction');
		disp('theta=180, phi=0 and rotation of E with respect to phi=90');
		Jexact = [];
    end

    if isempty(ang_e1) || isempty(ang_e2) || isempty(ang_h1) || isempty(ang_h2) || isempty (M_e) || isempty (M_h)
        error('All biestatic parameters must be specified');
    end
    disp('Sphere bistatic RCS:');
    ang_e = linspace(ang_e1,ang_e2,M_e);
    ang_h = linspace(ang_h1,ang_h2,M_h);
    [RCSbe_ex, RCSbh_ex] = sph_bist(param.R,lambda,ang_e,ang_h);


%     err_RCSbe = norm(RCSbe - RCSbe_ex)/norm(RCSbe_ex);
%     fprintf('Error RCS E-plane = %.2e\n', err_RCSbe)

    figure; 
    plot(ang_e, 10*log10([RCSbe_ex RCSbe_Je RCSbe_Jh]))
    title('Bistatic RCS E-plane')
    legend('Exact','EFIE','MFIE')
    
    figure; 
    plot(ang_h, 10*log10([RCSbh_ex RCSbh_Je RCSbh_Jh]))
    title('Bistatic RCS H-plane')
    legend('Exact','EFIE','MFIE')

end