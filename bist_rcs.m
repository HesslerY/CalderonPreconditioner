% Script to compute bistatic RCS

function [RCSbe, RCSbh] = bist_rcs(J, obj, k, eta, theta_i, phi_i, rot_EH_i, ang_e1, ang_e2, ang_h1, ang_h2, M_e, M_h)

if isempty(ang_e1) || isempty(ang_e2) || isempty(ang_h1) || isempty(ang_h2) || isempty(M_e) || isempty(M_h)
	error('All bistatic parameters must be specified');
end

disp('Computing bistatic RCS:');

th = theta_i*pi/180; ph = phi_i*pi/180; rot_EH = rot_EH_i*pi/180;
k_i = -[sin(th)*cos(ph); sin(th)*sin(ph); cos(th)];
e_i = cos(rot_EH)*[-sin(ph); cos(ph); 0] + sin(rot_EH)*[cos(th)*cos(ph); cos(th)*sin(ph); -sin(th)];
h_i = cross(k_i,e_i);

ang_e = linspace(ang_e1,ang_e2,M_e);
ang_h = linspace(ang_h1,ang_h2,M_h);

r_e = e_i*sin(pi*ang_e/180) - k_i*cos(pi*ang_e/180);
r_h = h_i*sin(pi*ang_h/180) - k_i*cos(pi*ang_h/180);

RCSbe = zeros(M_e,1);
for m=1:M_e
	RCSbe(m) = rcs(obj, J, k, eta, r_e(:,m));
end

RCSbh = zeros(M_h,1);
for m=1:M_h
	RCSbh(m) = rcs(obj, J, k, eta, r_h(:,m));
end


