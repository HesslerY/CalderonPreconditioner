% [Ei, Hi] = test_fields(obj, theta_i, phi_i, rot_EH_i, k)
% User function to compute incident field (linear system independent vector)
%
% Juan M. Rius, AntennaLab, Universitat Politecnica de Catalunya (Spain), v1.1, July 2008

function [Ei, Hi] = test_fields(obj, theta_i, phi_i, rot_EH_i, k, eta)

th = theta_i*pi/180; ph = phi_i*pi/180; rot_EH = rot_EH_i*pi/180;

k_i = -[sin(th)*cos(ph); sin(th)*sin(ph); cos(th)];
e_i = cos(rot_EH)*[-sin(ph); cos(ph); 0] + sin(rot_EH)*[cos(th)*cos(ph); cos(th)*sin(ph); -sin(th)];
h_i = cross(k_i,e_i);


Tp = obj.edges(1,:); Tm = obj.edges(2,:);	% T+ and T- triangles corresponding to vertex
pp = exp(-1j*k* k_i.'*obj.cent(:,Tp));	% Plane wave at cent of T+, 1xNe
pm = exp(-1j*k* k_i.'*obj.cent(:,Tm));	% Plane wave at cent of T-, 1xNe

rp = obj.cent(:,Tp) - obj.vertex(:,obj.edges(3,:));	% Rho of center in T+
rm = obj.vertex(:,obj.edges(4,:)) - obj.cent(:,Tm);	% Rho of center in T-

Ei = obj.ln .* (sum((e_i*pp).*rp) + sum((e_i*pm).*rm)) /2; Ei = Ei.';

% (un x Hi).rho = Hi.(rho x un)
Hi = obj.ln .* (sum((h_i*pp).*cross(rp,obj.un(:,Tp))) + sum((h_i*pm).*cross(rm,obj.un(:,Tm))) ) /(2*eta); Hi = Hi.';



