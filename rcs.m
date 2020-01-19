% RCSb = rcs(topol,vertex,edges,ln,cent,J,k,eta,r)
% Bistatic rcs, in square meters
%
% Input:
% topol = topology matrix (vertex of each triangle), 3 x Nt
% vertex = vertex matrix, 3 x Nv
% edges = Edges matrix, 4 x Ne. For each edge (column):
% 		Row 1: Triangle T+
%		Row 2: Triangle T-
%		Row 3: Global number of opposite vertex in T+
%		Row 4: Global number of opposite vertex in T-
% ln	= Length of edges, 1 x Ne
% cent 	= Centroid of each triangle, 3 x Nt
% J	= Normal current at edges (unknonws in MoM)
% k	= Wave number 2*pi/lam
% eta	= Wave impedance of the medium
% r	= Direction of observation (from origin to observer)
%
% IE-MEI v3.1, Juan M. Rius, January 1997

function RCSb = rcs(obj,J,k,eta,r)

Tp = obj.edges(1,:); Tm = obj.edges(2,:); % T+/- triangles for edges 1:Ne

rho_p = obj.cent(:,Tp) - obj.vertex(:,obj.edges(3,:));	% rho of center in T+
rho_m = obj.vertex(:,obj.edges(4,:)) - obj.cent(:,Tm);	% rho of center in T-

phase_p = exp(1j*k*r' * obj.cent(:,Tp)); % phase of T+ contributions
phase_m = exp(1j*k*r' * obj.cent(:,Tm)); % phase of T- contributions

% Compute ( rho_p x r exp(Tp) + rho_m x r exp(Tm) )
% tmp = cross(rho_p,r) .* (ones(3,1)*phase_p) + cross(rho_m,r) .* (ones(3,1)*phase_m);
% Faster: don't call 'cross'
rho_pxr = [rho_p(2,:).*r(3) - rho_p(3,:).*r(2)
           rho_p(3,:).*r(1) - rho_p(1,:).*r(3)
           rho_p(1,:).*r(2) - rho_p(2,:).*r(1)];

rho_mxr = [rho_m(2,:).*r(3) - rho_m(3,:).*r(2)
           rho_m(3,:).*r(1) - rho_m(1,:).*r(3)
           rho_m(1,:).*r(2) - rho_m(2,:).*r(1)];

tmp = rho_pxr .* (ones(3,1)*phase_p) + rho_mxr .* (ones(3,1)*phase_m);

% Multiply by ( J*ln ) and sum all the edges (J is col, ln is row)
tmp2 = tmp * (J .* obj.ln.');

RCSb = norm(tmp2)^2 * (k*eta)^2 / (16*pi);

