% D = D_mat(obj)
% D = <w_m,u_i> scalar producs of testing with basis functions
% Integration on testing functions: 1 point at cent
%
% Input:
% topol = topology matrix (vertex of each triangle), 3 x Nt
% vertex = vertex matrix, 3 x Nv
% trian = triangles matrix. For each triangle (column):
%		Row 1: Edge 1 (opposite is vertex 1)
%		Row 2: Edge 2 (opposite is vertex 2)
%		Row 3: Edge 3 (opposite is vertex 3)
%		If >0, T+ for that edge; if <0, T- for that edge
% cent 	= Centroid of triangles, 3 x Nt
% ds	= Area of triangles,	 1 x Nt
% edges = Edges matrix, 4 x Ne. For each edge (column):
% 		Row 1: Triangle T+
%		Row 2: Triangle T-
%		Row 3: Global number of opposite vertex in T+
%		Row 4: Global number of opposite vertex in T-
% ln	 = Length of edges, 1 x Ne
%
% IE-MEI v3.0, Juan M. Rius, January 1997

function D = D_mat(obj)

Ne = size(obj.edges,2);
Nt = size(obj.trian,2);
D = sparse(Ne,Ne);

for T=1:Nt
	ge = obj.trian(:,T);	% Global edges of this triangle
	le = (ge~=0);		% Local numbers of interior edges (~=0)
	ge = ge(le);		% Remove edges=0 (boundary edges)
	Nge = length(ge);	% Maybe < 3 if edges=0 have been removed
	si = sign(ge);		% Sign of T for this edges
	ge = abs(ge);

	% Aproximation: Rho at centroid -> Produces singular D matrix
	% Array of rho vectors for each edge, with +/- sign
	% rho = ones(3,1)*si.' .* (cent(:,T)*ones(1,Nge) - vertex(:,topol(le,T)));
	% D(ge,ge) = D(ge,ge) + ds(T) * (ln(ge).'*ln(ge)) .* (rho.'*rho) / (4*ds(T)^2);

	v1 = obj.vertex(:,obj.topol(1,T));
	v2 = obj.vertex(:,obj.topol(2,T));
	v3 = obj.vertex(:,obj.topol(3,T));

	% Gauss integration:
	% For quadratic functions it is exact with 3 points (error 5e-18)
	tmp = 0;

	r = v1*2/3 + v2/6 + v3/6;
	rho = (ones(3,1)*si.') .* (r*ones(1,Nge) - obj.vertex(:,obj.topol(le,T)));
	tmp = tmp + (rho.'*rho);

	r = v2*2/3 + v3/6 + v1/6;
	rho = (ones(3,1)*si.') .* (r*ones(1,Nge) - obj.vertex(:,obj.topol(le,T)));
	tmp = tmp + (rho.'*rho);

	r = v3*2/3 + v1/6 + v2/6;
	rho = (ones(3,1)*si.') .* (r*ones(1,Nge) - obj.vertex(:,obj.topol(le,T)));
	tmp = tmp + (rho.'*rho);

	tmp = tmp/6;	% weigth = 1/3 x area = 1/2

% 	% Exact computation: Slightly faster than integration
% 	rho_3m = [v1-v3 v2-v3 v3-v3]; rho_3n = rho_3m;
% 	r31 = rho_3m(:,1); r32 = rho_3m(:,2);
% 	l31 = norm(r31);
% 	l32 = norm(r32);
%
% 	tmp = (si*si.').*( l31^2 + l32^2 + r31.'*r32 + 6*rho_3m(:,le).'*rho_3n(:,le) - 2*(rho_3m(:,le).'*(r31+r32))*ones(1,length(ge)) - 2*ones(length(ge),1)*((r31+r32).'*rho_3n(:,le))  )/12;
	D(ge,ge) = D(ge,ge) + (obj.ln(ge).'*obj.ln(ge)) .* tmp / (2*obj.ds(T));
end

