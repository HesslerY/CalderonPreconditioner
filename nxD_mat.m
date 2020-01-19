% nxD = nxD_mat(obj)
% D = <w_m,f_n> scalar producs of testing with basis fobj.unctions
% Integration on testing fobj.unctions: 1 point at cent
%
% w_m = nxRWG, f_n = RWG
%
% Input:
% obj.topol = obj.topology matrix (obj.vertex of each obj.triangle), 3 x Nt
% obj.vertex = vertex matrix, 3 x Nv
% obj.trian = triangles matrix. For each obj.triangle (column):
%		Row 1: Edge 1 (opposite is obj.vertex 1)
%		Row 2: Edge 2 (opposite is obj.vertex 2)
%		Row 3: Edge 3 (opposite is obj.vertex 3)
%		If >0, T+ for that edge; if <0, T- for that edge
% obj.un 	= Unit normal to obj.triangles, 3 x Nt
% ds	= Area of obj.triangles,	 1 x Nt
% obj.edges = Edges matrix, 4 x Ne. For each edge (column):
% 		Row 1: Triangle T+
%		Row 2: Triangle T-
%		Row 3: Global number of opposite obj.vertex in T+
%		Row 4: Global number of opposite obj.vertex in T-
% obj.ln	 = Length of obj.edges, 1 x Ne
%
% IE-MEI v3.0, Juan M. Rius, January 1997

function nxD = nxD_mat(obj)

Ne = size(obj.edges,2);
Nt = size(obj.trian,2);
nxD = sparse(Ne,Ne);

for T=1:Nt
	ge = obj.trian(:,T);	% Global obj.edges of this obj.triangle
	le = (ge~=0);		% Local numbers of interior obj.edges (~=0)
	ge = ge(le);		% Remove obj.edges=0 (boobj.undary obj.edges)
	Nge = length(ge);	% Maybe < 3 if obj.edges=0 have been removed
	si = sign(ge);		% Sign of T for this obj.edges
	ge = abs(ge);

    unT = obj.un(:,T)*ones(1,Nge);  % Unit normal to this obj.triangle, for all the obj.edges

	% Aproximation: Rho at centroid -> Produces singular D matrix
	% Array of rho vectors for each edge, with +/- sign
	% rho = ones(3,1)*si.' .* (cent(:,T)*ones(1,Nge) - obj.vertex(:,obj.topol(le,T)));
	% D(ge,ge) = D(ge,ge) + ds(T) * (obj.ln(ge).'*obj.ln(ge)) .* (rho.'*rho) / (4*ds(T)^2);

	v1 = obj.vertex(:,obj.topol(1,T));
    v2 = obj.vertex(:,obj.topol(2,T));
    v3 = obj.vertex(:,obj.topol(3,T));

	% Gauss integration:
	% For quadratic fobj.unctions it is exact with 3 points (error 5e-18)
	tmp = 0;

	r = v1*2/3 + v2/6 + v3/6;
	rho = (ones(3,1)*si.') .* (r*ones(1,Nge) - obj.vertex(:,obj.topol(le,T)));
    nxrho = cross(unT,rho,1);
	tmp = tmp + (nxrho.'*rho);

	r = v2*2/3 + v3/6 + v1/6;
	rho = (ones(3,1)*si.') .* (r*ones(1,Nge) - obj.vertex(:,obj.topol(le,T)));
    nxrho = cross(unT,rho,1);
	tmp = tmp + (nxrho.'*rho);

	r = v3*2/3 + v1/6 + v2/6;
	rho = (ones(3,1)*si.') .* (r*ones(1,Nge) - obj.vertex(:,obj.topol(le,T)));
    nxrho = cross(unT,rho,1);
	tmp = tmp + (nxrho.'*rho);

	tmp = tmp/6;	% weigth = 1/3 x area = 1/2

	nxD(ge,ge) = nxD(ge,ge) + (obj.ln(ge).'*obj.ln(ge)) .* tmp / (2*obj.ds(T));
end

