% D2 = D4_mat(obj, obj2, trian_par)
% D2 = <w_m,f_n> scalar producs of testing RWG (original mesh) with basis RWG (refined mesh)
% Integration on testing functions: 3 Gauss points (on the smaller triangle -> refined mesh)
%
% w_m = RWG, f_n = RWG_refined
%
% Input: 
%   obj     = original testing mesh
%   obj2    = refined basis mesh
%   trian_par = parent triangles of each refined mesh triangle, output from function refine_mesh()
%
% With:
% obj.topol = obj.topology matrix (obj.vertex of each obj.triangle), 3 x Nt
% obj.vertex = vertex matrix, 3 x Nv
% obj.trian = triangles matrix. For each obj.triangle (column):
%		Row 1: Edge 1 (opposite is obj.vertex 1)
%		Row 2: Edge 2 (opposite is obj.vertex 2)
%		Row 3: Edge 3 (opposite is obj.vertex 3)
%		If >0, T+ for that edge; if <0, T- for that edge
% obj.un 	= Unit normal to obj.triangles, 3 x Nt
% obj.ds	= Area of obj.triangles,	 1 x Nt
% obj.edges = Edges matrix, 4 x Ne. For each edge (column):
% 		Row 1: Triangle T+
%		Row 2: Triangle T-
%		Row 3: Global number of opposite obj.vertex in T+
%		Row 4: Global number of opposite obj.vertex in T-
% obj.ln	 = Length of obj.edges, 1 x Ne
%
% IE-MEI v3.0, Juan M. Rius, January 1997

function D2 = D4_mat(obj, obj2, trian_pr)

Ne  = size(obj.edges,  2);
Ne2 = size(obj2.edges, 2);
Nt2 = size(obj2.trian, 2);

%D2 = sparse(Ne,Ne2);
D2 = zeros(Ne,Ne2);

for T2=1:Nt2  % Basis triangle, that belongs up to 3 refined basis functions
    
    % RWGs sharing this triangle, in the refined mesh obj2
	ge2 = obj2.trian(:,T2);	% Global obj.edges of this obj.triangle
	le2 = (ge2~=0);         % Local numbers of interior obj.edges (~=0)
	ge2 = ge2(le2);         % Remove obj.edges=0 (boundary obj.edges)
	Nge2 = length(ge2);     % Maybe < 3 if obj.edges=0 have been removed
	si2 = sign(ge2);		% Sign of T for this obj.edges
	ge2 = abs(ge2);
    
    % Parent triangle in original mesh, and RWGs sharing this triangle
    T = trian_pr(T2, 1); % Parent triangle, of testing original RWG
	ge = obj.trian(:,T);	% Global obj.edges of this obj.triangle
	le = (ge~=0);		% Local numbers of interior obj.edges (~=0)
	ge = ge(le);		% Remove obj.edges=0 (boobj.undary obj.edges)
	Nge = length(ge);	% Maybe < 3 if obj.edges=0 have been removed
	si = sign(ge);		% Sign of T for this obj.edges
	ge = abs(ge);    

    unT = obj.un(:,T)*ones(1,Nge);  % Unit normal to this obj.triangle, for all the obj.edges

	% 1-point integration: Rho at centroid -> Produces singular D matrix
	% Array of rho vectors for each edge, with +/- sign
	% rho = ones(3,1)*si.' .* (cent(:,T)*ones(1,Nge) - obj.vertex(:,obj.topol(le,T)));
	% D(ge,ge) = D(ge,ge) + ds(T) * (obj.ln(ge).'*obj.ln(ge)) .* (rho.'*rho) / (4*ds(T)^2);

	% Vertices of integration triangle: the small triangle of refined mesh (T2)
    v1 = obj2.vertex(:,obj2.topol(1,T2));
    v2 = obj2.vertex(:,obj2.topol(2,T2));
    v3 = obj2.vertex(:,obj2.topol(3,T2));

	% Gauss integration:
	% For quadratic functions it is exact with 3 points (error 5e-18)
	tmp = 0;

	r = v1*2/3 + v2/6 + v3/6;
	rho  = (ones(3,1)*si.' ) .* (r*ones(1,Nge ) - obj.vertex(:,  obj.topol( le, T )));
    rho2 = (ones(3,1)*si2.') .* (r*ones(1,Nge2) - obj2.vertex(:, obj2.topol(le2,T2)));
%    nxrho = cross(unT,rho,1);
	tmp = tmp + (rho.'*rho2);

	r = v2*2/3 + v3/6 + v1/6;
	rho  = (ones(3,1)*si.' ) .* (r*ones(1,Nge ) - obj.vertex(:,  obj.topol( le, T )));
    rho2 = (ones(3,1)*si2.') .* (r*ones(1,Nge2) - obj2.vertex(:, obj2.topol(le2,T2)));
%    nxrho = cross(unT,rho,1);
	tmp = tmp + (rho.'*rho2);

	r = v3*2/3 + v1/6 + v2/6;
	rho  = (ones(3,1)*si.' ) .* (r*ones(1,Nge ) - obj.vertex(:,  obj.topol( le, T )));
    rho2 = (ones(3,1)*si2.') .* (r*ones(1,Nge2) - obj2.vertex(:, obj2.topol(le2,T2)));
%    nxrho = cross(unT,rho,1);    
	tmp = tmp + (rho.'*rho2);

	tmp = tmp/6;	% weigth = 1/3 x area = 1/2

	D2(ge,ge2) = D2(ge,ge2) + (obj.ln(ge).'*obj2.ln(ge2)) .* tmp / (2*obj.ds(T)); %#ok<SPRIX>
end

