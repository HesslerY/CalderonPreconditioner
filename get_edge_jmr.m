% [obj] = get_edge(obj)
%
% Input:
% topol = topology matrix (vertex of each triangle), 3 x Nt
% vertex = vertex matrix, 3 x Nv
%
% Output:
% edges = Edges matrix, 4 x Ne. For each edge (column):
% 		Row 1: Triangle T+
%		Row 2: Triangle T-
%		Row 3: Global number of opposite vertex in T+
%		Row 4: Global number of opposite vertex in T-
% ln	= Length of edges, 1 x Ne
% trian = triangles matrix. For each triangle (column):
%		Row 1: Edge 1 (opposite is vertex 1)
%		Row 2: Edge 2 (opposite is vertex 2)
%		Row 3: Edge 3 (opposite is vertex 3)
%		If >0, T+ for that edge; if <0, T- for that edge
%		If==0, it is edge in open boundary
% un	= Unit normal to each triangle,	3 x Nt
% ds 	= Area of each triangle,	1 x Nt
% cent 	= Centroid of each triangle,	3 x Nt
%
% IE-MEI v3.0, Juan M. Rius, January 1997
% Modified in 2018 for input/output obj

function [obj] = get_edge(obj)

topol = obj.topol;
vertex = obj.vertex;

v1=vertex(:,topol(1,:));
v2=vertex(:,topol(2,:));
v3=vertex(:,topol(3,:));

cent =(v1+v2+v3)/3;
c    = cross(v3-v1,v2-v1);
un   = unitary(c);
ds   = sqrt(sum(c.^2))/2;

Nt = length(ds);
trian = zeros(3,Nt);
edges = zeros(4,ceil(Nt*3/2));

eg = 1;			% Global number of current edge
for Tp = 1:Nt		% Triangle T+
for el = [3 2 1]	% Local edge, scan order [3 2 1] for compatibility
	if ~trian(el,Tp)		% Edge not found yet
		% Find vertex of this edge
		ver = topol(find([1 2 3]~=el),Tp);
		if el == 2, ver = flipud(ver);
        end

		[tmp T1] = find(topol == ver(1)); % T1 = triangles that have ver(1)
		[tmp T2] = find(topol == ver(2)); % T2 = triangles that have ver(2)

		[TT1 TT2] = meshgrid(T1,T2);
		Tcom = TT1(TT1==TT2);	% Triangles with common edge

		if length(Tcom) == 2
			trian(el,Tp) = eg;
			edges(1,eg) = Tp;		% T+
			edges(3,eg) = topol(el,Tp);

			Tm = Tcom(find(Tcom~=Tp)); 	% T-, not equal to Tp
			edges(2,eg) = Tm;

			% Find the local number of eg in T-
			v = find(topol(:,Tm)~=ver(1) & (topol(:,Tm)~=ver(2)));
			trian(v,Tm) = -eg;
			edges(4,eg) = topol(v,Tm);

			ln(eg) = norm(vertex(:,ver(1))-vertex(:,ver(2)));
			eg = eg+1;
		elseif length(Tcom) > 2, error('More that 2 triangles share an edge');
		end
	end
end
end

edges = edges(:,1:eg-1);   % Remove void edges, in case of open object

obj.edges = edges;
obj.trian = trian;
obj.ln = ln;
obj.un = un;
obj.ds = ds;
obj.cent = cent;

