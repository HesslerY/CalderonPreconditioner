% function [obj2, trian_ch, trian_pr] = refine_mesh3(obj)
% Refone mesh obj
%
% New triangles: 
% Tn, where n is the position within the parent triangle
% T1 replaces local vertex 3 = 
% T2 replaces local vertex 1
% T3 replaces local vertex 2

function [obj2, trian_ch, trian_pr] = refine_mesh3(obj)

Nt = size(obj.trian,2);
Nv = size(obj.vertex,2);
obj2 = struct('topol', [], 'vertex', obj.vertex);

trian_ch = zeros(Nt, 3); % Children of each triangle
trian_pr = zeros(3*Nt, 2); % Parent of each triangle, and local position within parent

nt2 = 0; % Index of triangle in refined mesh

for t = 1:Nt			% For each triangle
	v1 = obj.topol(1,t); v2 = obj.topol(2,t); v3 = obj.topol(3,t);
    
    % Append barycenter as new vertex
    rb = ( obj.vertex(:,v1) + obj.vertex(:,v2) + obj.vertex(:,v3) )/3;
    obj2.vertex = [obj2.vertex rb];
    
    % Create new triangles and append, sorting opposite to v1, v2, v3
    vb = Nv+1; Nv = Nv+1; % Index of new vertex (barycenter)
	obj2.topol = [ obj2.topol [v2; v3; vb] [v3; v1; vb] [v1; v2; vb] ];
    
    trian_ch(t, :) = nt2 + [1 2 3];
    trian_pr(nt2 + [1 2 3], :) = [t 1; t 2; t 3];
    
    nt2 = nt2 + 3;
end