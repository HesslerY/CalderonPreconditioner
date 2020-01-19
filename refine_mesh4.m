function [obj2, trian_ch, trian_pr] = refine_mesh4(obj)

Nt = size(obj.trian,2);
Nv = size(obj.vertex,2);
obj2 = struct('topol', [], 'vertex', obj.vertex);

trian_ch = zeros(Nt, 4); % Children of each triangle, the first 3 located at corresponding vertex, the 4th is the center one
trian_pr = zeros(4*Nt, 2); % Parent of each triangle, and local position within parent

nt2 = 0; % Index of triangle in refined mesh

for t = 1:Nt			% For each triangle
	v1 = obj.topol(1,t); v2 = obj.topol(2,t); v3 = obj.topol(3,t);
	r12 = (obj.vertex(:,v1) + obj.vertex(:,v2))/2;
	r23 = (obj.vertex(:,v2) + obj.vertex(:,v3))/2;
	r31 = (obj.vertex(:,v3) + obj.vertex(:,v1))/2;

	% Check if new obj.vertex already exists. If not, append to obj2.vertex
	e12 = find(all([ r12(1)==obj2.vertex(1,:); r12(2)==obj2.vertex(2,:); r12(3)==obj2.vertex(3,:) ]));
	e23 = find(all([ r23(1)==obj2.vertex(1,:); r23(2)==obj2.vertex(2,:); r23(3)==obj2.vertex(3,:) ]));
	e31 = find(all([ r31(1)==obj2.vertex(1,:); r31(2)==obj2.vertex(2,:); r31(3)==obj2.vertex(3,:) ]));
	if e12, v12 = e12; 
	else
        v12 = Nv+1; Nv = Nv+1; obj2.vertex = [obj2.vertex r12];
	end
    
	if e23, v23 = e23; 
	else
        v23 = Nv+1; Nv = Nv+1; obj2.vertex = [obj2.vertex r23];
	end
    
	if e31, v31 = e31; 
	else
        v31 = Nv+1; Nv = Nv+1; obj2.vertex = [obj2.vertex r31];
	end

	obj2.topol = [ obj2.topol [v1; v12; v31] [v2; v23; v12] [v3; v31; v23] [v12; v23; v31]  ];
    
    trian_ch(t, :) = nt2 + [1 2 3 4];
    trian_pr(nt2 + [1 2 3 4], :) = [t 1; t 2; t 3; t 4];
    
    nt2 = nt2 + 4;
end