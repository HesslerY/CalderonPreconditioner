% function [obj2, trian_ch, trian_pr] = refine_mesh3(obj)
% Refine mesh obj to barycentric 6xN mesh obj2
%
% New triangles: 
% Tn, where n is the position within the parent triangle
% T1 replaces local vertex 3 = 
% T2 replaces local vertex 1
% T3 replaces local vertex 2

function [obj2, trian_ch, trian_pr] = refine_mesh6(obj)

Nt = size(obj.trian,2);
Nv = size(obj.vertex,2);
obj2 = struct('topol', [], 'vertex', obj.vertex);

trian_ch = zeros(Nt, 6); % Children of each triangle
trian_pr = zeros(6*Nt, 2); % Parent of each triangle, and local position within parent

nt2 = 0; % Index of triangle in refined mesh

for t = 1:Nt			% For each triangle
	vA = obj.topol(1,t); vF = obj.topol(3,t); vD = obj.topol(2,t);  % obj triangle vertices are clockwise
    
    rA = obj.vertex(:,vA);
    rF = obj.vertex(:,vF);
    rD = obj.vertex(:,vD);
    rB = (rD + rA     )/2;
    rC = (rA + rF     )/2;
	rE = (rF + rD     )/2;
    rL = (rA + rF + rD)/3;

	% Check if new vertex already exists. If not, append to obj2.vertex
	eB = find(all([ rB(1)==obj2.vertex(1,:); rB(2)==obj2.vertex(2,:); rB(3)==obj2.vertex(3,:) ]));
	eC = find(all([ rC(1)==obj2.vertex(1,:); rC(2)==obj2.vertex(2,:); rC(3)==obj2.vertex(3,:) ]));
	eE = find(all([ rE(1)==obj2.vertex(1,:); rE(2)==obj2.vertex(2,:); rE(3)==obj2.vertex(3,:) ]));
    eL = find(all([ rL(1)==obj2.vertex(1,:); rL(2)==obj2.vertex(2,:); rL(3)==obj2.vertex(3,:) ]));
    
	if eB, vB = eB; 
	else
        vB = Nv+1; Nv = Nv+1; obj2.vertex = [obj2.vertex rB];
	end
    
	if eC, vC = eC; 
	else
        vC = Nv+1; Nv = Nv+1; obj2.vertex = [obj2.vertex rC];
	end
    
	if eE, vE = eE; 
	else
        vE = Nv+1; Nv = Nv+1; obj2.vertex = [obj2.vertex rE];
    end    

	if eL, vL = eL; 
	else
        vL = Nv+1; Nv = Nv+1; obj2.vertex = [obj2.vertex rL];
	end        
    
    
    % Create new triangles and append, sorting starting T1 just right to v1=vA=obj.topol(1,t)
	obj2.topol = [ obj2.topol [vA; vL; vC] [vC; vL; vF] [vF; vL; vE] [vE; vL; vD] [vD; vL; vB] [vB; vL; vA] ];
    
    trian_ch(t, :) = nt2 + [1 2 3 4 5 6];
    trian_pr(nt2 + [1 2 3 4 5 6], :) = [t 1; t 2; t 3; t 4; t 5; t 6];
    
    nt2 = nt2 + 6;
end