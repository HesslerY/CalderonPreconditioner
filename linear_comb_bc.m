% lc_mat = linear_comb_bc(obj1, obj2, trian_ch)
% Linear combination matrix obj Buffa-Christiansen as linear combination of obj2 RWGs
%
% obj1 = original (large) RWGs
% obj2 = refined (6x) barycentric RWGs

function lc_mat_bc = linear_comb_bc(obj1, obj2, trian_ch)

Ne1 = length(obj1.ln);
Nv2 = size(obj2.vertex, 2);
%lc_mat = spalloc(Ne1, 6*Ne1, Ne1*10);
lc_mat_bc = zeros(Ne1, 6*Ne1);

for E = 1:Ne1    % for all original RWG, find refined RWG inside
   
    Tp = obj1.edges(1,E);
    Tm = obj1.edges(2,E);
    vp = obj1.edges(3,E);
    vm = obj1.edges(4,E);
    
    %% Barycentric edges
    vpg = obj1.topol(:,Tp);    % Golbal vertices of Tp, ramdom order
    vpli = circshift([1 2 3 4 5 6], -2+2*find(vpg==vp)); % Sorted local indices of vertices, starting from vp
    Tpch = trian_ch(Tp, vpli); % Sorted children of Tp, the one right to vp is 1st
    
    vmg = obj1.topol(:,Tm);    % Golbal vertices of Tp, ramdom order
    vmli = circshift([1 2 3 4 5 6], -2+2*find(vmg==vm)); % Sorted local indices of vertices, starting from vm
    Tmch = trian_ch(Tm, vmli); % Sorted children of Tm, the one right to vm is 1st

%     l1b  = nonzeros( intersect(obj2.trian(:,Tpch(6)),-obj2.trian(:,Tpch(1))) );
%     l2b  = nonzeros( intersect(obj2.trian(:,Tpch(6)),-obj2.trian(:,Tpch(5))) );
%     l3b  = nonzeros( intersect(obj2.trian(:,Tpch(1)),-obj2.trian(:,Tpch(2))) );
%     
%     l4b  = nonzeros( intersect(obj2.trian(:,Tpch(5)),-obj2.trian(:,Tpch(4))) );
    l5b  = nonzeros( intersect(-obj2.trian(:,Tpch(3)), obj2.trian(:,Tpch(4))) );  % edge sign in BC definition opposite to barycentric mesh
%     l6b  = nonzeros( intersect(obj2.trian(:,Tpch(2)),-obj2.trian(:,Tpch(3))) );
%     
     l7b  = nonzeros( intersect(obj2.trian(:,Tpch(4)),-obj2.trian(:,Tmch(3))) );
    l8b  = nonzeros( intersect(obj2.trian(:,Tpch(3)),-obj2.trian(:,Tmch(4))) );
    
%    l9b  = nonzeros( intersect(obj2.trian(:,Tmch(3)),-obj2.trian(:,Tmch(2))) );
    l10b = nonzeros( intersect(-obj2.trian(:,Tmch(3)), obj2.trian(:,Tmch(4))) );  % edge sign in BC definition opposite to barycentric mesh
%     l11b = nonzeros( intersect(obj2.trian(:,Tmch(4)),-obj2.trian(:,Tmch(5))) );
%     
%     l12b = nonzeros( intersect(obj2.trian(:,Tmch(2)),-obj2.trian(:,Tmch(1))) );
%     l13b = nonzeros( intersect(obj2.trian(:,Tmch(6)),-obj2.trian(:,Tmch(1))) );
%     l14b = nonzeros( intersect(obj2.trian(:,Tmch(5)),-obj2.trian(:,Tmch(6))) );
    
    %% Right vertex
    al_right = l5b; % Array of edges lr of right vertex
    ln = l8b;
    
    Tn = Tpch(3);
    aT_right = [ ];
    
    for n = 1:Nv2
        if n>1, al_right = [al_right ln]; %#ok<AGROW>
        end
        aT_right = [aT_right Tn]; %#ok<AGROW>

        % Find new ln, which is the edge before old ln in Tn of obj2
        ln3 = obj2.trian(:,Tn); % All 3 edges of Tn
        pos = find(abs(ln3) == abs(ln));  % Position of ln in Tn
        tmp = circshift([1 2 3],3-pos); % Bring old ln to 2nd position, new ln is then 1st
        ln = ln3(tmp(1));   % New ln
        
        % Find new Tn, which is the other trian sharing ln
        pos = find(obj2.edges(1:2,abs(ln)) == Tn); % Position of old Tn in ln
        Tn = obj2.edges(3-pos,abs(ln));  % New Tn
        
        if abs(ln) == abs(l8b), break;
        end
    end
    
    %% Left vertex
    al_left = l10b; % Array of edges ll of left vertex
    ln = l7b;
    
    Tn = Tmch(3);
    aT_left = [ ];
    
    for n = 1:Nv2
        if n>1, al_left = [al_left ln]; %#ok<AGROW>
        end
        aT_left = [aT_left Tn]; %#ok<AGROW>

        % Find new ln, which is the edge before old ln in Tn of obj2
        ln3 = obj2.trian(:,Tn); % All 3 edges of Tn
        pos = find(abs(ln3) == abs(ln));  % Position of ln in Tn
        tmp = circshift([1 2 3],3-pos); % Bring old ln to 2nd position, new ln is then 1st
        ln = ln3(tmp(1));   % New ln
        
        % Find new Tn, which is the other trian sharing ln
        pos = find(obj2.edges(1:2,abs(ln)) == Tn); % Position of old Tn in ln
        Tn = obj2.edges(3-pos,abs(ln));  % New Tn
        
        if abs(ln) == abs(l7b), break;
        end
    end
    
    %% Linear combination coefficients
%    L = obj1.ln(E); % Pura lógica
    L = 1; % Andriulli's paper
    Ncr = size(aT_right,2)/2; % Number of right triangles
    n = 0:2*Ncr -1;         % Indices of right edges
    cr = L * (Ncr-n) ./ ( 2*Ncr* obj2.ln(abs(al_right) ));  % Coefficients for right edges
    lc_mat_bc(E, abs(al_right)) = sign(al_right).*cr;
    
    Ncl = size(aT_left,2)/2; % Number of left triangles
    n = 0:2*Ncl -1;         % Indices of left edges
    cl = -L * (Ncl-n) ./ ( 2*Ncl* obj2.ln(abs(al_left) ));  % Coefficients for left edges
    lc_mat_bc(E, abs(al_left)) = sign(al_left).*cl;
    
% Como en el paper de Buffa-Christiansen, sin dividir por li
%     Ncr = size(aT_right,2)/2; % Number of right triangles
%     n = 0:2*Ncr -1;         % Indices of right edges
%     cr = (Ncr-n) ./ ( 2*Ncr);  % Coefficients for right edges
%     lc_mat_bc(E, abs(al_right)) = sign(al_right).*cr;
%     
%     Ncl = size(aT_left,2)/2; % Number of left triangles
%     n = 0:2*Ncl -1;         % Indices of left edges
%     cl = -(Ncl-n) ./ ( 2*Ncl);  % Coefficients for left edges
%     lc_mat_bc(E, abs(al_left)) = sign(al_left).*cl;    
    
end

end


