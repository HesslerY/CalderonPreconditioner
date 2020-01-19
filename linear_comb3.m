% lc_mat = linear_comb4(obj1, obj2, trian_ch)
% Linear combination matrix obj RWGs as linear combination of obj2 RWGs
%
% obj1 = original (large) RWGs
% obj2 = refined (3x) RWGs

function lc_mat = linear_comb3(obj1, obj2, trian_ch)

Ne1 = length(obj1.ln);
%lc_mat = spalloc(Ne1, 3*Ne1, Ne1*7);
lc_mat = zeros(Ne1, 3*Ne1);

for E = 1:Ne1    % for all original RWG, find refined RWG inside
   
    Tp = obj1.edges(1,E);
    Tm = obj1.edges(2,E);
    vp = obj1.edges(3,E);
    vm = obj1.edges(4,E);
    
    vpg = obj1.topol(:,Tp);    % Golbal vertices of Tp, ramdom order
%    vpli = circshift([1 2 3], 2-find(vpg==vp)); % Sorted local indices of vertices, vp at the middle
    vpli = circshift([1 2 3], 1-find(vpg==vp)); % Sorted local indices of vertices, starting from vp
    Tpch = trian_ch(Tp, vpli); % Sorted children of Tp, the one sharing edge E is last
    
    vmg = obj1.topol(:,Tm);    % Golbal vertices of Tp, ramdom order
%    vmli = circshift([1 2 3], 2-find(vmg==vm)); % Sorted local indices of vertices, vn at the middle
    vmli = circshift([1 2 3], 1-find(vmg==vm)); % Sorted local indices of vertices, starting from vm
    Tmch = trian_ch(Tm, vmli); % Sorted children of Tm, the one sharing edge E is last 

    l1 = nonzeros( intersect(obj2.trian(:,Tpch(2)),-obj2.trian(:,Tpch(1))) );
    l2 = nonzeros( intersect(obj2.trian(:,Tpch(2)),-obj2.trian(:,Tpch(3))) );
    l3 = nonzeros( intersect(obj2.trian(:,Tpch(3)),-obj2.trian(:,Tpch(1))) );
    l4 = nonzeros( intersect(obj2.trian(:,Tpch(1)),-obj2.trian(:,Tmch(1))) );
    l5 = nonzeros( intersect(obj2.trian(:,Tmch(1)),-obj2.trian(:,Tmch(2))) );
    l6 = nonzeros( intersect(obj2.trian(:,Tmch(2)),-obj2.trian(:,Tmch(3))) );
    l7 = nonzeros( intersect(obj2.trian(:,Tmch(1)),-obj2.trian(:,Tmch(3))) );
    
    L = obj1.ln(E);
    lc_mat(E, abs(l1)) = sign(l1)*L/(3*obj2.ln(abs(l1)));
    lc_mat(E, abs(l2)) = 0;
    lc_mat(E, abs(l3)) = sign(l3)*L/(3*obj2.ln(abs(l3)));
    lc_mat(E, abs(l4)) = sign(l4)*1;
    lc_mat(E, abs(l5)) = sign(l5)*L/(3*obj2.ln(abs(l5)));
    lc_mat(E, abs(l6)) = 0;
    lc_mat(E, abs(l7)) = sign(l7)*L/(3*obj2.ln(abs(l7)));
end

end


