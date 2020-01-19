% lc_mat = linear_comb4(obj1, obj2, trian_ch)
% Linear combination matrix obj RWGs as linear combination of obj2 RWGs
%
% obj1 = original (large) RWGs
% obj2 = refined (4x) RWGs

function lc_mat = linear_comb4(obj1, obj2, trian_ch)

Ne1 = length(obj1.ln);
%lc_mat = spalloc(Ne1, 4*Ne1, Ne1*8);
lc_mat = zeros(Ne1, 4*Ne1);

for E = 1:Ne1    % for all original RWG, find refined RWG inside
   
    Tp = obj1.edges(1,E);
    Tm = obj1.edges(2,E);
    vp = obj1.edges(3,E);
    vm = obj1.edges(4,E);
    
    vpg = obj1.topol(:,Tp);    % Golbal vertices of Tp, ramdom order
    vpli = circshift([1 2 3], 1-find(vpg==vp)); % Local indices of vertices, starting from vp
    Tpch = trian_ch(Tp, [vpli 4]); % Sorted children of Tp
    
    vmg = obj1.topol(:,Tm);    % Golbal vertices of Tp, ramdom order
    vmli = circshift([1 2 3], 1-find(vmg==vm)); % Local indices of vertices, starting from vm
    Tmch = trian_ch(Tm, [vmli 4]); % Sorted children of Tp    
    

    l1 = nonzeros( intersect(obj2.trian(:,Tpch(1)),-obj2.trian(:,Tpch(4))) );
    l2 = nonzeros( intersect(obj2.trian(:,Tpch(4)),-obj2.trian(:,Tpch(2))) );
    l3 = nonzeros( intersect(obj2.trian(:,Tpch(4)),-obj2.trian(:,Tpch(3))) );
    l4 = nonzeros( intersect(obj2.trian(:,Tpch(2)),-obj2.trian(:,Tmch(3))) );
    l5 = nonzeros( intersect(obj2.trian(:,Tpch(3)),-obj2.trian(:,Tmch(2))) );
    l6 = nonzeros( intersect(obj2.trian(:,Tmch(3)),-obj2.trian(:,Tmch(4))) );
    l7 = nonzeros( intersect(obj2.trian(:,Tmch(2)),-obj2.trian(:,Tmch(4))) );
    l8 = nonzeros( intersect(obj2.trian(:,Tmch(4)),-obj2.trian(:,Tmch(1))) );
    
    lc_mat(E, abs(l1)) = sign(l1)*0.5;
    lc_mat(E, abs(l2)) = sign(l2)*obj1.ln(E)/(4*obj2.ln(abs(l2)));
    lc_mat(E, abs(l3)) = sign(l3)*obj1.ln(E)/(4*obj2.ln(abs(l3)));
    lc_mat(E, abs(l4)) = sign(l4)*1;
    lc_mat(E, abs(l5)) = sign(l5)*1;
    lc_mat(E, abs(l6)) = sign(l6)*obj1.ln(E)/(4*obj2.ln(abs(l6)));
    lc_mat(E, abs(l7)) = sign(l7)*obj1.ln(E)/(4*obj2.ln(abs(l7)));
    lc_mat(E, abs(l8)) = sign(l8)*0.5;
end

end


