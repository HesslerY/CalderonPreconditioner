% lc_mat = linear_comb6(obj1, obj2, trian_ch)
% Linear combination matrix obj RWGs as linear combination of obj2 RWGs
%
% obj1 = original (large) RWGs
% obj2 = refined (6x) barycentric RWGs

function lc_mat = linear_comb6(obj1, obj2, trian_ch)

Ne1 = length(obj1.ln);
%lc_mat = spalloc(Ne1, 6*Ne1, Ne1*10);
lc_mat = zeros(Ne1, 6*Ne1);

for E = 1:Ne1    % for all original RWG, find refined RWG inside
   
    Tp = obj1.edges(1,E);
    Tm = obj1.edges(2,E);
    vp = obj1.edges(3,E);
    vm = obj1.edges(4,E);
    
    vpg = obj1.topol(:,Tp);    % Golbal vertices of Tp, ramdom order
    vpli = circshift([1 2 3 4 5 6], -2+2*find(vpg==vp)); % Sorted local indices of subtriangles, starting from vp
    Tpch = trian_ch(Tp, vpli); % Sorted children of Tp, the one right to vp is 1st
    
    vmg = obj1.topol(:,Tm);    % Golbal vertices of Tm, ramdom order
    vmli = circshift([1 2 3 4 5 6], -2+2*find(vmg==vm)); % Sorted local indices of subtriangles, starting from vp
    Tmch = trian_ch(Tm, vmli); % Sorted children of Tm, the one right to vm is 1st

    l1  = nonzeros( intersect(obj2.trian(:,Tpch(6)),-obj2.trian(:,Tpch(1))) );
    l2  = nonzeros( intersect(obj2.trian(:,Tpch(6)),-obj2.trian(:,Tpch(5))) );
    l3  = nonzeros( intersect(obj2.trian(:,Tpch(1)),-obj2.trian(:,Tpch(2))) );
    
    l4  = nonzeros( intersect(obj2.trian(:,Tpch(5)),-obj2.trian(:,Tpch(4))) );
    l5  = nonzeros( intersect(obj2.trian(:,Tpch(3)),-obj2.trian(:,Tpch(4))) );
    l6  = nonzeros( intersect(obj2.trian(:,Tpch(2)),-obj2.trian(:,Tpch(3))) );
    
    l7  = nonzeros( intersect(obj2.trian(:,Tpch(4)),-obj2.trian(:,Tmch(3))) );
    l8  = nonzeros( intersect(obj2.trian(:,Tpch(3)),-obj2.trian(:,Tmch(4))) );
    
    l9  = nonzeros( intersect(obj2.trian(:,Tmch(3)),-obj2.trian(:,Tmch(2))) );
    l10 = nonzeros( intersect(obj2.trian(:,Tmch(3)),-obj2.trian(:,Tmch(4))) );
    l11 = nonzeros( intersect(obj2.trian(:,Tmch(4)),-obj2.trian(:,Tmch(5))) );
    
    l12 = nonzeros( intersect(obj2.trian(:,Tmch(2)),-obj2.trian(:,Tmch(1))) );
    l13 = nonzeros( intersect(obj2.trian(:,Tmch(6)),-obj2.trian(:,Tmch(1))) );
    l14 = nonzeros( intersect(obj2.trian(:,Tmch(5)),-obj2.trian(:,Tmch(6))) );
    
    L = obj1.ln(E);
    lc_mat(E, abs( l1)) =  0;
    lc_mat(E, abs( l2)) =  sign( l2)*L/(6*obj2.ln(abs( l2)));
    lc_mat(E, abs( l3)) =  sign( l3)*L/(6*obj2.ln(abs( l3)));
    lc_mat(E, abs( l4)) =  sign( l4)*L/(3*obj2.ln(abs( l4)));
    lc_mat(E, abs( l5)) =  0;
    lc_mat(E, abs( l6)) =  sign( l6)*L/(3*obj2.ln(abs( l6)));
    lc_mat(E, abs( l7)) =  sign( l7)*1;
    lc_mat(E, abs( l8)) =  sign( l8)*1;
    lc_mat(E, abs( l9)) =  sign( l9)*L/(3*obj2.ln(abs( l9)));
    lc_mat(E, abs(l10)) =  0;
    lc_mat(E, abs(l11)) =  sign(l11)*L/(3*obj2.ln(abs(l11)));
    lc_mat(E, abs(l12)) =  sign(l12)*L/(6*obj2.ln(abs(l12)));
    lc_mat(E, abs(l13)) =  0;
    lc_mat(E, abs(l14)) =  sign(l14)*L/(6*obj2.ln(abs(l14)));
end

end


