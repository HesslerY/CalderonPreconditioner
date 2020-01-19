% geom = 'sq_plate';
% param = struct('Nx',1, 'Lx',1, 'Nz',2, 'Lz',2, 'x',0.5, 'y',0, 'z',1, 'cor',0);
% 
% cd objects
% obj = feval(geom,param);
% obj.vertex(2,:) = obj.vertex(3,:);
% obj.vertex(3,:) = 0;
% cd ..

%obj = get_edge(obj);
%user_plot_geom3d(obj1); view(2);

%[obj2, trian_ch, trian_pr] = refine_mesh3(obj);
%[obj2, trian_ch, trian_pr] = refine_mesh4(obj);
[obj2, trian_ch, trian_pr] = refine_mesh6(obj);
obj2 = get_edge(obj2);
obj2.N = length(obj2.ln);
N2 = obj2.N;
user_plot_geom3d(obj2); 
plot_obj_numbers(obj2);

%lc_mat = linear_comb3(obj, obj2, trian_ch);
%lc_mat = linear_comb4(obj, obj2, trian_ch);
lc_mat = linear_comb6(obj, obj2, trian_ch);
R = lc_mat.';
piR = pinv(R);

%% Test Gramm matrix
D2 = D2_mat(obj, obj2, trian_pr);
D1 = full(D_mat(obj));
D1ref = D2*R;
errD = comp_mat(D1ref, D1,       'Error refined D    ');

nxD2 = nxD2_mat(obj, obj2, trian_pr);
nxD1 = full(nxD_mat(obj));
nxD1ref = nxD2*R;
errnxD = comp_mat(nxD1ref, nxD1, 'Error refined nxD  ');

Db = D_mat(obj2);
Dbref = R.'*Db*R;   % Funciona
errDb = comp_mat(Dbref, D1,       'Error refined Db   ');

nxDb = nxD_mat(obj2);
nxDbref = R.'*nxDb*R;   % Funciona
errnxDb = comp_mat(nxDbref, nxD1, 'Error refined nxDb ');

%% Test error in exact current
if strcmp(geom,'sphere')
    D3 = D_mat(obj2);
    Jex1 = j_ex(param.R,D, obj ,k,eta); % Exact current in original mesh
    user_plot_obj_current(Jex1,obj, EM_data);
    title('Jex1 exact original')
    
    Jex2 = j_ex(param.R,D3,obj2,k,eta); % Exact current in refined mesh
    user_plot_obj_current(Jex2,obj2, EM_data);
    title('Jex2 exact refined')
    
    Jex1ref = piR*Jex2;
    user_plot_obj_current(Jex1ref,obj, EM_data);
    title('Jex1re LC of refined exact Jex2')
    
    comp_mat(Jex1ref, Jex1, 'Error Jex1re LC of refined exact Jex2      ');
end

%% Test error in Ze
EM_data.field = 1;
Ze2 = user_impedance(1:N2, 1:N2, obj2, EM_data); 
Zeref = R.'*Ze2*R;   % Funciona

Jerefmat = -Zeref\Ei;
comp_mat(Jerefmat, Je,      'Error Je computed from LC of refined matrix');

user_plot_obj_current(Jerefmat,obj, EM_data);
title('Je computed from LC of refined matrix')

%% Test results of refined mesh
[Ei2, Hi2] = test_fields(obj2, th_i, ph_i, rot_EH, k, eta);
Je2 = -Ze2\Ei2;
user_plot_obj_current(Je2,obj2, EM_data);
title('Je2 computed refined')

Jeref = piR*Je2;
comp_mat(Jeref, Je,          'Error rJe from LC refined Je2              ');

user_plot_obj_current(Jeref,obj, EM_data);
title('Je from LC refined Je2')

errZ = comp_mat(Zeref, Ze, 'Error refined Z');

%% Test pseudo inverse
pinxD2 = pinv(nxD2);
GGp = nxD2*pinxD2;
GGp(GGp<1e-14) = 0;
errGGp = comp_mat(GGp, eye(N), 'Error GG^{+}-I ');

% GpG = pinxD2*nxD2;
% GpG(GpG<1e-14) = 0;
% errGpG = comp_mat(GpG, eye(N*6));
% fprintf('Error G^{+}G-I = %.3e\n', errGpG);
% 
% RRp = R*piR;
% RRp(RRp<1e-14) = 0;
% errRRp = comp_mat(RRp, eye(N*6));
% fprintf('Error RR^{+}-I = %.3e\n', errRRp);

RpR = piR*R;
RpR(RpR<1e-14) = 0;
errRpR = comp_mat(RpR, eye(N), 'Error R^{+}R-I ');



