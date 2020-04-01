% Calderon precondiioner with refined mesh

disp(' ')

%% Problem data

N = size(Ze,1);
EM_data.field = 1;

% Refined mesh and associated objects
[obj2, trian_ch, trian_pr] = refine_mesh6(obj);
obj2 = get_edge(obj2);
obj2.N = length(obj2.ln);
N2 = obj2.N;
%user_plot_geom3d(obj2); 
%plot_obj_numbers(obj2);

%% Preliminary definitions
% In this section the matrices for the rest of the code are defined. This
% includes Gram matrices, inclusion operators and global operators

Ze2 = user_impedance(1:N2, 1:N2, obj2, EM_data); % Impedance matrix

% Inclusion and change of basis operators
R = linear_comb6(obj, obj2, trian_ch).'; % Inclusion matrix of RWG in RWGb
P = linear_comb_bc(obj, obj2, trian_ch).'; % Inclusion matrix of BC in RWGb

% Gram matrices
G_f_f = D_mat(obj); % Base RWG, test RWG %% Alternative definition via R
G_fb_fb = D_mat(obj2); % Base RWGb, test RWGb
G_nxf_f = nxD_mat(obj); % Base RWG, test nxRWG
G_nxfb_fb = nxD_mat(obj2); % Base RWGb, test nxRWGb
G_nxf_fb = nxD2_mat(obj, obj2, trian_pr); % Base RWGb, test nxRWG
G_bc_bc = P.'*G_fb_fb*P; % Base BC, test BC
G_nxf_bc = R.'*G_nxfb_fb*P; % Base BC, test nxRWG
G_y1_y1 = [G_fb_fb, G_nxfb_fb'; G_nxfb_fb, G_fb_fb]; % Base Y1, test Y1
G_y2_y2 = [G_fb_fb, G_nxf_fb'; G_nxf_fb, G_f_f]; % Base Y2, test Y2

% Matrices for connection operators and Z operators
Z_nxf_f = inv(G_f_f)*R.'*Ze2*R; % T operator + projection on nxRWG
Z_nxbc_bc = inv(G_bc_bc)*P.'*Ze2*P; % T operator + projection on nxBC
Z_nxfb_fb = inv(G_fb_fb)*Ze2; % T operator + projection on nxRWGb
C_nxf_bc = inv(G_f_f)*G_nxf_bc; % Projection from BC to nxRWG
C_fb_nxf = inv(G_fb_fb)*G_nxf_fb'; % Projection from nxRWG to RWGb
C_pi = G_nxf_fb\G_f_f; % Pseudoinverse connection operator

%% Matrix chains for different procedures and condition numbers
% The total chain for a MoM EFIE procedure is defined. Different numbers of
% condition are also written

fprintf('------------------------------------------------\nCONDITION NUMBERS FOR EACH APPROACH \n')

O_andr = Z_nxbc_bc*inv(C_nxf_bc)*Z_nxf_f; % Andriulli method
O_jmr1 = Z_nxfb_fb*C_fb_nxf*Z_nxf_f; % JMR method with orthogonal projection
O_jmr2 = C_pi.'*Ze2*C_pi*Z_nxf_f; % JMR method with pseudoinverse

fprintf('Condition number for JMR projection method:    %.5e,\n', cond(O_jmr1));
fprintf('Condition number for JMR pseudoinverse method: %.5e,\n', cond(O_jmr2));
fprintf('Condition number for Andriulli method:         %.5e.\n\n', cond(O_andr));

%% Norm of the error vectors in JM and Andriulli methods.
% Method exposed in the article. The Gram matrix for vector space Y is
% defined and so are the resulting vectors using each approximation. Notice
% that the definition of the difference of vectors uses thar the resulting
% vector of each approximation lies in RWGb and the result of exact
% computation lies in nxRWGb, adn both spaces are disjoint.

% Define the original vector and projections for every approximation
S_e = Z_nxf_f*Je;
S_jmr1 = C_fb_nxf*S_e;
S_jmr2 = C_pi*S_e;
S_andr = P*inv(C_nxf_bc)*S_e;

% Derivation by the first superspace
fprintf('------------------------------------------------\nFIRST VERSION: EXTENDED BOTH BASES \n')

% Difference vectors on superspace Y1
D1_jmr1 = [-S_jmr1; R*S_e];
D1_jmr2 = [-S_jmr2; R*S_e];
D1_andr = [-S_andr; R*S_e];

% Error norm calculation on superspace Y1
fprintf('Error with JMR projection method:    %.5e,\n', sqrt(D1_jmr1.'*G_y1_y1*D1_jmr1));
fprintf('Error with JMR pseudoinverse method: %.5e,\n', sqrt(D1_jmr2.'*G_y1_y1*D1_jmr2));
fprintf('Error with Andriulli approach:       %.5e.\n\n', sqrt(D1_andr.'*G_y1_y1*D1_andr));

% Derivation by the second superspace
fprintf('------------------------------------------------\nSECOND VERSION: EXTENDED ONLY DOMAIN BASE \n');

% Vectors on direct sum superspace and error norm calculation
D2_jmr1 = [-S_jmr1; S_e];
D2_jmr2 = [-S_jmr2; S_e];
D2_andr = [-S_andr; S_e];

% Error norm calculation on superspace Y2
fprintf('Error with JMR projection method:    %.5e,\n', sqrt(D2_jmr1.'*G_y2_y2*D2_jmr1));
fprintf('Error with JMR pseudoinverse method: %.5e,\n', sqrt(D2_jmr2.'*G_y2_y2*D2_jmr2));
fprintf('Error with Andriulli approach:       %.5e.\n\n', sqrt(D2_andr.'*G_y2_y2*D2_andr));
