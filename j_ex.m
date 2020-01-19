% J = j_ex(a,D,obj,k,eta)
% Exact current for spheres expressed in RWG basis functions,
% incident plane wave, incidence direction (theta=180,phi=0),
% rotation of E with respect to phi = 90.
% REF: R.F. Harrington. Time-harmonic electromagnetic fields.
% McGraw-Hill,1971, pp. 293-295.
%
%
% Input:
% a	= sphere radius
% topol = topology matrix (vertex of each triangle), 3 x Nt
% vertex = vertex matrix, 3 x Nv
% trian = triangles matrix. For each triangle (column):
%		Row 1: Edge 1 (opposite is vertex 1)
%		Row 2: Edge 2 (opposite is vertex 2)
%		Row 3: Edge 3 (opposite is vertex 3)
%		If >0, T+ for that edge; if <0, T- for that edge
% cent 	= Centroid of each triangle,	3 x Nt
% un	= Unit normal to each triangle,	3 x Nt
% edges = Edges matrix, 4 x Ne. For each edge (column):
% 		Row 1: Triangle T+
%		Row 2: Triangle T-
%		Row 3: Global number of opposite vertex in T+
%		Row 4: Global number of opposite vertex in T-
% ln	 = Length of edges, 1 x Ne
% k	 = Wave number 2*pi/lam
% eta	 = Wave impedance of the medium
%
% IE-MEI v3.0, Josep Parr�n, March 1997
%
% WARNING: be careful with legendre function, if you are working with earlier version than
%	11-october-1996, you must change Pn1=Pn(2,:) by Pn1=-Pn(2,:) in order to obtain right
%	results.
% 	
%		

function J = j_ex(a,D,obj,k,eta)

topol   = obj.topol;
vertex  = obj.vertex;
trian   = obj.trian;
cent    = obj.cent;
un      = obj.un;
edges   = obj.edges;
ln      = obj.ln;

Tp=edges(1,:); Tm = edges(2,:);			% T+ and T- triangles corresponding to vertex
rp=cent(:,Tp) - vertex(:,edges(3,:));		% Rho of center in T+
rm=vertex(:,edges(4,:)) - cent(:,Tm);		% Rho of center in T-

d1=(cent(1,:).^2+cent(2,:).^2).^.5;			% Transform centroid to spherical coordinates
d2=(d1.^2+cent(3,:).^2).^.5;				%
cosph=cent(1,:)./d1; 					% 
sinph=cent(2,:)./d1;					%
costh=cent(3,:)./d2;					%
sinth=d1./d2;						%

res=10000;							% Resolution to stop iterations
n=1;
Jth=zeros(size(d1));termth=ones(size(d1));
Jph=zeros(size(d1));termph=ones(size(d1));

H2n1=sqrt(pi*k*a/2)*besselh(0.5,2,k*a);

while any(abs(termth) > abs(Jth)/res) & any(abs(termph) > abs(Jph)/res),

	H2n =sqrt(pi*k*a/2)*besselh(n+0.5,2,k*a);	% Hankel function order n+0.5
	H2nd=H2n1-(n*H2n/(k*a));			% Hankel function order n+0.5 derivated

	Pn=legendre(n,costh);				% Legendre function degree n order 1:n
%BE CAREFUL!!!!!
%	Pn1=-Pn(2,:);					% Legendre function degree n order 1 <11-10-96
	Pn1=Pn(2,:);					% Legendre function degree n order 1 >11-10-96
	if n==1,						% Legendre funciton degree n order 2
		Pn2=zeros(size(Pn1));			%
	else							%
		Pn2=Pn(3,:);				%
	end;							%

   Pn1sinth=(Pn1./sinth);				% Legendre function degree n order 1
                        % divided by sinth
	sinthPn1d=(-Pn1.*(costh./sinth))-Pn2;	% Legendre function degree n order 1 derived
								% multiplied by sinth

	an=(j)^(-n)*(2*n+1)/(n*(n+1));		% n-coefficient

	termth=an*(((sinthPn1d)/H2nd)+j*(Pn1sinth/H2n));
	termph=an*((Pn1sinth/H2nd)+j*((sinthPn1d)/H2n));
	Jth=Jth+termth;
	Jph=Jph+termph;

	n=n+1;
	H2n1=H2n;
end

disp('Number of terms computed:');disp(n);

Jth=j*cosph.*Jth/(eta*k*a);
Jph=j*sinph.*Jph/(eta*k*a);

Jla(1,:)=Jth.*costh.*cosph-Jph.*sinph;		% Transform current in spherical 
Jla(2,:)=Jth.*costh.*sinph+Jph.*cosph;		% coordinates to rectangular coordinates
Jla(3,:)=-Jth.*sinth;					%


Jl=-ln.*sum((Jla(:,Tp).*rp)+(Jla(:,Tm).*rm))/2; % Obtain RWG triangle coefficients 
J=D\(Jl.');							%

