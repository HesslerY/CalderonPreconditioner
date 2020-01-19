% [RCS_e, RCS_h] = bistat(R,lam,th_e,th_h)
%
% Bistatic RCS of a PEC sphere, from Balanis advanced EM, pp. 651..
% Incidence:	propagating from th=0 (+z) to th=-180 (-z)
%		polarized E=Ex, H=Hy
%
% R = radius in meters
% lam = wavelength in meters
% th_e = theta in E plane (degrees) Can be a vector
% th_h = theta in H plane (degrees) Can be a vector
%
% RCS_e = observation en E plane (xz) (column vector) in square meters
% RCS_h = observation en H plane (yz) (column vector) in square meters
%
% Juan M. Rius, november 1996

function [RCS_e, RCS_h] = bistat(a,lam,th_e,th_h)

k = 2*pi/lam;
ka = k*a;
res = 10000;
n = 1;

M_e = length(th_e); M_h = length(th_h);
sth_e = reshape(sin(th_e*pi/180),1,M_e); % Must be row vector for Legendre
sth_h = reshape(sin(th_h*pi/180),1,M_h);
cth_e = reshape(cos((180-th_e)*pi/180),1,M_e);
cth_h = reshape(cos((180-th_h)*pi/180),1,M_h);

term_e = ones(size(sth_e));
term_h = ones(size(sth_h));
E_e = zeros(size(sth_e));
E_h = zeros(size(sth_h));

while any(abs(term_e) > abs(E_e)/res) || any(abs(term_h) > abs(E_h)/res)

	an = (2*n+1) / (n*(n+1));
	bn = -an * jep(n,ka) ./ hep(n,ka);
	cn = -an * je(n,ka) ./ he(n,ka);

	Pn1cth_e  = pn1(n,cth_e); Pn1pcth_e = pn1p(n,cth_e);
	Pn1cth_h  = pn1(n,cth_h); Pn1pcth_h = pn1p(n,cth_h);

	term_e = (  bn * sth_e .* Pn1pcth_e - cn * Pn1cth_e ./ sth_e );	
	term_h = ( -cn * sth_h .* Pn1pcth_h + bn * Pn1cth_h ./ sth_h );

	E_e = E_e + term_e;
	E_h = E_h + term_h;
	n = n+1;
end

RCS_e = 4*pi*abs(E_e'/k).^2;	% Column vector
RCS_h = 4*pi*abs(E_h'/k).^2;	% Column vector

view_var('Number of terms in series summation',n);



