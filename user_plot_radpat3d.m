% user_plot_radpat3d(J, OG_data, EM_data)
% Draws normalized radiation pattern in 3D
% User post-processing function
%
% OG_data   = struct containing Object Geometry data 
% EM_data   = struct containing ElectroMagnetic data 
%
% Juan M. Rius, AntennaLab, Universitat Politecnica de Catalunya (Spain), v1.0, August 2007

% /***************************************************************************
%  *   Copyright (C) 2007 by Juan M. Rius                                    *
%  *   AntennaLab, Universitat Politecnica de Catalunya, rius@tsc.upc.edu    *
%  *                                                                         *
%  *   This program is free software; you can redistribute it and/or modify  *
%  *   it under the terms of the GNU General Public License as published by  *
%  *   the Free Software Foundation; either version 2 of the License, or     *
%  *   (at your option) any later version.                                   *
%  *                                                                         *
%  *   This program is distributed in the hope that it will be useful,       *
%  *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
%  *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
%  *   GNU General Public License for more details.                          *
%  *                                                                         *
%  *   You should have received a copy of the GNU General Public License     *
%  *   along with this program; if not, write to the                         *
%  *   Free Software Foundation, Inc.,                                       *
%  *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
%  ***************************************************************************/

function user_plot_radpat3d(J, OG_data, EM_data)

%% BEGIN Parameters

dth = 10;   % Step in theta direction
dph = 10;   % Step in phi direction
mdin = 40;  % Dynamic range of 3D plot (dB)
dbcolor = 0;   % Flag (1/0). 1-> Show dB level in color scale, 0-> Lighted surface

%% END Parameters

% Draws radiaton pattern cuts phi=0, phi=90, theta=90
th = linspace(0,180,(180/dth)+1);
ph = linspace(0,360,(360/dph)+1);

% Compute radiation vector transverse components
Nth = [ ]; Nph = [ ];
for angth=1:length(th),
    for angph=1:length(ph),
        [ Nth(angph,angth), Nph(angph,angth) ] = r_pat(th(angth)*pi/180, ph(angph)*pi/180, OG_data, EM_data, J);
    end
end

Ntot = sqrt(abs(Nth).^2 + abs(Nph).^2);
N_max = max(abs(Ntot(:)));
Ntot = Ntot/N_max;

figure;

if dbcolor==0,
    d3patern(th,ph,Ntot,mdin,-37.5,30,'copper',0,0);
else
    d3patern(th,ph,Ntot,mdin,-37.5,30,'jet',0,1);
end;
title(sprintf('Normalized radiation pattern (dB), dynamic range = %d dB',mdin));
   
end

% [Nth,Nph] = r_pat(th, ph, OG_data,EM_data, J)
%
% Computes radiation vector transverse components for an observation angle
%
% Input:
% th, ph    = observation angle (theta, phi in radians)
% OG_data   = struct containing Object Geometry data 
% EM_data   = struct containing ElectroMagnetic data 
% J         = Normal current at RWG edges (unknonws in MoM)
%
% Output:
% Nth, Nph = radiation pattern 
%
% Josep Parron, December 1998

function [Nth,Nph] = r_pat(th, ph, OG_data, EM_data, J)

k = EM_data.k;

Tp = OG_data.edges(1,:); Tm = OG_data.edges(2,:);				% T+ and T- triangles corresponding to vertex
rho_p = OG_data.cent(:,Tp) - OG_data.vertex(:,OG_data.edges(3,:));	% Rho of center in T+
rho_m = OG_data.vertex(:,OG_data.edges(4,:)) - OG_data.cent(:,Tm);	% Rho of center in T-
un_r = [sin(th)*cos(ph) sin(th)*sin(ph) cos(th)];       % direction of observation
In = OG_data.ln.*(J.')/2;
N = rho_p*(In.*exp(j*k*un_r*OG_data.cent(:,Tp))).' + rho_m*(In.*exp(j*k*un_r*OG_data.cent(:,Tm))).';

Nth = N(1)*cos(th)*cos(ph) + N(2)*cos(th)*sin(ph) - N(3)*sin(th);
Nph = N(2)*cos(ph) - N(1)*sin(ph);

end

% function d3patern(th,ph,pattern,dB,az,el,mapacolor,rotacio,iluminacio)
%
% th:		vector fila amb els punts mesurats de theta
% ph:		vector fila amb els punts mesurats de phi
%
% pattern:	matriu amb els valors complexes de camp:
%		cada columna es correspon a un valor de theta,
%		cada fila es correspon a un valor de phi,
%
% dB:		marge dinamic per visualitzar el diagrama
%
% az:		perspectiva de visualitzacio en azimut		
%
% el:		perspectiva de visualitzacio en elevacio
%
% mapacolor:	especifica l'escala de color (vegeu 'colormap')
%
% rotacio:	rota el diagrama 90 graus al voltant de l'eix x en		
%		cas de que aquest parametre sigui a 1.
%
% iluminacio:	parametre que determina el tipus de visualitzacio:		
%		il.luminacio=0  : el diagrama s'il.lumina com si fos un objecte 3D
%		il.luminacio=1  : cada lobul s'il.lumina independentment segons el
%		nivell de potencia. Es recomana utilitzar el mapacolor 'copper' quan
%		la opcio escollida sigui 0, i el mapacolor 'jet' quan es vulgui que
%		el color indiqui el nivell de potencia dels lobuls
%
%
%		Carles Puente Baliarda,
%		D3-Enginyeria Electromagnetica i Fotonica,
%		UPC
%		12/1/1997
%
%		Modified by J.M. Rius, 8-1-98

function d3patern(th,ph,pattern,dB,az,el,mapacolor,rotacio,iluminacio)

% Normalize
r=20*log10(abs(pattern));
maxpat = 5*round(max(r(:))/5);
r = r-maxpat;

r(r<-dB)=-dB*ones(size(r(r<-dB)));
r=r+dB;
the=ones(size(1:size(pattern,1)))'*th/180*pi;
phi=(ones(size(1:size(pattern,2)))'*ph)'/180*pi;
x1=r.*sin(the).*cos(phi);
y1=r.*sin(the).*sin(phi);
z1=r.*cos(the);

if rotacio==1
y=y1;x=-z1;z=x1;
else
x=x1;y=y1;z=z1;
end

if iluminacio==0
	surfl(x,y,z),shading interp,colormap(mapacolor),axis('image'),view([az el]);
else
	surf(x,y,z,r),shading interp,colormap(mapacolor),axis('image'),view([az el]);
	hc = colorbar;
	yt = get(hc,'YTickLabel');
	set(hc,'YTickLabel',num2str(str2num(yt)+maxpat-dB));
end

xlabel('x'); ylabel('y'); zlabel('z');
set(gca,'xticklabel',[]); set(gca,'yticklabel',[]); set(gca,'Zticklabel',[])
set(gca,'ticklength',[0 0]);
box

end
