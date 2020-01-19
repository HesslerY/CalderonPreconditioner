% user_plot_obj_current(J,OG_data,EM_data)
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

function user_plot_obj_current(J,OG_data,EM_data)

%% BEGIN Parameters

if OG_data.N < 5000,    arrows = 1;
else                    arrows = 0; %#ok<*SEPEX>
end

dB = 1;
%% END Parameters

figure

obj = OG_data;
eta = EM_data.eta;

mincoor=1e10; maxcoor=-1e10;
mincoor = min([mincoor; obj.vertex(:)]); maxcoor = max([maxcoor; obj.vertex(:)]);

Nv = length(obj.vertex);
Nt = length(obj.trian);

Jv = zeros(3,Nv);	% Jx,Jy,Jz for vertex 1:Nv
Jc = zeros(3,Nt);	% J at center of triangles
count = zeros(1,Nv);	% Number of edges for each vertex

for T = 1:Nt			% For all triangles

    for lv = [1 2 3]		% For all local vertex
        gv  = obj.topol(lv,T);	% Global vertex analyzed
        

        for le = find([1 2 3]~= lv)	% The local edges in this vertex
            ov = obj.topol(le,T);	% Opposite vertex to edge
            ge = obj.trian(le,T);	% Global edges
            if ge			% If the edge is not boundary
                si = sign(ge);	% Sign of triangle for this vertex
                ge = abs(ge);
                rho = si * (obj.vertex(:,gv)-obj.vertex(:,ov));	% rho vector
                Jt =  ( J(ge)*obj.ln(ge) / (2*obj.ds(T)) ) * rho;
                Jv(:,gv) = Jv(:,gv) + Jt;	% Add to total current
                count(gv) = count(gv)+1;
            end
        end

        % Compute current at center of triangle
        ge = obj.trian(lv,T);		% Edge opposite
        if ge				% If the edge is not boundary
            si = sign(ge);
            ge = abs(ge);
            rho = si * (obj.cent(:,T) - obj.vertex(:,gv));
            Jc(:,T) = Jc(:,T) + ( J(ge)*obj.ln(ge) / (2*obj.ds(T)) ) * rho;
        end
    end % For lv
end	% For T

count(count==0)=1;

% Normalize by count/2 (each edge has 2 vertex)
Jabs = eta * sqrt(sum( abs(Jv).^2 )) ./ (count/2);

maxJabs = max(Jabs); maxJ = maxJabs;
%minJabs = min(Jabs(find(Jabs)));
minJabs = min(Jabs(Jabs~=0));
minJ = minJabs;

if dB
    % minJ = maxJ/10000;

    Jabs(Jabs==0) = minJ;   % if there are zeros, set them to minimum value
    Jabs=20*log10(Jabs/maxJ);
end

polx = reshape(obj.vertex(1,obj.topol(:)),3,Nt);	% x coordinates, 3xNt
poly = reshape(obj.vertex(2,obj.topol(:)),3,Nt);	% y coordinates, 3xNt
polz = reshape(obj.vertex(3,obj.topol(:)),3,Nt);	% z coordinates, 3xNt

Cabs = reshape(Jabs(obj.topol(:)),3,Nt);		% Color Jabs, 3xNt

fill3(polx,poly,polz,Cabs);

%colormap(jet);
%set(gca,'color',[0.8 0.8 0.8]);
%xlabel('x'); ylabel('y'); zlabel('z');

%if feeding(1:5) == 'plane',
%    if dB, title('20log(|J/Hi|)');else title('|J/Hi|');end;
%else
%    if dB, title('20log(|J/Jmax|)');else title('|J|');end;
%end;
%colorbar;
%axis([mincoor maxcoor mincoor maxcoor mincoor maxcoor]); axis('equal');
%ca = caxis;

shading('interp');

if arrows
    % Plot current vectors
    % First scale Jc
    maxlen = mean(obj.ln)/3;			% Maximum length
    for T=1:Nt
        % Remove phase of Jc(T) relative to larger component of Jc(T)
%        imax = Jc(:,T) == max(Jc(:,T));
        tmp = find( Jc(:,T)==max(Jc(:,T)) ); 
        imax = tmp(1);
        Jc(:,T) = Jc(:,T) / sign(Jc(imax,T));
    end
    Jc = real(Jc);			% Remove circular polarization components
    nJc = sqrt(sum(abs(Jc).^2));	% Norm of Jc vectors
    Jc = Jc * maxlen / max(nJc);	% Scale maximum length to maxlen

    % Vectors to plot3
    cent = obj.cent + obj.un * maxlen/20;	% Add 'un' to plot just outside surface
    x = cent(1,:);
    y = cent(2,:);
    z = cent(3,:);
    u = Jc(1,:);	v = Jc(2,:);	w = Jc(3,:);
    xu = [x-u/2; x+u/2; NaN*ones(1,Nt)];
    yv = [y-v/2; y+v/2; NaN*ones(1,Nt)];
    zw = [z-w/2; z+w/2; NaN*ones(1,Nt)];

    hold on
    plot3(xu(:),yv(:),zw(:),'w');

    % Make arrow heads and plot them

    % Relative phase between different points is lost
    % -> vector sign has no meaning -> don't draw arrow heads

    % alpha = 1/3; % Size of arrow head relative to the length of the vector
    % beta = 0.5;  % Width of the base of the arrow head relative to the length
    %
    % norJ = cross(Jc,obj.un);	% Direction normal to Jc and tangent to surface
    %
    % ax = [x+u/2 - alpha*(u + beta*(norJ(1,:))) + obj.un(1,:)*maxlen/20 + eps; ...
    %       x+u/2; ...
    %       x+u/2 - alpha*(u - beta*(norJ(1,:))) + obj.un(1,:)*maxlen/20 + eps; ...
    % 	NaN*ones(1,Nt) ];
    %
    % ay = [y+v/2 - alpha*(v + beta*(norJ(2,:))) + obj.un(2,:)*maxlen/20 + eps; ...
    % 	y+v/2; ...
    %       y+v/2 - alpha*(v - beta*(norJ(2,:))) + obj.un(2,:)*maxlen/20 + eps; ...
    % 	NaN*ones(1,Nt) ];
    %
    % az = [z+w/2 - alpha*(w + beta*(norJ(3,:))) + obj.un(3,:)*maxlen/20 + eps; ...
    % 	z+w/2; ...
    %       z+w/2 - alpha*(w - beta*(norJ(3,:))) + obj.un(3,:)*maxlen/20 + eps; ...
    % 	NaN*ones(1,Nt) ];
    %
    % plot3(ax(:),ay(:),az(:),'w');

    hold off
end	% if arrows

colormap(jet);
set(gca,'color',[0.8 0.8 0.8]);
xlabel('x'); ylabel('y'); zlabel('z');

if dB, title('20log(|J/Hi|)'); 
else title('|J/Hi|');
end

axis([mincoor maxcoor mincoor maxcoor mincoor maxcoor]); axis('equal');

mdin = 20*log10(minJ/maxJ);
if mdin < -60
    fprintf('Dynamic margin in Js was %g dB',-mdin);
    disp('It has been changed to 60 dB');
    disp('If you want to display another dynamic margin, use ''caxis[-dyn_mar 0]; colorbar;''');
    mdin = -60; 
end

if dB, caxis( [mdin 0] ); 
else caxis([0 maxJ]); 
end

colorbar
set(gcf,'renderer','zbuffer')
