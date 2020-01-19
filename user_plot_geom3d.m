% h = user_plot_geom3d(OG_data)
% Draws boundary of object
%
% Input:
% OG_data   = struct containing Object Geometry data
%
% Output:
% h         = handle to figure surface patches
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

function ret = user_plot_geom3d(OG_data)

figure;
color='g';

Nt = size(OG_data.topol,2);
polx = reshape(OG_data.vertex(1,OG_data.topol(:)),3,Nt);
poly = reshape(OG_data.vertex(2,OG_data.topol(:)),3,Nt);
polz = reshape(OG_data.vertex(3,OG_data.topol(:)),3,Nt);

if isfield(OG_data,'feed') && ~isempty(OG_data.feed)
    colormap([0 1 0; 1 0 0]);
    T1 = OG_data.edges(1,OG_data.feed);
    T2 = OG_data.edges(2,OG_data.feed);
    color = ones(1,Nt);
    color([T1 T2]) = 2;
end

ret = fill3(polx,poly,polz,color);

mincoor=1e10; maxcoor=-1e10;
mincoor = min([mincoor; OG_data.vertex(:)]); maxcoor = max([maxcoor; OG_data.vertex(:)]);

axis([mincoor maxcoor mincoor maxcoor mincoor maxcoor]);

%plot_feeding;

box;
xlabel('x'); ylabel('y'); zlabel('z');
axis equal
grid
set(gcf,'renderer','zbuffer')

end