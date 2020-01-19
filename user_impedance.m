% Z = user_impedance(m, n, OG_data, EM_data)
%
% Impedance submatrix elements
% Green's function is EFIE, MFIE or CFIE 3-D
%
%
% m         = indices of testing functions
% n         = indices of basis functions
% OG_data   = struct containing Object Geometry data 
% EM_data   = struct containing ElectroMagnetic data 
%
% Fields of EM_data struct: all fields MUST BE DOUBLE data type
%    field		= 1->Ze (EFIE), 2->Zh (MFIE), 3->Ze+eta*Zh (CFIE)
%    k			= Wave number 2*pi/lam
%    eta			= Wave impedance of the medium
%    Rint_s		= Radius of 4-point source-integration.
%  	   			  For triangle centers farther that Rinteg, 1-point integration
%    Ranal_s      = Radius of analytical source-integration of the 1/R term in the Kernel
%    Rint_f       = Radius of 4-point field-integration              
%    solid_angle  = 0 : no solid angle correction in MFIE
%                   1 : solid angle correction in MFIE
%    flag         = 0: Compute scattered field = principal value + Js*solid_angle/(4*pi)
%                   1: Compute MFIE matrix = principal value + Js*solid_angle/(4*pi) - Js
%
% Juan M. Rius, E. Ubeda, AntennaLab, Universitat Politecnica de Catalunya (Spain), v1.0, August 2007

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

function Z = user_impedance(m, n, OG_data, EM_data)

% Computation of MoM impedance matrices
% Basis and test functions R,W,G, Galerkin
% Integration on basis functions:	Analytical if Ts==Tf
%                                   4 points (Gauss) if R < Rint
%                                   1 point at cent if R > Rint
% Integration on testing functions:	3 points (Gauss) if Tf==Ts
%                                   1 point at cent if Ts ~= Tf
%
% Output:
% Z = Impedance matrices for E or H radiated by J
%
% Input: all input arguments MUST BE DOUBLE data type
% f_elem		= Vector with indices of rows of Ze, Zh (test. func.) to be computed
% s_elem		= Vector with indices of cols of Ze, Zh (basis func.) to be computed
% object_f 	= Field testing object (=[ ] if same than source object)
% object_s 	= Source object
%   			 	If object_f = object_s, set obj_f=[] to compute correctly self-impedances
%				 	For information on object geometry, see object.m
% field		= 1->Ze (EFIE), 2->Zh (MFIE), 3->Ze+eta*Zh (CFIE)
% k			= Wave number 2*pi/lam
% eta			= Wave impedance of the medium
% Rint_s		= Radius of 4-point source-integration.
%	   			  For triangle centers farther that Rinteg, 1-point integration
% Ranal_s       = Radius of analytical source-integration of the 1/R term in the Kernel
% Rint_f        = Radius of 4-point field-integration              
% solid_angle   = 0 : no solid angle correction in MFIE
%                 1 : solid angle correction in MFIE
% flag      = 0 (default): Compute scattered field = principal value + Js*solid_angle/(4*pi)
%             1: Compute MFIE matrix = principal value + Js*solid_angle/(4*pi) - Js

field       = EM_data.field;
k           = EM_data.k;
eta         = EM_data.eta;
Rint_s      = EM_data.Rint_s;
Rint_f      = EM_data.Rint_f;
Ranal_s     = EM_data.Ranal_s;
corr_solid  = EM_data.corr_solid;
flag        = EM_data.flag;

Z = oper_3d_2_free(m,n,[ ],OG_data,field,k,eta,Rint_s,Ranal_s,Rint_f,corr_solid,flag);




