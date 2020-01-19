/*
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
*/

/***************************************************************************
 *   Copyright (C) 2007 by Juan M. Rius, Eduard Ubeda, Josep Parr�n        *
 *   AntennaLab, Universitat Politecnica de Catalunya, rius@tsc.upc.edu    *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
 
#include <math.h>
#include <omp.h>
#include "mex.h"

#define _PI     3.14159265358979

/* Arguments */
#define IN_ARG   4
#define p_rows      prhs[0]  /* m = indices of testing functions (arrays of double) */
#define p_cols      prhs[1]  /* n = indices of basis functions (arrays of double) */
#define p_OG_data   prhs[2]  /* OG_data   = struct containing Object Geometry data */
#define p_EM_data   prhs[3]  /* EM_data   = struct containing ElectroMagnetic data */

#define OUT_ARG  1
#define p_Z plhs[0] /* Impedance matrix */

#define norm(V)            sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]) 
#define rrdot(V1,V2)       (V1[0]*V2[0] + V1[1]*V2[1] + V1[2]*V2[2])
#define rcdot_r(V1,V2)     (V1[0]*V2[0].r + V1[1]*V2[1].r + V1[2]*V2[2].r)
#define rcdot_i(V1,V2)     (V1[0]*V2[0].i + V1[1]*V2[1].i + V1[2]*V2[2].i)
#define cross_rc(C,A,B)    { C[0].r = A[1]*B[2].r - A[2]*B[1].r; C[0].i = A[1]*B[2].i - A[2]*B[1].i;\
                             C[1].r = A[2]*B[0].r - A[0]*B[2].r; C[1].i = A[2]*B[0].i - A[0]*B[2].i;\
                             C[2].r = A[0]*B[1].r - A[1]*B[0].r; C[2].i = A[0]*B[1].i - A[1]*B[0].i;\
                           }
#define cross_rr(C,A,B)    { C[0] = A[1]*B[2] - A[2]*B[1];\
                             C[1] = A[2]*B[0] - A[0]*B[2];\
                             C[2] = A[0]*B[1] - A[1]*B[0];\
                           }


#define Gp 3

typedef double triad[3]; /* Vector or triangle information */
typedef double cuad[4];  /* Edge information */
typedef struct { double r,i; } dcomplex;

void autoz(dcomplex *pI, dcomplex Ii[3][3], int Ts, triad *topol, triad *vertex, triad *trian, double rt[3], triad *un, double k);
void anal_1divR(double *I_sc, double I_vc[3], int Ts, triad *topol, triad *vertex, triad *trian, double rt[3], triad *un, double *ds, double k, double signe );
int check_password(void);
void decode_password(unsigned int key, unsigned int p1, unsigned int p2, unsigned int *psys, unsigned int *psn, unsigned int *pdia, unsigned int *pmes, unsigned int *pano);
int compara_date(unsigned int dia, unsigned int mes, unsigned int ano, unsigned int dia2, unsigned int mes2, unsigned int ano2);

static int checked = 0;

/*extern "C" {*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ double *Zr, *Zi;
  /* Do not declare and initialise classes inside the loops !!!*/
  /*vector Rv, cTf, cTs, rt1, rt2, rt3, v0s, v1s, v2s, v0f, v1f, v2f;*/
  double Rv[3], Rn[3], cTf[3], cTs[3], rt1[3], rt2[3], rt3[3];
  double v0s[3], v1s[3], v2s[3], v0f[3], v1f[3], v2f[3];
  /*vector ovs, ovf, rho_s, rho_f, rho_1, rho_2, rho_3, rv;*/
  double ovs[3], ovf[3], rho_s[3], rho_f[3], rho_1[3], rho_2[3], rho_3[3], rv[3];
  double Rn_intf[3], unR_intf[3];
  double Rv_intf[3][3], rt_intf[3][3], A_intf[3][3];
  double rho_f_intf[3][3];
  /*cvector Iv, Ii1[3], Ii2[3], Ii3[3]; */ 
  dcomplex Iv[3], Ii1[3][3], Ii2[3][3], Ii3[3][3]; 
  /*cvector Ip, Ipe, Ipn, Ipz, nHs; */
  dcomplex Ip, Ipe, Ipn, Ipz, nHs[3], prod2[3];
  dcomplex  aux1[3], aux2[3], aux3[3], aux4[3], aux5[3], aux6[3], a1,a2,a3;
  dcomplex nHs_intf[3][3], prod2_intf[3][3], Iv_intf[3][3], aux1_intf[3][3];
  dcomplex Ip_intf[3], aux2_intf[3], aux3_intf[3];

  int se, si_s, ge_s, fe, si_f, ge_f, elem;
  double cteG, w0, w1;
  double r31[3], r32[3], cteD2[3], rho_3s[3], rho_3f[3], sum2[3];
  triad ren[3];
  double l31, l32, cteD1, sol1, sol2, sol3, solid, d1, d2, d3;
  int T1, T2, T3, e1, e2, e3;
  double *rows, *cols, *ds_s, *ln_s, *ds_f, *ln_f;
  int Nr, Nc, field, Ef, Hf, Nt_s, Ne_s, Nt_f, Ne_f;
  double k, eta, Rinteg_s, Rinteg_f, Ranal_s;
  triad *topol_s, *vertex_s, *trian_s, *un_s, *cent_s, *topol_f, *vertex_f, *trian_f, *un_f, *cent_f;
  cuad *edges_s, *edges_f;
  char *Tcomp_s, *Tcomp_f;
  double tmp, variabl;
  int flag, cor_solid;
  
  /* Fields of EM_data */
  mxArray *p_k, *p_eta, *p_field, *p_Rint_s, *p_Rint_f, *p_Ranal_s, *p_cor_solid, *p_flag;

  /* Fields of object */
  mxArray *p_topol_s, *p_vertex_s, *p_edges_s, *p_trian_s, *p_ln_s, *p_ds_s, *p_un_s, *p_cent_s;
  mxArray *p_topol_f, *p_vertex_f, *p_edges_f, *p_trian_f, *p_ln_f, *p_ds_f, *p_un_f, *p_cent_f;
  int singleobject;

  /* Basic type variables used inside the triangle loops. Faster declare here*/
  int r, c, i, iT_f, iT_s, T, Tf, Ts, nTf, nTs, integ_s, integ_f , j, p, anal_s;
  dcomplex G, Gc, Ge, Gh, Ae, Ah, Phi, I, Ie, In, Iz, I1, I2, I3;
  /* new for oper_3d_intf_anal */
  dcomplex Gc_intf[3], Ge_intf[3], Gh_intf[3];
  dcomplex I_intf[3];
  dcomplex iu;
  double R, kR, rhorho, sisi, cte, unR, unRho_s;
  double rhorho_intf[3];
  double R_intf[3], kR_intf[3], cte_intf;
  double rf1[3], rf2[3], rf3[3];
  double Ic_1divR_sc, I1_1divR_sc, I2_1divR_sc, I3_1divR_sc;
  double Ic_1divR_vc[3], I1_1divR_vc[3], I2_1divR_vc[3], I3_1divR_vc[3];
  double rintf_rv_prj_un, rc_rv_prj_un, rhof_A_c, rhof_A_intf;
  double vert_x_un, nrm_10_x_20;
  double rc_rv[3], rhoc_rhov[3], vec_edge_c[3], A_c[3];
  double rintf_rv[3], rhointf_rhov[3], vec_edge_intf[3], dif_vrt_1_0[3], dif_vrt_2_0[3];
  double vec_nrm_10_x_20[3], vec_10_x_20[3];
  int *Tlist_s, *Tlist_f, *Elist_s, *Elist_f;

  /* Integration*/
  iu.r=0;
  iu.i=1;
  /*triad ren[3] = { { 0.6, 0.2, 0.2 }, { 0.2, 0.6, 0.2 }, { 0.2, 0.2, 0.6 } };*/
  /* Integration points */
  /*triad ren[3]; */
  ren[0][0] = 0.6;
  ren[0][1] = 0.2;
  ren[0][2] = 0.2;
  ren[1][0] = 0.2;
  ren[1][1] = 0.6;
  ren[1][2] = 0.2;
  ren[2][0] = 0.2;
  ren[2][1] = 0.2;
  ren[2][2] = 0.6;
  w1 = 25.0/48;   /* Weight of symetric 3 points */
  w0 = -9.0/16;   /* Weight of centroid */

  /* Basic type variables used inside the edge loops. Faster declare here */
  /* int se, si_s, ge_s, fe, si_f, ge_f, elem_s, elem_f; */
  /* double cteG; */

  /* Computation of D/2 */
  /*double r31[3], r32[3], cteD2[3], rho_3s[3], rho_3f[3]; */
  /*double l31, l32, cteD1, sol1, sol2, sol3, solid, d1, d2, d3; */
  /*int T1, T2, T3, e1, e2, e3; */

  /* Check number of arguments */
  if(nrhs != IN_ARG)  { mexErrMsgTxt("Number of input arguments incorrect");}
  
  if(!mxIsStruct(p_OG_data)) { mexErrMsgTxt("OG_data is not struct"); }
  if(mxGetNumberOfElements(p_OG_data)!=1) { mexErrMsgTxt("OG_data is array with more than one element"); }
  
  if(!mxIsStruct(p_EM_data)) { mexErrMsgTxt("EM_data is not struct"); }
  if(mxGetNumberOfElements(p_EM_data)!=1) { mexErrMsgTxt("EM_data is array with more than one element"); }

  p_topol_s   = mxGetField(p_OG_data,0,"topol");
  p_vertex_s  = mxGetField(p_OG_data,0,"vertex");
  p_edges_s   = mxGetField(p_OG_data,0,"edges");
  p_trian_s   = mxGetField(p_OG_data,0,"trian");
  p_ln_s      = mxGetField(p_OG_data,0,"ln");
  p_ds_s      = mxGetField(p_OG_data,0,"ds");
  p_un_s      = mxGetField(p_OG_data,0,"un");
  p_cent_s    = mxGetField(p_OG_data,0,"cent");

  singleobject = 1;
  p_topol_f   = p_topol_s;
  p_vertex_f  = p_vertex_s;
  p_edges_f   = p_edges_s;
  p_trian_f   = p_trian_s;
  p_ln_f      = p_ln_s;
  p_ds_f      = p_ds_s;
  p_un_f      = p_un_s;
  p_cent_f    = p_cent_s;
  
  /* Asign pointers */
  rows  = mxGetPr(p_rows); Nr = mxGetN(p_rows); if(Nr==1) Nr = mxGetM(p_rows);
  cols  = mxGetPr(p_cols); Nc = mxGetN(p_cols); if(Nc==1) Nc = mxGetM(p_cols);
  
  p_k       	= mxGetField(p_EM_data,0,"k");
  p_eta     	= mxGetField(p_EM_data,0,"eta");
  p_field   	= mxGetField(p_EM_data,0,"field");
  p_Rint_s  	= mxGetField(p_EM_data,0,"Rint_s");
  p_Rint_f  	= mxGetField(p_EM_data,0,"Rint_f");
  p_Ranal_s  	= mxGetField(p_EM_data,0,"Ranal_s");
  p_cor_solid 	= mxGetField(p_EM_data,0,"corr_solid");
  p_flag	= mxGetField(p_EM_data,0,"flag");

  field    	= (int)*mxGetPr(p_field);
  k        	= *mxGetPr(p_k);
  eta      	= *mxGetPr(p_eta);
  Rinteg_s 	= *mxGetPr(p_Rint_s);
  Ranal_s  	= *mxGetPr(p_Ranal_s);
  Rinteg_f 	= *mxGetPr(p_Rint_f);
  cor_solid 	= (int)*mxGetPr(p_cor_solid);
  flag 		= (int)*mxGetPr(p_flag);
  
  topol_s  = (triad*) mxGetPr(p_topol_s);
  vertex_s = (triad*) mxGetPr(p_vertex_s);
  trian_s  = (triad*) mxGetPr(p_trian_s);  Nt_s = mxGetN(p_trian_s);
  un_s     = (triad*) mxGetPr(p_un_s);
  ds_s     = mxGetPr(p_ds_s);
  cent_s   = (triad*) mxGetPr(p_cent_s);
  edges_s  = (cuad*) mxGetPr(p_edges_s);
  ln_s     = mxGetPr(p_ln_s);              Ne_s = mxGetN(p_ln_s);

  topol_f  = (triad*) mxGetPr(p_topol_f);
  vertex_f = (triad*) mxGetPr(p_vertex_f);
  trian_f  = (triad*) mxGetPr(p_trian_f);  Nt_f = mxGetN(p_trian_f);
  un_f     = (triad*) mxGetPr(p_un_f);
  ds_f     = mxGetPr(p_ds_f);
  cent_f   = (triad*) mxGetPr(p_cent_f);
  edges_f  = (cuad*) mxGetPr(p_edges_f);
  ln_f     = mxGetPr(p_ln_f);              Ne_f = mxGetN(p_ln_f);
  
  Ef = (int)(field) % 2; /* Compute E field */
  Hf = field > 1;        /* Compute H field */

  p_Z = mxCreateDoubleMatrix(Nr,Nc,mxCOMPLEX);
  Zr  = mxGetPr(p_Z); Zi = mxGetPi(p_Z);

  /* Lists of triangles to compute*/
  Tlist_f = (int*) mxCalloc(2*Nr,sizeof(int)); /* Field*/
  Tlist_s = (int*) mxCalloc(2*Nc,sizeof(int)); /* Source*/

  /* 1->Tf computed, 0->Tf not computed*/
  Tcomp_f = (char*) mxCalloc(Nt_f,sizeof(char)); /* Field */
  Tcomp_s = (char*) mxCalloc(Nt_s,sizeof(char)); /* Source */

  /* Indices of edges in row. 0->Edge not in row*/
  Elist_f = (int*) mxCalloc(Ne_f,sizeof(int)); /* Field*/
  Elist_s = (int*) mxCalloc(Ne_s,sizeof(int)); /* Source*/

  for(i=0, nTf=0; i<Nr; i++)  /* nTf is current number of field triangles*/
  { r = rows[i]-1;
    T = edges_f[r][0]-1; if(!Tcomp_f[T]) { Tlist_f[nTf++] = T; Tcomp_f[T] = 1; }
    T = edges_f[r][1]-1; if(!Tcomp_f[T]) { Tlist_f[nTf++] = T; Tcomp_f[T] = 1; }
    Elist_f[r] = i+1;   /* Elist_f contains MATLAB indices: 0 -> Edge not in row*/
  } 

  for(i=0, nTs=0; i<Nc; i++)  /* nTs is current number of source triangles*/
  { c = cols[i]-1;
    T = edges_s[c][0]-1; if(!Tcomp_s[T]) { Tlist_s[nTs++] = T; Tcomp_s[T] = 1; }
    T = edges_s[c][1]-1; if(!Tcomp_s[T]) { Tlist_s[nTs++] = T; Tcomp_s[T] = 1; }
    Elist_s[c] = i+1;   /* Elist_s contains MATLAB indices: 0 -> Edge not in row*/
  }

 
#define INTERMED_VARS Ip, Ipe, Ipn, Ipz, a1, a2, a3, se, si_s, ge_s, fe, si_f, ge_f, elem, cteG, \
l31, l32, cteD1, sol1, sol2, sol3, solid, d1, d2, d3, T1, T2, T3, e1, e2, e3, tmp, \
r, c, i, iT_f, iT_s, T, Tf, Ts, integ_s, integ_f , j, p, anal_s, \
G, Gc, Ge, Gh, Ae, Ah, Phi, I, Ie, In, Iz, I1, I2, I3, \
R, kR, rhorho, sisi, cte, unR, unRho_s, cte_intf, Ic_1divR_sc, I1_1divR_sc, I2_1divR_sc, I3_1divR_sc, \
rintf_rv_prj_un, rc_rv_prj_un, rhof_A_c, rhof_A_intf, vert_x_un, nrm_10_x_20

#define INTERMED_ARRAYS Rv, Rn, cTf, cTs, rt1, rt2, rt3, v0s, v1s, v2s, v0f, v1f, v2f, \
ovs, ovf, rho_s, rho_f, rho_1, rho_2, rho_3, rv, Rn_intf, unR_intf, Rv_intf, rt_intf, A_intf, \
rho_f_intf, Iv, Ii1, Ii2, Ii3, nHs, prod2, aux1, aux2, aux3, nHs_intf, prod2_intf, Iv_intf, aux1_intf, \
Ip_intf, aux2_intf, aux3_intf, r31, r32, cteD2, rho_3s, rho_3f, sum2, \
Gc_intf, Ge_intf, Gh_intf, I_intf, rhorho_intf, R_intf, kR_intf, rf1, rf2, rf3, \
Ic_1divR_vc, I1_1divR_vc, I2_1divR_vc, I3_1divR_vc, rc_rv, rhoc_rhov, vec_edge_c, A_c, \
rintf_rv, rhointf_rhov, vec_edge_intf, dif_vrt_1_0, dif_vrt_2_0, vec_nrm_10_x_20, vec_10_x_20

/* Sets the number of threads in subsequent parallel regions equal to the number of processors that are available */
omp_set_num_threads(omp_get_num_procs());

  /* Computation */
#pragma omp parallel private(INTERMED_VARS, INTERMED_ARRAYS)
{

#pragma omp for schedule(dynamic) nowait
  for(iT_f=0; iT_f<nTf; iT_f++)
  { Tf = Tlist_f[iT_f];
    
  for(iT_s=0; iT_s<nTs; iT_s++)
  { Ts = Tlist_s[iT_s];

/*  First, compute all the factors that don't depend on edges*/
/*  if((Ts != Tf) || !singleobject)        /* No singular integration */
/*  { /* cTf = cent[Tf]; cTs = cent[Ts]; Rv = cTf-cTs; R = norm(Rv); */

    for(j=0; j<3; j++)
    { cTf[j] = cent_f[Tf][j]; cTs[j] = cent_s[Ts][j];
      Rv[j] = cTf[j] - cTs[j]; 
    }
    R = norm(Rv);
    
    if ((R > 1e-14) || !singleobject)	/* No self-term integration */
    { 
      integ_s = ( R < Rinteg_s );
      anal_s  = ( R < Ranal_s );
      integ_f = ( R < Rinteg_f );

      kR = k*R;

      Gc.i = -sin(kR)/(4*_PI*R); 

      if(Hf) 
      { /* Gh = Gc * complex(1,kR)/(R*R); */
        Gc.r = cos(kR)/(4*_PI*R);  	/* G at cent */

        Gh.r = (Gc.r - Gc.i*kR)/(R*R);
        Gh.i = (Gc.i + Gc.r*kR)/(R*R);

      } /* if (Hf) */

      if (Ef)
      {
	if (!anal_s)  
        { 
	   Gc.r = cos(kR)/(4*_PI*R); 
        }
	else 	
	{  

	   Gc.r = (cos(kR) - 1)/(4*_PI*R);  /* G at cent */

	   for(j=0; j<3; j++)
	   {
              dif_vrt_1_0[j] = vertex_s[ (int)topol_s[Ts][1]-1 ][j] - vertex_s[ (int)topol_s[Ts][0]-1 ][j];
              dif_vrt_2_0[j] = vertex_s[ (int)topol_s[Ts][2]-1 ][j] - vertex_s[ (int)topol_s[Ts][0]-1 ][j]; 
	   }

	   cross_rr(vec_10_x_20 , dif_vrt_1_0 , dif_vrt_2_0);
	   nrm_10_x_20 = norm(vec_10_x_20);
           for(j=0; j<3; j++)
           {
	      vec_nrm_10_x_20[j] = vec_10_x_20[j]/nrm_10_x_20;
	   }

	   vert_x_un = -rrdot( vec_nrm_10_x_20 , un_s[Ts]);

	} /* else if (!anal_s) */

      } /* if (Ef) */	

      if(!integ_s)  /* No integration (1 point at centroid) */
      { if(Ef) /*Ge = iu*eta * Gc;*/
        { 
 	  Ge.r = eta*(- iu.i*Gc.i);
	  Ge.i = eta*(  iu.i*Gc.r);

	  /* Analytical integration of the 1/R term for EFIE */
	  if(anal_s) 
	  {
	    anal_1divR(&Ic_1divR_sc, Ic_1divR_vc, Ts, topol_s, vertex_s, trian_s, cTf, un_s, ds_s, k, vert_x_un);
	  }

        } /* if (Ef) */

        if(Hf) unR = rrdot(un_f[Tf],Rv);
      }
      else        /* Integrate centroid + 3 other points */
      { 

	for(j=0; j<3; j++)
        {  v0s[j] = vertex_s[ (int)topol_s[Ts][0]-1 ][j];
           v1s[j] = vertex_s[ (int)topol_s[Ts][1]-1 ][j];
           v2s[j] = vertex_s[ (int)topol_s[Ts][2]-1 ][j];
        } 

        if(Hf)
        { 
          Ip.r = w0*Gh.r;
          Ip.i = w0*Gh.i;
          Ipe.r = Ip.r/3; Ipn.r = Ipe.r;
          Ipe.i = Ip.i/3; Ipn.i = Ipe.i;
        } /* if(Hf) */

        if(Ef) 
        { 
          I.r  = w0*Gc.r;
          I.i  = w0*Gc.i;

          Ie.r = I.r/3; In.r = Ie.r;
          Ie.i = I.i/3; In.i = Ie.i;
        }  /* if(Ef) */

        for(i=0; i<3; i++)  /* Integration loop: 3 points */
        { for(j=0; j<3; j++)
          { rv[j] = ren[i][0]*v0s[j] + ren[i][1]*v1s[j] + ren[i][2]*v2s[j];
            Rv[j] = cTf[j]-rv[j];
          }  

          R = norm(Rv); kR = k*R;
          G.i = w1*(-sin(kR)/(4*_PI*R)); 

          if(Ef)
          { 
	    if (!anal_s) G.r = w1*cos(kR)/(4*_PI*R); 
	    else  G.r = w1*(cos(kR) - 1)/(4*_PI*R); 

            I.r  += G.r; I.i  += G.i;
            Ie.r += G.r *ren[i][0]; Ie.i += G.i *ren[i][0];
            In.r += G.r *ren[i][1]; In.i += G.i *ren[i][1];
          } /* if(Ef) */

          if(Hf)
          {
	    G.r = w1*cos(kR)/(4*_PI*R); 
	
            Gh.r= (G.r - G.i*kR)/(R*R);
            Gh.i= (G.i + G.r*kR)/(R*R);

            Ip.r  += Gh.r;           Ip.i  += Gh.i;
            Ipe.r += Gh.r*ren[i][0]; Ipe.i += Gh.i*ren[i][0];
            Ipn.r += Gh.r*ren[i][1]; Ipn.i += Gh.i*ren[i][1];
          } /* if(Hf) */

        } /* for(i) integration loop */

        if(Ef)
        { 
	  /* Analytical integration of the 1/R term for EFIE */
	  if(anal_s) 
          {
	    anal_1divR(&Ic_1divR_sc, Ic_1divR_vc, Ts, topol_s, vertex_s, trian_s, cTf, un_s, ds_s, k, vert_x_un);
	  }

	  /* I = iu*eta*I; Ie = iu*eta*Ie; In = iu*eta*In; */
          tmp = -eta*I.i;  I.i  = eta*I.r;  I.r  = tmp;
          tmp = -eta*Ie.i; Ie.i = eta*Ie.r; Ie.r = tmp;
          tmp = -eta*In.i; In.i = eta*In.r; In.r = tmp;

          /*Iz = I - Ie - In;*/
          Iz.r = I.r - Ie.r - In.r;
          Iz.i = I.i - Ie.i - In.i;

          for(j=0; j<3; j++)
          { Iv[j].r = v0s[j]*Ie.r+ v1s[j]*In.r + v2s[j]*Iz.r; 
            Iv[j].i = v0s[j]*Ie.i+ v1s[j]*In.i + v2s[j]*Iz.i; 
          } 

        } /* if(Ef) */

        if(Hf)
        { /*Ipz = Ip - Ipe - Ipn;*/
          Ipz.r = Ip.r - Ipe.r - Ipn.r;
          Ipz.i = Ip.i - Ipe.i - Ipn.i;

          for(j=0; j<3; j++)
          { aux1[j].r = v0s[j]*Ipe.r + v1s[j]*Ipn.r + v2s[j]*Ipz.r;
            aux1[j].i = v0s[j]*Ipe.i + v1s[j]*Ipn.i + v2s[j]*Ipz.i;
          }

        } /* if (Hf) */

      }	/* else if (!integ_s) */

      if (integ_f)	/* field integration */
      {
         /* setups */
         for(j=0; j<3; j++)
         { 
	   v0f[j] = vertex_f[ (int)topol_f[Tf][0]-1 ][j];
           v1f[j] = vertex_f[ (int)topol_f[Tf][1]-1 ][j];
           v2f[j] = vertex_f[ (int)topol_f[Tf][2]-1 ][j]; 
 	 }

         for(j=0; j<3; j++)
         {
            for(i=0; i<3; i++)
            {
              rt_intf[j][i] = ren[j][0]*v0f[i] + ren[j][1]*v1f[i] + ren[j][2]*v2f[i]; 
              Rv_intf[j][i] = rt_intf[j][i] - cTs[i];
            }
            R_intf[j] = sqrt( pow(Rv_intf[j][0],2)+pow(Rv_intf[j][1],2)+pow(Rv_intf[j][2],2) ); 
            kR_intf[j] = k*R_intf[j];      

	    Gc_intf[j].i = -sin(kR_intf[j])/(4*_PI*R_intf[j]); 

 	    if (Ef)
	    {
	      if (!anal_s) 
	      {	
    	        Gc_intf[j].r =  cos(kR_intf[j])/(4*_PI*R_intf[j]); 
	      }
	      else
	      {	
	        Gc_intf[j].r =  (cos(kR_intf[j]) - 1)/(4*_PI*R_intf[j]); 
	      }
	    } /* if (Ef) */

            if(Hf) 
            {
 	      Gc_intf[j].r =  cos(kR_intf[j])/(4*_PI*R_intf[j]); 

	      Gh_intf[j].r = (Gc_intf[j].r - Gc_intf[j].i*kR_intf[j])/(R_intf[j]*R_intf[j]);
              Gh_intf[j].i = (Gc_intf[j].i + Gc_intf[j].r*kR_intf[j])/(R_intf[j]*R_intf[j]);

            } /* if (Hf) */

	  } /* for(j=0; j<3; j++) */

          if(!integ_s)  /* No integration (1 point at centroid) */
          {

	      if(Ef) 
              { 
		for (j=0; j<3; j++)
		{
		   Ge_intf[j].r = eta*(- iu.i*Gc_intf[j].i);
                   Ge_intf[j].i = eta*( iu.i*Gc_intf[j].r);
		}

		if (anal_s) 
		{

		   for (j=0; j<3; j++)
		   {	
		     rf1[j] = rt_intf[0][j];
  		     rf2[j] = rt_intf[1][j];
		     rf3[j] = rt_intf[2][j];	
		   }  

  		   anal_1divR( &I1_1divR_sc, I1_1divR_vc, Ts, topol_s, vertex_s, trian_s, rf1, un_s, ds_s, k , vert_x_un);
		   anal_1divR( &I2_1divR_sc, I2_1divR_vc, Ts, topol_s, vertex_s, trian_s, rf2, un_s, ds_s, k , vert_x_un);
		   anal_1divR( &I3_1divR_sc, I3_1divR_vc, Ts, topol_s, vertex_s, trian_s, rf3, un_s, ds_s, k , vert_x_un);

		}   /* if(anal_s) */

               }    /* if(Ef) */

               if(Hf) 
               { 
		 unR_intf[0] = rrdot(un_f[Tf],Rv_intf[0]);
                 unR_intf[1] = rrdot(un_f[Tf],Rv_intf[1]);
                 unR_intf[2] = rrdot(un_f[Tf],Rv_intf[2]);
               }

          }
          else
          {
            for(p=0; p<3; p++)   /* Integrate centroid + 3 other points for 3 field points */
            {  
	       if(Hf)
               { 
                 Ip_intf[p].r = w0*Gh_intf[p].r;
                 Ip_intf[p].i = w0*Gh_intf[p].i;
                 Ipe.r = Ip_intf[p].r/3; Ipn.r = Ipe.r;
                 Ipe.i = Ip_intf[p].i/3; Ipn.i = Ipe.i;
               } /* if (Hf) */

               if(Ef) 
               { 
                  I_intf[p].r  = w0*Gc_intf[p].r;
                  I_intf[p].i  = w0*Gc_intf[p].i;

                  Ie.r = I_intf[p].r/3; In.r = Ie.r;
                  Ie.i = I_intf[p].i/3; In.i = Ie.i;
               } /* if (Ef) */

               for(i=0; i<3; i++)  /* Integration loop: 3 points */
               { 
		   for(j=0; j<3; j++)
                   { 
		      rv[j] = ren[i][0]*v0s[j] + ren[i][1]*v1s[j] + ren[i][2]*v2s[j];
                      Rv[j] = rt_intf[p][j]-rv[j]; 
		   }  

                   R = norm(Rv); kR = k*R;
 		   G.i =  w1*(-sin(kR)/(4*_PI*R));		

                   if (Ef)				
		   {	
		      if (!anal_s )  G.r = w1*cos(kR)/(4*_PI*R);
 		      else 	  G.r = w1*(cos(kR) - 1)/(4*_PI*R);

                      I_intf[p].r += G.r; I_intf[p].i += G.i;
                      Ie.r += G.r *ren[i][0]; Ie.i += G.i *ren[i][0];
                      In.r += G.r *ren[i][1]; In.i += G.i *ren[i][1];

                   } /* if(Ef) */

                   if(Hf)
                   { 	
		      G.r = w1*cos(kR)/(4*_PI*R);

		      Gh.r= (G.r - G.i*kR)/(R*R);
                      Gh.i= (G.i + G.r*kR)/(R*R);

                      Ip_intf[p].r  += Gh.r;   Ip_intf[p].i  += Gh.i;
                      Ipe.r += Gh.r*ren[i][0]; Ipe.i += Gh.i*ren[i][0];
                      Ipn.r += Gh.r*ren[i][1]; Ipn.i += Gh.i*ren[i][1];

                   } /* if(Hf) */

               } /* for(i) integration loop */

	       /* Analytical integration of the 1/R term for EFIE */
	       if(Ef)
	       {

		   if (anal_s)
		   {

   		      for (j=0; j<3; j++)
		      {	
		       rf1[j] = rt_intf[0][j];
  		       rf2[j] = rt_intf[1][j];
		       rf3[j] = rt_intf[2][j];	
		      } 

  		      anal_1divR(&I1_1divR_sc, I1_1divR_vc, Ts, topol_s, vertex_s, trian_s, rf1, un_s, ds_s, k , vert_x_un);
		      anal_1divR(&I2_1divR_sc, I2_1divR_vc, Ts, topol_s, vertex_s, trian_s, rf2, un_s, ds_s, k , vert_x_un);
		      anal_1divR(&I3_1divR_sc, I3_1divR_vc, Ts, topol_s, vertex_s, trian_s, rf3, un_s, ds_s, k , vert_x_un);

		   }	/* if (anal_s) */

                   tmp = -eta*I_intf[p].i; I_intf[p].i  = eta*I_intf[p].r; I_intf[p].r  = tmp;
                   tmp = -eta*Ie.i; Ie.i = eta*Ie.r; Ie.r = tmp;
                   tmp = -eta*In.i; In.i = eta*In.r; In.r = tmp;
                                
                   Iz.r = I_intf[p].r - Ie.r - In.r;
                   Iz.i = I_intf[p].i - Ie.i - In.i;
        
                   for(j=0; j<3; j++)
                   {
		       Iv_intf[p][j].r = v0s[j]*Ie.r+ v1s[j]*In.r + v2s[j]*Iz.r; 
                       Iv_intf[p][j].i = v0s[j]*Ie.i+ v1s[j]*In.i + v2s[j]*Iz.i; 
                   }

               }   /* if(Ef) */

               if(Hf)
               {     
		   Ipz.r = Ip_intf[p].r - Ipe.r - Ipn.r;
                   Ipz.i = Ip_intf[p].i - Ipe.i - Ipn.i;

                   for(j=0; j<3; j++)
                   { 
                      aux1_intf[p][j].r = v0s[j]*Ipe.r + v1s[j]*Ipn.r + v2s[j]*Ipz.r;
                      aux1_intf[p][j].i = v0s[j]*Ipe.i + v1s[j]*Ipn.i + v2s[j]*Ipz.i;
                   }
               } /* if (Hf) */                                                                       

            } /* for(p=0; p<3; p++) */

         }  /*  else if (!integ_s)  integrate centroid + 3 points  */

       }  /* if integ_f */

    }
    else /* Self-triangle*/
    { /* Vertex of triangle*/
      for(j=0; j<3; j++)
      { 
        v0f[j] = vertex_f[ (int)topol_f[Tf][0]-1 ][j];
        v1f[j] = vertex_f[ (int)topol_f[Tf][1]-1 ][j];
        v2f[j] = vertex_f[ (int)topol_f[Tf][2]-1 ][j];
      } 
        
      if(Ef) 	/* Singular integration for EFIE*/
      { /* Test points*/
        for(j=0; j<3; j++)
        { /* Test points */
          rt1[j]= v0f[j]*2/3 + v1f[j]*1/6 + v2f[j]*1/6;
          rt2[j]= v0f[j]*1/6 + v1f[j]*2/3 + v2f[j]*1/6;
          rt3[j]= v0f[j]*1/6 + v1f[j]*1/6 + v2f[j]*2/3;
        }

        autoz(&I1,Ii1,Ts,topol_s,vertex_s,trian_s,rt1,un_s,k);
        autoz(&I2,Ii2,Ts,topol_s,vertex_s,trian_s,rt2,un_s,k);
        autoz(&I3,Ii3,Ts,topol_s,vertex_s,trian_s,rt3,un_s,k);

      } /* if (Ef) */

      if(Hf)  	/* MFIE: compute D/2*/
      { 
	for(j=0; j<3; j++)
        { 
          r31[j] = v0f[j] - v2f[j];
          r32[j] = v1f[j] - v2f[j];
          cteD2[j] = r31[j] + r32[j];
        }
        l31 = norm(r31);
        l32 = norm(r32);
        cteD1 = l31*l31 + l32*l32 + rrdot(r31,r32);
        
        solid = 1;	
        
        if (cor_solid)
        {        /* Compute solid angle for neighbor triangles*/
                e1 = abs((int)trian_s[Ts][0])-1;
                e2 = abs((int)trian_s[Ts][1])-1;
                e3 = abs((int)trian_s[Ts][2])-1;

                T1 = (int)edges_s[e1][0]-1; if(T1==Ts) T1 = (int)edges_s[e1][1]-1;
                T2 = (int)edges_s[e2][0]-1; if(T2==Ts) T2 = (int)edges_s[e2][1]-1;
                T3 = (int)edges_s[e3][0]-1; if(T3==Ts) T3 = (int)edges_s[e3][1]-1;

                d1 = rrdot(un_s[Ts],un_s[T1]); if(d1>1) d1=1; else if(d1<-1) d1=-1;
                d2 = rrdot(un_s[Ts],un_s[T2]); if(d2>1) d2=1; else if(d2<-1) d2=-1;
                d3 = rrdot(un_s[Ts],un_s[T3]); if(d3>1) d3=1; else if(d3<-1) d3=-1;

                sol1 = 1+acos(d1)/_PI; sol2 = 1+acos(d2)/_PI; sol3 = 1+acos(d3)/_PI;
                solid = (sol1 + sol2 + sol3)/3;
        }   /* if (cor_solid) */

        if(flag) solid -= 2;
      } /* if(Hf) */
    } /* if else self-triangle */
      
    for(se=0; se<3; se++)      /* Edges of source triangle*/
    { ge_s = trian_s[Ts][se];    /* Global edges of source triangles*/

      if(ge_s==0) continue;   /* Boundary edge, not interior edge*/
      si_s = (ge_s>0)? 1: -1;  /* Sign of T for this edges*/
      ge_s = abs(ge_s)-1;      /* Edge number in trian 1:Ne*/
      if(Elist_s[ge_s]==0) continue;

      /* Create rho_s vector for source edge, with +/- sign */
      for(j=0; j<3; j++)
      { 
        ovs[j] = vertex_s[ (int)topol_s[Ts][se]-1 ][j];   /* vertex of source triangle */
        rho_s[j] = si_s * (cTs[j] - ovs[j]);              /* rho vector for centroid*/
      }

      if(Ef) 
      {
	if(anal_s) 
	{
 	   for (j=0; j<3; j++)
           {
             rc_rv[j] = cTf[j] - ovs[j];
           }	

           rc_rv_prj_un = rrdot( rc_rv , un_s[Ts] );
              
           for (j=0; j<3; j++)
           {
             rhoc_rhov[j] = rc_rv[j] - rc_rv_prj_un*un_s[Ts][j];	 
           }

           for (j=0; j<3; j++)
           {
             vec_edge_c[j] = rhoc_rhov[j] * Ic_1divR_sc;
           }

           for (j=0; j<3; j++)
           {
             A_c[j] = vec_edge_c[j] + Ic_1divR_vc[j];
           }


	   if (integ_f) 
	   {

 	      for(p=0; p<3; p++)
	      {

		for (j=0; j<3; j++)
	        {
	          rintf_rv[j] = rt_intf[p][j] - ovs[j];
	        }	

                rintf_rv_prj_un = rrdot( rintf_rv , un_s[Ts] );

	        for (j=0; j<3; j++)
                {
                  rhointf_rhov[j] = rintf_rv[j] - rintf_rv_prj_un * un_s[Ts][j];
                }

	        for(j=0; j<3; j++)
	        {
	          if (p==0) vec_edge_intf[j] = rhointf_rhov[j] * I1_1divR_sc;
	          if (p==1) vec_edge_intf[j] = rhointf_rhov[j] * I2_1divR_sc;
	          if (p==2) vec_edge_intf[j] = rhointf_rhov[j] * I3_1divR_sc;
                } 

	        for(j=0; j<3; j++)
                { 
                  if (p==0) A_intf[p][j] = vec_edge_intf[j] + I1_1divR_vc[j];
		  if (p==1) A_intf[p][j] = vec_edge_intf[j] + I2_1divR_vc[j];
		  if (p==2) A_intf[p][j] = vec_edge_intf[j] + I3_1divR_vc[j];
                } 

	      }	/* for(p=0; p<3; p++) */

	   }  /* if (integ_f) */	

	}  /* if(anal_s) */	

      } /*if (Ef)*/


      if(Hf)
      { if(!integ_s) unRho_s = rrdot(un_f[Tf],rho_s);
        else
        { /* nHs = n x Hs */
          /* nHs = -si_s * cross(un_f[Tf], cross(Rn, (aux1 - ovs*Ip) )); */ 
          for(j=0; j<3; j++)
          { aux2[j].r = -si_s * (aux1[j].r - ovs[j]*Ip.r);
            aux2[j].i = -si_s * (aux1[j].i - ovs[j]*Ip.i);
            Rn[j] = cTf[j] - ovs[j]; /* Rn = cTf - ovs; */
          }
          cross_rc(aux3,Rn,aux2);
          cross_rc(nHs,un_f[Tf],aux3);

          if (integ_f)
          {
            for(p=0; p<3; p++)
            {  for(j=0; j<3; j++)
               { aux2_intf[j].r = -si_s * (aux1_intf[p][j].r - ovs[j]*Ip_intf[p].r);
                 aux2_intf[j].i = -si_s * (aux1_intf[p][j].i - ovs[j]*Ip_intf[p].i);
                 Rn_intf[j] = rt_intf[p][j] - ovs[j]; /* Rn = cTf - ovs; */
               }
             cross_rc(aux3_intf,Rn_intf,aux2_intf);
             cross_rc(nHs_intf[p],un_f[Tf],aux3_intf);
            } 
          } /* if (integ_f) */
        }   /* else if (!integ_s) */
      }  /* if (Hf) */

      for(fe=0; fe<3; fe++)     /* Edges of field triangle */
      { 
	ge_f = trian_f[Tf][fe];    /* Global edges of field triangles */

        if(ge_f==0) continue;   /* Boundary edge, not interior edge*/
        si_f = (ge_f>0)? 1: -1;
        ge_f = abs(ge_f)-1;

        if(Elist_f[ge_f]==0) continue;
        
        for(j=0; j<3; j++)
        { 
	  ovf[j] = vertex_f[ (int)topol_f[Tf][fe]-1 ][j]; /* vertex of field triangle */
        }           

        if((Ts != Tf) || !singleobject)  /* Green's function tested at centroid */
        { 
       	  for(j=0; j<3; j++)
          {
           rho_f[j] = si_f * (cTf[j] - ovf[j]); /* rho vector for centroid */
          }
          
          if (integ_f)
          {
            for(p=0; p<3; p++)
            {
              for(j=0; j<3; j++)
              {
               rho_f_intf[p][j] = si_f * (rt_intf[p][j] - ovf[j]);        /* rho vector for rt_intf*/
              }
            }
           } /* if integ_f */

           if(!integ_s)
           { 
	     rhorho = rrdot(rho_s,rho_f);

             if (integ_f)
             {
               for(j=0; j<3; j++)
               {
                 rhorho_intf[j] = rrdot(rho_s,rho_f_intf[j]);
               }
             } /* if integ_f */ 

             if(Ef)
             { 
	       sisi = si_s * si_f / k;
               Ae.r = -k * Ge.r * rhorho/4;
               Ae.i = -k * Ge.i * rhorho/4;
               Phi.r = Ge.r * sisi;
               Phi.i = Ge.i * sisi;

               /* Analytical integration of 1/R ; field point in the centroid */
	       if (anal_s) 
               {  
   	         rhof_A_c = rrdot(rho_f,A_c);

  	         Ae.i  += -0.25 * si_s *  k * eta * rhof_A_c /(4*_PI);  
                 Phi.i +=  eta * sisi * Ic_1divR_sc /(4*_PI) ; 
	       }  /* if (anal_s) */

               if (integ_f)
               {  
		   tmp = Ae.r;  Ae.r  = w0*tmp;
		   tmp = Ae.i;  Ae.i  = w0*tmp;
                   tmp = Phi.r; Phi.r = w0*tmp;
                   tmp = Phi.i; Phi.i = w0*tmp;

                   for(j=0; j<3; j++)
                   {
                     Ae.r +=  -k * w1 * Ge_intf[j].r * rhorho_intf[j]/4;
                     Ae.i +=  -k * w1 * Ge_intf[j].i * rhorho_intf[j]/4;
                     Phi.r +=  w1 * Ge_intf[j].r * sisi;
                     Phi.i +=  w1 * Ge_intf[j].i * sisi;
                   }

                   /* Analytical integration of 1/R ; + 3 field points */
    	      	   if (anal_s) 
		   {	

		     for (p=0; p<3; p++)	
                     { 
 	               rhof_A_intf = rrdot(rho_f_intf[p],A_intf[p]);

  	               Ae.i +=  -w1 * 0.25 * si_s * k * eta * rhof_A_intf / (4*_PI) ; 
                       if (p==0) Phi.i +=  w1 * eta * sisi * I1_1divR_sc / (4*_PI) ;
                       if (p==1) Phi.i +=  w1 * eta * sisi * I2_1divR_sc / (4*_PI) ;
                       if (p==2) Phi.i +=  w1 * eta * sisi * I3_1divR_sc / (4*_PI) ;

                     } /* for (p=0; p<3; p++) */

                   } /* if (anal_s)  */

                 }  /* if integ_f */
               
               } /* if (Ef) */

               if(Hf)
               { cte = ( unR * rhorho - unRho_s * rrdot(Rv,rho_f) )/4;

                 Ah.r = Gh.r * cte;
                 Ah.i = Gh.i * cte;

                 if (integ_f)
                 {    
		      Ah.r = w0 * Ah.r;
                      Ah.i = w0 * Ah.i;

                      for(j=0; j<3; j++)
                      {
                       cte_intf = ( unR_intf[j] * rhorho_intf[j] - unRho_s * rrdot(Rv_intf[j],rho_f_intf[j]) )/4;
                       Ah.r +=  w1 * Gh_intf[j].r * cte_intf;
                       Ah.i +=  w1 * Gh_intf[j].i * cte_intf;
                      }

                 } /* if integ_f */

               } /* if (Hf) */

          }
          else	/* if (integ_s) */
          { 
	    if(Ef)
            { 
	      sisi = si_s * si_f /k;
              /*Ae = -k * si_s * dot(rho_f, Iv - ovs*I)/4;*/

              for (j=0; j<3; j++)
              { prod2[j].r = Iv[j].r - ovs[j]*I.r;
                prod2[j].i = Iv[j].i - ovs[j]*I.i;
              }
              Ae.r= -k * si_s * rcdot_r(rho_f, prod2)/4;
              Ae.i= -k * si_s * rcdot_i(rho_f, prod2)/4;
              
              /* Phi = I * sisi;*/
              Phi.r = I.r * sisi;
              Phi.i = I.i * sisi;
              
              /* Analytical integration of 1/R ; field point in the centroid */
	      if (anal_s) 
              {  
   	         rhof_A_c = rrdot(rho_f,A_c);

  	         Ae.i  += -0.25 * si_s *  k * eta * rhof_A_c /(4*_PI);  
                 Phi.i +=  eta * sisi * Ic_1divR_sc /(4*_PI) ; 

	      }  /* if (anal_s) */
              

              if (integ_f)
              {   
		  tmp = Ae.r; Ae.r = w0 * tmp;
                  tmp = Ae.i; Ae.i = w0 * tmp;

                  tmp = Phi.r; Phi.r = w0 * tmp;
                  tmp = Phi.i; Phi.i = w0 * tmp;

                  for(p=0; p<3; p++)
                  {
                    for(j=0; j<3; j++)
                    {
                         prod2_intf[p][j].r = Iv_intf[p][j].r - ovs[j]*I_intf[p].r;
                         prod2_intf[p][j].i = Iv_intf[p][j].i - ovs[j]*I_intf[p].i;
                    }
                  }
                  
                  for(j=0; j<3; j++)
                  {
                     Ae.r += - w1 * k * si_s * rcdot_r(rho_f_intf[j], prod2_intf[j])/4;
                     Ae.i += - w1 * k * si_s * rcdot_i(rho_f_intf[j], prod2_intf[j])/4;

                     Phi.r += w1 * I_intf[j].r * sisi;
                     Phi.i += w1 * I_intf[j].i * sisi;
                   }

                  /* Analytical integration of 1/R ; + 3 field points*/
		  if (anal_s) 
		  {	
		    for (p=0; p<3; p++)	
                    { 

 	              rhof_A_intf = rrdot(rho_f_intf[p],A_intf[p]);

  	              Ae.i +=  -w1 * 0.25 * si_s * k * eta * rhof_A_intf /(4*_PI) ; 
                      if (p==0) Phi.i +=  w1 * eta * sisi * I1_1divR_sc /(4*_PI) ;
                      if (p==1) Phi.i +=  w1 * eta * sisi * I2_1divR_sc /(4*_PI) ;
                      if (p==2) Phi.i +=  w1 * eta * sisi * I3_1divR_sc /(4*_PI) ;

                     } /* for (p=0; p<3; p++) */

                   } /* if (anal_s)  */

	        } /* if integ_f */

            } /* if (Ef) */

            if(Hf) 
            { 
	      /*Ah = dot(rho_f, nHs)/4; */ 
              Ah.r = rcdot_r(rho_f, nHs)/4;
              Ah.i = rcdot_i(rho_f, nHs)/4;

              if (integ_f)
              {  Ah.r = w0 * Ah.r;
                 Ah.i = w0 * Ah.i;

                 for(j=0; j<3; j++)
                 {
                   Ah.r +=  w1 * rcdot_r(rho_f_intf[j],nHs_intf[j])/4;
                   Ah.i +=  w1 * rcdot_i(rho_f_intf[j],nHs_intf[j])/4;
                 }

               } /* if (integ_f) */
 
             } /* if (Hf) */

           } /*  else if (integ_s)  */

        } /* if((Ts != Tf) || !singleobject) */
        else      /* Self-triangle integration */
        { if(Ef)  /* Analytical integrals tested at 3 points, only EFIE*/
          { /* Rho of the 3 test points in the field triangle */
            for(j=0; j<3; j++)
            { rho_1[j] = si_f * (rt1[j] - ovf[j]);
              rho_2[j] = si_f * (rt2[j] - ovf[j]);
              rho_3[j] = si_f * (rt3[j] - ovf[j]); 
            }                   
 
            cteG = eta/(4*_PI*ds_s[Ts]);    /* General constant */
            cte = -cteG * k/4;
            /*Ae = cte * iu * (dot(rho_1,Ii1[se]) + dot(rho_2,Ii2[se]) + dot(rho_3,Ii3[se])) /3; */
            a1.r = rcdot_r(rho_1,Ii1[se]); 
            a1.i = rcdot_i(rho_1,Ii1[se]); 
            a2.r = rcdot_r(rho_2,Ii2[se]);
            a2.i = rcdot_i(rho_2,Ii2[se]);
            a3.r = rcdot_r(rho_3,Ii3[se]); 
            a3.i = rcdot_i(rho_3,Ii3[se]); 
            Ae.r = cte*((iu.r*a1.r-iu.i*a1.i)+(iu.r*a2.r-iu.i*a2.i)+(iu.r*a3.r-iu.i*a3.i))/3;
            Ae.i = cte*((iu.i*a1.r+iu.r*a1.i)+(iu.i*a2.r+iu.r*a2.i)+(iu.i*a3.r+iu.r*a3.i))/3;
            cte = cteG * (si_f*si_s) / k;
            /*Phi = cte * iu * (I1+I2+I3) /3; */
            Phi.r = cte * ( iu.r * (I1.r+I2.r+I3.r) - iu.i * (I1.i+I2.i+I3.i) ) /3;
            Phi.i = cte * ( iu.i * (I1.r+I2.r+I3.r) + iu.r * (I1.i+I2.i+I3.i) ) /3;
          } /* if (Ef) */
          if(Hf) /* Compute D/2 */
          { /* rho_3f = ovf-v2f; rho_3s = ovs-v2f; */ 
            for (j=0;j<3;j++)
            { rho_3f[j] = ovf[j]-v2f[j];
              rho_3s[j] = ovs[j]-v2f[j];
              sum2[j] = rho_3s[j] + rho_3f[j]; 
            }
            Ah.r = solid * si_s*si_f * (cteD1 + 6*rrdot(rho_3f,rho_3s) - 2*rrdot(sum2,cteD2))/(48*ds_f[Tf]);
            Ah.i = 0;
            /* 1/48 = 1/12 from I, 1/2 from 1/(2*ds), 1/2 from D/2 */

          } /* if (Hf) */

        }  /* else if((Ts != Tf) || !singleobject) */
 
        cte = ln_s[ge_s] * ln_f[ge_f];
        elem = Elist_f[ge_f]-1 + Nr*(Elist_s[ge_s]-1);

        if(Ef) { 
	        #pragma omp atomic
	        Zr[elem] += (Ae.r + Phi.r) * cte; /* EFIE or CFIE */
	        #pragma omp atomic
                 Zi[elem] += (Ae.i + Phi.i) * cte; /* EFIE or CFIE */
               }

        if(Hf) {         
		 #pragma omp atomic
		 Zr[elem] += (field==2)? Ah.r*cte : Ah.r*cte*eta; /* MFIE or CFIE */
		 #pragma omp atomic
                 Zi[elem] += (field==2)? Ah.i*cte : Ah.i*cte*eta; /* MFIE or CFIE */
               }

      } /* for fe */
    } /* for se */
  } /* for Ts */
  } /* for Tf */
}/*omp parallel zone*/

}
/*}  extern "C" */


void autoz(dcomplex *pI, dcomplex Ii[3][3], int Ts, triad *topol, triad *vertex, triad *trian, double rt[3], triad *un, double k)
{ double nrc1, nrc2, nrc3, ph12, ph23, ph31, p12, p23, p31, a12, a23, a31;
  triad r1, r2, r3, rc1, rc2, rc3;
  triad dif1, dif2, dif3;
  int j;
  double ph, kr;
  int i;
  dcomplex I_1; 
  dcomplex I_2;   
  dcomplex z, z_1, z_2;
  double k2;
  dcomplex tmp;
  double rc1un[3];
  triad si_s;
  dcomplex tmpv[3]; 
  double x[Gp], w[Gp];

  for (j=0; j<3; j++)
  /*double r1[3] = vertex[ (int)topol[Ts][0]-1 ]; */
  { r1[j] = vertex[(int)topol[Ts][0]-1][j]; /*vertex 1*/
    /*vector r2 = vertex[ (int)topol[Ts][1]-1 ]; */ 
    r2[j] = vertex[ (int)topol[Ts][1]-1][j]; /* Vertex */
    /*vector r3 = vertex[ (int)topol[Ts][2]-1 ]; */ 
    r3[j] = vertex[ (int)topol[Ts][2]-1][j];  /* Vertex 3 */
    /*vector rc1 = r1 - rt; */  /* Vector from center to vertex 1 */
    rc1[j] = r1[j] - rt[j];
    /*vector rc2 = r2 - rt;   */ /* Vector from center to vertex 2 */
    rc2[j] = r2[j] - rt[j];
    /* zvector rc3 = r3 - rt;   */  /*Vector from center to vertex 3 */
    rc3[j] = r3[j] - rt[j];
  }
  
  nrc1 = norm(rc1); nrc2 = norm(rc2); nrc3 = norm(rc3);
  
  ph12 = acos(rrdot(rc1,rc2)/(nrc1*nrc2));
  ph23 = acos(rrdot(rc2,rc3)/(nrc2*nrc3));
  ph31 = 2*_PI - ph12 - ph23;

  /* Parameters p and a for polar coordinates: r = p / cos(ph - a) */
  for (j=0; j<3; j++)
  { dif1[j] = r1[j] - r2[j];
    dif2[j] = r2[j] - r3[j];
    dif3[j] = r3[j] - r1[j];
  } 
  /* p12 = nrc1*nrc2*sin(ph12) / norm(r1-r2); */
  p12 = nrc1*nrc2*sin(ph12) / norm(dif1); 
  /* p23 = nrc2*nrc3*sin(ph23) / norm(r2-r3); */
  p23 = nrc2*nrc3*sin(ph23) / norm(dif2);
  /* p31 = nrc3*nrc1*sin(ph31) / norm(r3-r1); */
  p31 = nrc3*nrc1*sin(ph31) / norm(dif3);

  a12 = acos(p12 / nrc1);
  a23 = acos(p23 / nrc2);
  a31 = acos(p31 / nrc3);

  /* Gauss integration points in [0 1] and weigths */
  /* double x[Gp], w[Gp]; */ /* 3 Gauss points */
  x[0] = (1-sqrt(0.6))/2;
  x[1] = 0.5;
  x[2] = (1+sqrt(0.6))/2;
  w[0] = 5.0/18;
  w[1] = 8.0/18;
  w[2] = w[0];
  /* double ph, kr; int i; complex z; */

  pI->r = 0;   /* Scalar potential */
  pI->i = 0;  

  /* for(i=0, z=0; i<Gp; i++) */
  z.r=0; z.i=0;
  for(i=0 ; i<Gp; i++) 
  { ph = x[i] * ph12;       kr = k*p12/cos(ph-a12);
    /*z += w[i]*complex(sin(kr),cos(kr)-1); */
    z.r += w[i]*sin(kr);
    z.i += w[i]*(cos(kr)-1);
  } /* *pI += z * ph12/k; */
  pI->r += z.r * ph12/k; 
  pI->i += z.i * ph12/k; 


  /* for(i=0, z=0; i<Gp; i++) */
  z.r=0; z.i=0;
  for(i=0; i<Gp; i++) 
  { ph = x[i] * ph23;       kr = k*p23/cos(ph-a23);
    /*z += w[i] * complex(sin(kr),cos(kr)-1); */
    z.r += w[i] * sin(kr); 
    z.i += w[i] * (cos(kr) -1);
  } /* *pI += z * ph23/k; */
  pI->r += z.r * ph23/k; 
  pI->i += z.i * ph23/k; 

  /* for(i=0, z=0; i<Gp; i++) */
  z.r=0; z.i=0;
  for(i=0; i<Gp; i++) 
  { ph = x[i] * ph31;       kr = k*p31/cos(ph-a31);
  /* z += w[i] * complex(sin(kr),cos(kr)-1); */
    z.r += w[i] * sin(kr);
    z.i += w[i] * (cos(kr)-1);
  } /* *pI += z * ph31/k; */
  pI->r += z.r * ph31/k; 
  pI->i += z.i * ph31/k; 
        
  /* Vector potential */
  k2 = k*k;

  /*complex I_1=0; */ /* Integral of ru, 1st componet */ 
  /* dcomplex I_1 */
  I_1.r=0;
  I_1.i=0;
  /*complex I_2=0;   */ /* Integral of ru, 2nd component */
  /* dcomplex I_2 */
  I_2.r=0;
  I_2.i=0;
  /*dcomplex z_1, z_2; */

  /* for(i=0, z_1=0, z_2=0; i<Gp; i++) */
  z_1.r=0; z_1.i=0; z_2.r=0; z_2.i=0;
  for(i=0; i<Gp; i++)
  { ph = x[i] * ph12;       kr = k*p12/cos(ph-a12);
    /*tmp = w[i] * ( complex(cos(kr),-sin(kr)) * complex(1,kr) -1); */
    tmp.i = w[i] * ( cos(kr)*kr - sin(kr) );
    tmp.r = w[i] * ( cos(kr) + sin(kr)*kr -1);
    /*z_1 += tmp*cos(ph); z_2 += tmp*sin(ph); */
    z_1.r += tmp.r*cos(ph); z_2.r += tmp.r*sin(ph);
    z_1.i += tmp.i*cos(ph); z_2.i += tmp.i*sin(ph);
  }
  /*I_1 += z_1 * ph12/k2; I_2 += z_2 * ph12/k2; */
  I_1.r += z_1.r * ph12/k2; I_2.r += z_2.r * ph12/k2;
  I_1.i += z_1.i * ph12/k2; I_2.i += z_2.i * ph12/k2;

  /*for(i=0, z_1=0,z_2=0; i<Gp; i++) */
  z_1.r=0; z_1.i=0; z_2.r=0; z_2.i=0;
  for(i=0; i<Gp; i++) 
  { ph = x[i] * ph23;       kr = k*p23/cos(ph-a23);
    /*tmp = w[i] * ( complex(cos(kr),-sin(kr)) * complex(1,kr) -1); */
    tmp.i = w[i] * ( cos(kr)*kr - sin(kr) );
    tmp.r = w[i] * ( cos(kr) + sin(kr)*kr -1);
    /*z_1 += tmp*cos(ph+ph12); z_2 += tmp*sin(ph+ph12); */
    z_1.r += tmp.r*cos(ph+ph12); z_2.r += tmp.r*sin(ph+ph12);
    z_1.i += tmp.i*cos(ph+ph12); z_2.i += tmp.i*sin(ph+ph12);
  }
  /*I_1 += z_1 * ph23/k2; I_2 += z_2 * ph23/k2; */
  I_1.r += z_1.r * ph23/k2; I_2.r += z_2.r * ph23/k2;
  I_1.i += z_1.i * ph23/k2; I_2.i += z_2.i * ph23/k2;

  /*for(i=0, z_1=0, z_2=0; i<Gp; i++)*/
  z_1.r=0; z_1.i=0; z_2.r=0; z_2.i=0;
  for(i=0; i<Gp; i++)
  { ph = x[i] * ph31;       kr = k*p31/cos(ph-a31);
    /*tmp = w[i] * ( complex(cos(kr),-sin(kr)) * complex(1,kr) -1); */
    tmp.i = w[i] * ( cos(kr)*kr - sin(kr) );
    tmp.r = w[i] * ( cos(kr) + sin(kr)*kr -1);
    /*z_1 += tmp*cos(ph+ph12+ph23); z_2 += tmp*sin(ph+ph12+ph23); */
    z_1.r += tmp.r*cos(ph+ph12+ph23);
    z_1.i += tmp.i*cos(ph+ph12+ph23);
    z_2.r += tmp.r*sin(ph+ph12+ph23); 
    z_2.i += tmp.i*sin(ph+ph12+ph23); 
  }
  /*I_1 += z_1 * ph31/k2; I_2 += z_2 * ph31/k2; */
  I_1.r += z_1.r * ph31/k2; I_2.r += z_2.r * ph31/k2;
  I_1.i += z_1.i * ph31/k2; I_2.i += z_2.i * ph31/k2;

  /* vector rc1un = cross(rc1,un[Ts]); */ 
  /* triad rc1un;*/ /* Reference */
  rc1un[0] = rc1[1]*un[Ts][2]-rc1[2]*un[Ts][1];
  rc1un[1] = rc1[2]*un[Ts][0]-rc1[0]*un[Ts][2];
  rc1un[2] = rc1[0]*un[Ts][1]-rc1[1]*un[Ts][0];

  /* Sign of rho vectors for source edges */ 
  /*triad si_s;*/
  for(i=0; i<3; i++) si_s[i] = (trian[Ts][i]>0)? 1: -1;
  
  /*cvector tmpv = (rc1*I_1 + rc1un*I_2) / nrc1; */
  /*dcomplex tmpv[3]; */
  for(j=0; j<3; j++)
  { /* tmpv[j] = (rc1[j]*I_1 + rc1un[j]*I_2) / nrc1; */
    tmpv[j].r = (rc1[j]*I_1.r + rc1un[j]*I_2.r) / nrc1;
    tmpv[j].i = (rc1[j]*I_1.i + rc1un[j]*I_2.i) / nrc1;
  }
  
  /* Ii[0] = si_s[0] * (-rc1*(*pI) + tmpv); */
  for(j=0; j<3; j++)
  { /*  Ii[0][j] = si_s[0] * (-rc1[j]*(*pI) + tmpv[j]); */
    Ii[0][j].r = si_s[0] * (-rc1[j]*(pI->r) + tmpv[j].r);
    Ii[0][j].i = si_s[0] * (-rc1[j]*(pI->i) + tmpv[j].i);

    /* Ii[1] = si_s[1] * (-rc2*(*pI) + tmpv); */
    /*Ii[1][j] = si_s[1] * (-rc2[j]*(*pI) + tmpv[j]);*/
    Ii[1][j].r = si_s[1] * (-rc2[j]*(pI->r) + tmpv[j].r); 
    Ii[1][j].i = si_s[1] * (-rc2[j]*(pI->i) + tmpv[j].i); 
    /* Ii[2] = si_s[2] * (-rc3*(*pI) + tmpv); */
    /* Ii[2][j] = si_s[2] * (-rc3[j]*(*pI) + tmpv[j]); */
    Ii[2][j].r = si_s[2] * (-rc3[j]*(pI->r) + tmpv[j].r);
    Ii[2][j].i = si_s[2] * (-rc3[j]*(pI->i) + tmpv[j].i);
  }
}

void anal_1divR(double *I_sc, double I_vc[3], int Ts, triad *topol, triad *vertex, triad *trian, double rt[3], triad *un, double *ds, double k, double signe)
{ int i, j, pp, pm, p;
  double rt_prj_un;
  double norm_d, Rs_pp, Rs_pm, l_pp, l_pm, Po, Po_sc, Ro;
  double Rsl_pp, Rsl_pm, Ro_dpp, Ro_dpm, Po_prj_u, norm_vec_l;
  double Rs[3], d[3], rv_prj_un[3], rho[3], vec_l[3], vec_u[3], vec_Po[3];
  double dif_rho_pp[3], dif_rho_pm[3];
  double dif[3][3], rv[3][3], rho_v[3][3];
  double un_anal[3];

  for (j=0; j<3; j++)
  {
   un_anal[j] = signe * un[Ts][j]; 
  }
  
  for (i=0; i<3; i++)
  {
    for (j=0; j<3; j++)
    { 
      rv[i][j] = vertex[ (int)topol[Ts][i] - 1][j]; 
      dif[i][j] = rt[j] - rv[i][j];
    }
  }

  for (i=0; i<3; i++)
  {
  Rs[i] = norm(dif[i]);
  }

  rt_prj_un = rrdot(rt,un_anal);
  
  for (i=0; i<3; i++)
  {
   rv_prj_un[i] = rrdot(rv[i],un_anal);
  }
  
  for (i=0; i<3; i++)
  {
  d[i] = ( rt_prj_un - rv_prj_un[0] ) * un_anal[i];
  }
  
  for (i=0; i<3; i++)
  {
    for (j=0; j<3; j++)
    { 
     rho_v[i][j] = rv[i][j] - rv_prj_un[i] * un_anal[j];
     }
  }

  for (i=0; i<3; i++)
  {
    rho[i] = (rt[i] - d[i]) - (rv[0][i] - rho_v[0][i]);
  }

  norm_d = norm(d);

  *I_sc = 0;

  for (j=0; j<3; j++)
  {
    I_vc[j] = 0;
  }   

  for (pp=0; pp<3; pp++) 
  {

    if (pp==0) pm=1;	
    if (pp==1) pm=2;	
    if (pp==2) pm=0; 	

    Rs_pp  =  Rs[pp];
    Rs_pm  =  Rs[pm];

    for (j=0; j<3; j++)
    {
      vec_l[j] = rho_v[pp][j] - rho_v[pm][j]; 
    }
    norm_vec_l = norm(vec_l);

    for (j=0; j<3; j++)
    {
      vec_l[j] = vec_l[j]/norm_vec_l;
    }     
    
    cross_rr(vec_u,vec_l,un_anal);

    for (j=0; j<3; j++)
    {
      dif_rho_pp[j] = rho_v[pp][j] - rho[j];
      dif_rho_pm[j] = rho_v[pm][j] - rho[j]; 
    }

    l_pp = rrdot(dif_rho_pp,vec_l);
    l_pm = rrdot(dif_rho_pm,vec_l);

    Po_sc = rrdot(dif_rho_pp,vec_u);

    Po = sqrt( pow(Po_sc,2) );  
    Ro = sqrt( pow(Po,2) +  pow(norm_d,2) );

    for(j=0; j<3; j++)
    {
    vec_Po[j] = 0;
    }

    if (Po > 1e-10)
    {
      for(j=0; j<3; j++)
      {  
       vec_Po[j] = ( dif_rho_pp[j] - l_pp*vec_l[j] )/Po;
      }     
    }

    Rsl_pp = Rs_pp + l_pp;
    Rsl_pm = Rs_pm + l_pm;

    Ro_dpp = Ro*Ro + norm_d*Rs_pp;
    Ro_dpm = Ro*Ro + norm_d*Rs_pm;

    Po_prj_u = rrdot(vec_Po,vec_u);

    /* I_sc update */
    if ( (Rsl_pp>1e-10) && (Rsl_pm>1e-10) && (Ro_dpp>1e-10) && (Ro_dpm>1e-10) ) 
    {
    *I_sc = *I_sc +  Po_prj_u * ( Po*log(Rsl_pp/Rsl_pm) - norm_d*( atan((Po*l_pp)/Ro_dpp) - atan((Po*l_pm)/Ro_dpm) ) )/ds[Ts];
    }

    /* I_vc update */
    if ( (Rsl_pp >1e-10) && (Rsl_pm>1e-10) ) 
    {
      for(j=0; j<3; j++)
      {
	I_vc[j]  =  I_vc[j]  +  0.5*vec_u[j]*( (Ro*Ro)*log(Rsl_pp/Rsl_pm) + l_pp*Rs_pp - l_pm*Rs_pm )/ds[Ts]; 
      }
    } 

  } /* for (pp=0; pp<3; pp++) */

}

