// =============================================================================
//                                      GHW1 
// =============================================================================
// 
// 2D delta-f gyrofluid (modified) Hasegawa-Wakatani Code with periodic b.c.
// (by Alexander Kendl, Innsbruck; ghw1 version V1.25 of 10.01.2025)
//
// ( changes to V1.24: - speedup of FFT by using r2c/c2r instead of complex FFT
//                     - put initial conditions into functions for clarity )
//
// Physics: 
// Simulation of 2D drift wave turbulence and vortex dynamics in magnetised plasma.
// In cold ion (taui=0) limit: comparable to results by Numata in arXiv:0708.4317.
//
// Examples for usage (V1.0) : arXiv:1708.06213
// Normalisation and numerical methods: in similar full-f version see arXiv:1502.05494.
// More background on gyrofluid theory and models see references therein.
// Parameters are set in 'ghw1.inp' input file. Explanations given within the code below.
//
// Dependencies:
// Needs some libraries (e.g. fftw3, omp) available in every Linux distribution.
//
// For command line execution under Linux first apply these settings:
//   ulimit -s unlimited
//   export OMP_STACKSIZE=8G
// Compilation with attached makefile:
//   make ghw1
// Execution with run time control:
//   time ./ghw1
//
// Parallelisation with OpenMP is not extremely efficient, but good speed-up up to np=8.
// Run time depends on problem size: a few minutes for 64^2, many hours for 1024^2.
//
// In this code version most subroutines are inlined instead of calling functions:
// allthough not elegant in coding, this considerably speeds up the code execution.
//
// Time step and hyperviscosity have to be refined with the grid resolution
// to prevent crashing through numerical instability. (Set "dt" and "diff" by trial and error.)
// Time step (and typical time scales) also sensitive to choice of "chat" and "delta".
// 
// Running in a RAM disc folder (here: at location '/ramdisc/') saves some output time:
// 'mkdir /ramdisc', 'sudo mount -t tmpfs none /ramdisc -o size=2048m'.
// 
// Visualisation of 1D graphs with any plotting programme (e.g. xmgrace),
// visualisation of 2D contour plots with attached gnuplot script file "gpl-ctr.plt".
// Watch the simulation live by uncommenting "reread". In gnuplot: "load 'gpl-ctr.plt'".
//
// Makefile example:
// exec = ghw1
// objs = ghw1.o 
// LFLAGS = -mcmodel=large -march=native -lm -O2 -fopenmp -lfftw3_omp -lfftw3 
// CFLAGS = -I.
// LIBS = /usr/lib64 
// all: $(exec)
// $(exec): $(objs)
//	g++  $(exec).C -o $@  $(LFLAGS) -L$(LIBS)
//	rm *.o
// ghw.o: ghw1.C 
// 	g++  -c $< -o $@
//
// ==============================================================================

// Include libraries: available in standard Linux distributions (tested: OpenSuse)
#ifdef _OPENMP 
#include <omp.h>
#endif

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fftw3.h>

// Define static constants:
const static double TwoPi = 2.*M_PI, r12 = 1./12.;

// nxm, nym: maximum nx and ny sizes, to initialise arrays:
const static int nxm=512+2, nym=512+2;  // nx + 2 = 2^n + 2  ideal for fft

 // number of time trace recording points (e.g. 8) per x and per y axis each:
const static int kt_maxsqrt = 8; 

static unsigned int fflag = FFTW_ESTIMATE;

typedef double AXY[nxm][nym];

static int nx, ny, nxh, nyh, nx1, ny1, ict, jct;
static int itstp, itmax, npar, ntrace, jnull_old;
static int ipade, ihype;
static int maxtrace;

static double p_amp = 0., p_old = 0., t_amp = 0., t_old = 0.;
static double t_1, t_2;
static int jsgn, jsgn_old;

static double pkxavg[nxm/2+1], pkyavg[nym/2+1];
static double pkyavgn[nxm/2+1], pkyavgw[nym/2+1], pkyavge[nym/2+1];

static double hy, hy2, hysq, xyz, dr, amp, wsrc, sigma, xlwl;
static double diff, n00, ly, dt, dtt, ddtt, incon;
static double chat, diss, neavg, niavg, peavg, vorfree;
static double aspect, eno, enwo, delta, delinv, g_n;
static double pdist, vorset, clin, vde_old, vde_avg;
static double aae, aai,  mue, mui, taue, taui, zze, zzi;

// electrostatic potential pe, gyroaveraged pi, ww:
static AXY pe, pi, ww;

// electron and ion gyrocenter densities ne and ni, and stored values:
static AXY nn, ne, ni, ne0, ne1, ne2, ni0, ni1, ni2;

// gyroscreened ion density:
static AXY gni;

// terms defining r.h.s. of ne and ni continuity equations
static AXY fne0, fne1, fne2, fni0, fni1, fni2;

// arrays to store Laplacians (e.g. for viscosity)
static AXY vis_e, vis_i, hyve, hyvi;

// FFT kernels
static AXY cvort, cpoti, clap, chyv, cvfis;

static char *wisdom_sf;

static bool printtraces, b_mhw, b_4th;

// Functions first defined here, see further below main part:

void init_parameters( void );
void init_add_blob( AXY ne, AXY ni );
void init_add_turb( AXY ne, AXY ni );
void init_add_dual( AXY ne, AXY ni );
void init_add_flow( AXY ne, AXY ni );
void init_add_mode( AXY ne, AXY ni );

void arakawa( AXY uuu, AXY vvv, AXY www );
void arakaw4( AXY uuu, AXY vvv, AXY www );
void poisson( AXY pe, AXY cp, AXY pz );
void laplace( AXY pp, AXY ww );
void diagnose(double ttt, double t00, int it, int itmax,
	      AXY nn, AXY ne, AXY ni, AXY pp, AXY gni );

static inline double DFDY( AXY arrin, int i, int j );
static inline double DFDX( AXY arrin, int i, int j );
static inline void nextp( int i, int i0, int i1, int& im, int& ip, int& im2, int& ip2 );
static inline void f_copy2darray( AXY arrinp, AXY arrout, int iend, int jend);

static double gamma0( double bkk );
inline double timer_start();
inline double timer_stop(double t_1);


// --------------------------- Main programme part ---------------------------------

int main(void)
{ 
  int i, j, k, l, m, n, ik, jk, it, is, itn, itend;
  int iout[kt_maxsqrt*kt_maxsqrt], jout[kt_maxsqrt*kt_maxsqrt], k_max;
  int im, ip, jm, jp, im1, ip1, im2, ip2, jm2, jp2, icc, ixloc, itrace, jtrace;

  double c0 = 18./11., c1 = 9./11., c2 = 2./11., cf = 6./11.;
  double ttt, t00, phase, zuf, nue, nui, nuzf;
  double neprof, pbar, nbar, dve, gre, gri;
 
  // variables to measure time consumption of programme parts
  double td_init=0., td_update=0., td_pois1=0., td_pois2=0., td_out=0.;
  double td_bnd=0.,td_mem=0., td_end=0., td_pol=0., td_tot;

  AXY ane, ani, dummy;
  double bndys[nxm];
  double ne_za[nxm], ni_za[nxm];
  char s[80], str[80]; 
  char const *fn, *fw;
  FILE *f, *g, *h, *g1, *g2, *g3, *g4, *g5, *ft;
  
  t_1 = timer_start(); // start measure td_init (analyse relative performance of code parts)
  
  // read input file and initialise parameters:
  init_parameters(); 
  
  //  initialise parallelized FFTW:
#ifdef _OPENMP 
  fftw_plan_with_nthreads(npar); 
#endif
  int nyp = ny/2+1;
  fftw_complex wwk[nx][nyp];
  double wwx[nx][ny];
  
  fftw_plan hinp, herp;
  hinp=fftw_plan_dft_r2c_2d(nx,ny,*wwx,*wwk,fflag);
  herp=fftw_plan_dft_c2r_2d(nx,ny,*wwk,*wwx,fflag);
  wisdom_sf = fftw_export_wisdom_to_string(); 
  
  // optional: define Dirichlet boundary condition (if chosen):
  for (i=0; i<=nx1; i++) {
    int nx_i = double(nx1-i); double dx0 = 2.*2.;
    bndys[i] = 1. - exp(-nx_i*nx_i/dx0) - exp(-double((i-0)*(i-0))/dx0);
  }
  
  // set initial density perturbation (blob, vortices, turbulence, flow, ...)
  for (i=0; i<=nx1; ++i) for (j=0; j<=ny1; j++) { ne[i][j] = 0.; ni[i][j] = 0.;} 

  if (incon==0.) init_add_blob(ne,ni); // Gaussian density blob initialisation
  if (incon==1.) init_add_turb(ne,ni); // pseudo-turbulent bath initialisation
  // (incon==2:  restart -> see below...)  
  if (incon==3.) init_add_dual(ne,ni); // dual vortex merger initialisation
  if (incon==4.) init_add_flow(ne,ni); // shear flow initialisation
  if (incon==5.) init_add_mode(ne,ni); // "linear" single ky mode initialisation

  // compute initial potential
  poisson(ni,cpoti,gni);
  for (i=0; i<=nx1; ++i) for (j=0; j<=ny1; j++) ww[i][j] = - aae*ne[i][j] - aai*gni[i][j];
  poisson(ww,cvort,pe);

 
  // set artificial "previous" time values (for multi-step solver)
  #pragma omp parallel for private(i,j) 
  for (i=0; i<=nx1; ++i) for (j=0; j<=ny1; ++j) {
      ne0[i][j] = ne[i][j]; ne1[i][j] = ne[i][j]; ne2[i][j] = ne[i][j];
      ni0[i][j] = ni[i][j]; ni1[i][j] = ni[i][j]; ni2[i][j] = ni[i][j];
      fne1[i][j] = 0.; fne2[i][j] = 0.; fni1[i][j] = 0.; fni2[i][j] = 0.;
    }

  // initialise energy:
  eno = 0.;
  #pragma omp parallel for private(i,j) reduction (+:eno)
  for (i=0; i<=nx1; i++) for (j=0; j<=ny1; j++) eno += exp(2.*ne[i][j]); 
  eno *= .5*xyz;
  
  if (incon==2.) { // restart from saved end output of previous run
    printf("| continued data set... \n");
    int rs_check = 1;   
    g = fopen( "restart.dat", "r" );
    if (g == NULL) { printf("\n File 'restart.dat' not found. Exiting ... \n"); exit(1); }
  
    fscanf(g,"%s ",str);  t00  = atof(str); 
    while (!feof(g)) { 
      rs_check *= fscanf(g,"%s ",str);  i = atoi(str); 
      rs_check *= fscanf(g,"%s ",str);  j = atoi(str); 
      rs_check *= fscanf(g,"%s ",str);  eno = atof(str);
      rs_check *= fscanf(g,"%s ",str);  pe[i][j]  = atof(str); 
      rs_check *= fscanf(g,"%s ",str);  ne[i][j]  = atof(str); 
      rs_check *= fscanf(g,"%s ",str);  ne1[i][j] = atof(str); 
      rs_check *= fscanf(g,"%s ",str);  ne2[i][j] = atof(str); 
      rs_check *= fscanf(g,"%s ",str);  ni[i][j]  = atof(str); 
      rs_check *= fscanf(g,"%s ",str);  ni1[i][j] = atof(str); 
      rs_check *= fscanf(g,"%s ",str);  ni2[i][j] = atof(str); 	   
      rs_check *= fscanf(g,"%s ",str);  fne1[i][j] = atof(str); 
      rs_check *= fscanf(g,"%s ",str);  fne2[i][j] = atof(str); 
      rs_check *= fscanf(g,"%s ",str);  fni1[i][j] = atof(str); 
      rs_check *= fscanf(g,"%s ",str);  fni2[i][j] = atof(str); 
    }
    fclose( g ); 
    poisson(pe,cpoti,pi);       
    for (i=0; i<=nx1; i++) for (j=0; j<=ny1; j++) { ne0[i][j]=ne[i][j]; ni0[i][j]=ni[i][j]; }
  }

  // set the traces output locations
  double trace_n[kt_maxsqrt*kt_maxsqrt][maxtrace], trace_p[kt_maxsqrt*kt_maxsqrt][maxtrace];
  int kk = 0;
  k_max =  kt_maxsqrt*kt_maxsqrt;
  for (int ki=0; ki<kt_maxsqrt; ++ki)
    for (int kj=0; kj<kt_maxsqrt; ++kj)
      { iout[kk] = nx1*ki/kt_maxsqrt; jout[kk] = ny1*kj/kt_maxsqrt; kk++;  }


  // initialise time stepping:
  printf("| start time stepping ...\n");
  itn = 0;
  if (itstp!=0) itend = itmax / itstp;
  if (itstp==0) itend==0;
  dtt = dt*cf;
  ddtt = double(itstp)*dt;
  itrace = 0; jtrace = 0;
  if (itmax == 0) return(0);
  
  // control output of initial profiles
  diagnose(0.,0.,it,itmax,nn,ne,ni,pe,gni); 

  td_init = timer_stop(t_1); // stop measure td_init
 
  // Time step  --------------------------------------------------------------
  // 3rd order Karniadakis time scheme; Arakawa brackets for convective term

  for (it=itn/itstp; it<itend; ++it) { // total time loop (output step)
    for (is=0; is<itstp; ++is) {    // inner time loop (without outputs)
      ttt = t00 + (it+0)*ddtt + (is+1)*dt;
	  
      // time step update of gyrofluid densities ne and ni:
      t_1 = timer_start(); // start measure td_update
 
      // ExB convection term: compute Poisson brackets [n,p] by Arakawa scheme
      if (b_4th) { arakaw4(ne,pe,ane); arakaw4(ni,pi,ani); }
      else { arakawa(ne,pe,ane); arakawa(ni,pi,ani); }; // 2nd order more stable for larger dt

      // choose how to calculate the (hyper)viscosity Laplacians
      if (ihype==-1) { // evaluate del^p by (p/2) laplace on Laplace operators
	laplace(ne,hyve); 
	laplace(ni,hyvi);
      }
      else if (ihype==0) { // evaluate del^p by (p/2) laplace on Laplace operators
	laplace(ne,vis_e); laplace(vis_e,hyve); 
	laplace(ni,vis_i); laplace(vis_i,hyvi);
      }
      else if (ihype>0) { // evaluate FFT of del^p as k^p: recommended (p==2)

#pragma omp parallel for private(j)
	for (i=0; i<=nx1; ++i) for (j=0; j<=ny1; ++j) { wwx[i][j] = ne[i][j]; }	  
	fftw_execute(hinp);
#pragma omp parallel for private(j,k)
	for (i=0; i<=nx1; ++i) for (j=0; j<=nyp-1; ++j) {
	    for (k=0; k<=1; k++) wwk[i][j][k] *= chyv[i][j];
	  }
	fftw_execute(herp);
#pragma omp parallel for private(j)
	for (i=0; i<=nx1; ++i) for (j=0; j<=ny1; ++j) { hyve[i][j] = wwx[i][j]; }
	
#pragma omp parallel for private(j)
	for (i=0; i<=nx1; ++i) for (j=0; j<=ny1; ++j) { wwx[i][j] = ni[i][j]; }
	fftw_execute(hinp);
#pragma omp parallel for private(j,k)
	for (i=0; i<=nx1; ++i) for (j=0; j<=nyp-1; ++j) {
	    for (k=0; k<=1; k++) wwk[i][j][k] *= chyv[i][j];
	  }
	fftw_execute(herp);
#pragma omp parallel for private(j)
	for (i=0; i<=nx1; ++i) for (j=0; j<=ny1; ++j) { hyvi[i][j] = wwx[i][j]; }

      } // end hyperviscosity calculations
       
      

      // update ne and ni:
#pragma omp parallel for private(j,neavg,peavg,nue,nui,dve)
      for (i=0; i<=nx1; i++) {
	
	// MHW ("modified Hasegawa-Wakatani") average terms:
	neavg = 0.; peavg = 0.;
	if (b_mhw) {
	  for (j=0; j<=ny1; j++) { neavg += ne[i][j]; peavg += pe[i][j]; }
	  neavg/=ny; peavg/=ny;
	}

	for (j=0; j<=ny1; j++) {

	  // MHW drive: nonadiabatic (dissipative) coupling term
	  dve = chat*( (ne[i][j]-neavg)-(pe[i][j]-peavg) );

	  // densities (continuity equations) r.h.s.
	  fne0[i][j]= clin*delinv*ane[i][j] - g_n*DFDY(pe,i,j) - dve;
	  fni0[i][j]= clin*delinv*ani[i][j] - g_n*DFDY(pi,i,j);
	  
	  // hyperviscosities
	  nue = - diff*hyve[i][j]; 
	  nui = - diff*hyvi[i][j];
	  
	  // time step update of gyro densities
	  ne[i][j] = c0*ne0[i][j] - c1*ne1[i][j] + c2*ne2[i][j] 
	    + dtt*(3.*fne0[i][j] - 3.*fne1[i][j] + fne2[i][j] + nue);
	  ni[i][j] = c0*ni0[i][j] - c1*ni1[i][j] + c2*ni2[i][j] 
	    + dtt*(3.*fni0[i][j] - 3.*fni1[i][j] + fni2[i][j] + nui);
	}
      }

      td_update += timer_stop(t_1); // stop measure td_update

      // test: damping of zonal densities:
      /*
      nuzf = 0.0; // 0.02
      nuzf *= dtt/double(ny);
#pragma omp parallel for private(j, ne_za, ni_za) 
      for (i=0; i<=nx1; ++i) {
	ne_za[i] = 0.; ni_za[i] = 0.;
	// zonal averages of ne and ni:
	for (j=0; j<=ny1; ++j) { ne_za[i] += ne[i][j]; ni_za[i] += ni[i][j]; };
	// reduce ne and ni by a fraction of the zonal densities:
	for (j=0; j<=ny1; ++j) { ne[i][j] -= nuzf*ne_za[i]; ni[i][j] -= nuzf*ni_za[i]; };
     }
      */


       // gyro-averaging of ni to gni: -------------------------------
      t_1 = timer_start(); // start measure td_pois1
      if (taui>0.) { 
#pragma omp parallel for private(j) 
	for (i=0; i<=nx1; ++i) for (j=0; j<=ny1; ++j) wwx[i][j] = ni[i][j]; 
	fftw_execute(hinp); // threaded FFTW call
#pragma omp parallel for private(j,k) 
	for (i=0; i<=nx1; ++i) for (j=0; j<=nyp-1; ++j) {
	    for (k=0; k<=1; k++) wwk[i][j][k] *= cpoti[i][j];
	  }
	fftw_execute(herp); // threaded FFTW call
#pragma omp parallel for private(j) 
	for (i=0; i<=nx1; ++i) for (j=0; j<=ny1; ++j) gni[i][j] = wwx[i][j];
      }
      else { f_copy2darray(ni,gni,nx1,ny1); } // shortcut if taui == 0 
      td_pois1 += timer_stop(t_1); // stop measure td_pois1


      // boundary condition (channel geometry consistent with TIFF code local b.c.):
      t_1 = timer_start(); // start measure td_bnd
      if (vorfree==2.) {
#pragma omp parallel for private(j)
	for (i=0; i<=nx1; i++) for (j=0; j<=ny1; j++) {
	    ne[i][j]*=bndys[i]; ni[i][j]*=bndys[i]; gni[i][j]*=bndys[i];
	  }
      }
      td_bnd += timer_stop(t_1); // stop measure td_bnd


      // solve delta-f gyrofluid polarisation equation: ----------------------------
      // computes electric potential phi from vorticity by FFT in k-space
      t_1 = timer_start(); // start measure td_pol

#pragma omp parallel for private(j) 
      for (i=0; i<=nx1; ++i) for (j=0; j<=ny1; ++j) {
	  wwx[i][j]   = - zze*ne[i][j] - zzi*gni[i][j]; // (equal to vorticity for taui=0)	  
	}
      fftw_execute(hinp);
#pragma omp parallel for private(j,k) 
      for (i=0; i<=nx1; ++i) for (j=0; j<=nyp-1; ++j)  {
	  for (k=0; k<=1; k++) wwk[i][j][k] *= cvort[i][j];
	}
      fftw_execute(herp);
#pragma omp parallel for private(j) 
      for (i=0; i<=nx1; ++i) for (j=0; j<=ny1; ++j) pe[i][j] = wwx[i][j];

      td_pol += timer_stop(t_1); // stop measure td_pol


      // gyro-screened potential pi: ---------------------------------------------
      t_1 = timer_start(); // start measure td_pois2
      if (taui>0.) {	
#pragma omp parallel for private(j) 
	for (i=0; i<=nx1; ++i) for (j=0; j<=ny1; ++j) wwx[i][j] = pe[i][j]; 
	fftw_execute(hinp);
#pragma omp parallel for private(j,k) 
	for (i=0; i<=nx1; ++i) for (j=0; j<=nyp-1; ++j) {
	    for (k=0; k<=1; k++) wwk[i][j][k] *= cpoti[i][j];
	  }
	fftw_execute(herp);
#pragma omp parallel for private(j) 
	for (i=0; i<=nx1; ++i) for (j=0; j<=ny1; ++j) pi[i][j] = wwx[i][j]; 
      }
      else { f_copy2darray(pe,pi,nx1,ny1); } // shortcut for taui == 0
      td_pois2 += timer_stop(t_1); // stop measure td_pois2

            
      // multi-step memory: remember values two steps backwards -------------------
      t_1 = timer_start(); // start measure td_mem
#pragma omp parallel for private(j) collapse(2)
      for (i=0; i<=nx1; ++i) for (j=0; j<=ny1; ++j) { 
	  ne2[i][j] = ne1[i][j]; ne1[i][j] = ne0[i][j]; ne0[i][j] = ne[i][j]; 
	  ni2[i][j] = ni1[i][j]; ni1[i][j] = ni0[i][j]; ni0[i][j] = ni[i][j];
	  fne2[i][j] = fne1[i][j]; fne1[i][j] = fne0[i][j];
	  fni2[i][j] = fni1[i][j]; fni1[i][j] = fni0[i][j];
	}
      td_mem += timer_stop(t_1); // stop measure td_mem
	  
      // high-res file output of local traces (optional, frequent i/o slows down)
      // for time-trace analysis (e.g. to compute frequency spectra from pe(t) at some x,y)
      // /*
      if (printtraces) {
      t_1 = timer_start(); // start measure td_out

	itrace++; 
	if (itrace >= ntrace) {
	  // ft = fopen("traces.dat","a");
	  // fprintf(ft,"%.6f  %.6e  %.6e  \n", ttt, ne[iout][jout], pe[iout][jout]); 
	  itrace = 0;
	  for (k=0; k<k_max; k++) {
	    trace_n[k][jtrace] = ne[iout[k]][jout[k]];
	    trace_p[k][jtrace] = pe[iout[k]][jout[k]];
	  }
	  // fclose(ft);
	  jtrace ++;
	}
	td_out += timer_stop(t_1); // stop measure td_out
      }
      // */
	  
    } // ... end inner time loop ..............................................


      // low-time resolution diagnostics: energies and snapshot outputs
    t_1 = timer_start(); // start measure td_out
    diagnose(ttt,t00,it,itmax,nn,ne,ni,pe,gni);
    td_out += timer_stop(t_1); // stop measure td_out

  } // ... end outer time loop ..........................................

  // End time step --------------------------------------------------------
  
  t_1 = timer_start(); // start measure td_out
 
  // write complete restart file
  g = fopen( "restart.dat", "w" );
  fprintf(g,"%5e  \n", ttt);
  for (i=0; i<=nx1; i++)
    for (j=0; j<=ny1; j++)
      fprintf(g,"%d  %d  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e\n",
	      i,j,eno,pe[i][j],ne[i][j],ne1[i][j],ne2[i][j],ni[i][j],ni1[i][j],ni2[i][j],
	      fne1[i][j],fne2[i][j],fni1[i][j],fni2[i][j]);
  fclose( g );

   
  // Time trace / frequency spectrum analysis:
  if (printtraces) {
    printf("|\n| computing sum over %d frequency spectra ... \n", k_max);

    double pt[maxtrace],  kt_sum[maxtrace/2], func_window[maxtrace];
    // double dx02 = double(maxtrace/10)*double(maxtrace/10);  
    double dx0 = double(maxtrace)/10.;
    fftw_complex kt[maxtrace/2];
    fftw_plan hindftt;
    hindftt = fftw_plan_dft_r2c_1d(maxtrace, &pt[0], &kt[0], FFTW_ESTIMATE);
    
    for (i=0; i<maxtrace; i++) {
      func_window[i] = 1.;
      func_window[i] -= exp(-(double(i)/dx0)*(double(i)/dx0)) + exp(-(double((maxtrace-i)/dx0)*(double(maxtrace-i)/dx0)));
    }
    
    // Traces output:
    
    ft = fopen("traces.dat","w");
    for (i=0; i<maxtrace; i++)
      // fprintf(ft,"%d  %.6e  %.6e  %.6e \n",
      //     i, trace_n[1][i]*func_window[i], trace_p[1][i]*func_window[i], func_window[i]); 
      fprintf(ft,"%.6e  %.6e  %.6e  %.6e \n",
    	      (t00 + double(i*ntrace)*dt), trace_n[1][i]*func_window[i], trace_p[1][i]*func_window[i], func_window[i]); 
    fclose(ft);
    
    // Fourier frequency spectra for ne and phi:
    
    for (i=0; i<=maxtrace/2; i++) kt_sum[i] = 0.;
    for (k=0; k<k_max; k++) {
      for (i=0; i<=maxtrace; i++) pt[i] = trace_n[k][i]*func_window[i];
      fftw_execute(hindftt);
      for (i=0; i<=maxtrace/2; i++) kt_sum[i] += fabs(kt[i][0]);
    }
    ft = fopen( "ptt_n.dat", "w" );      
    for (i=1; i<maxtrace/2; i++) 
      fprintf(ft,"%.6e  %.6e\n",double(i)*TwoPi/(double(itmax)*dt),kt_sum[i]/double(k_max*maxtrace));
    fclose( ft );
    
    
    for (i=0; i<=maxtrace/2; i++) kt_sum[i] = 0.;
    for (k=0; k<k_max; k++) {
      for (i=0; i<=maxtrace; i++) pt[i] = trace_p[k][i]*func_window[i];
      fftw_execute(hindftt);
      for (i=0; i<=maxtrace/2; i++) kt_sum[i] += fabs(kt[i][0]);
    }      
    ft = fopen( "ptt_p.dat", "w" );
    for (i=1; i<maxtrace/2; i++) 
      fprintf(ft,"%.6e  %.6e\n",double(i)*TwoPi/(double(itmax)*dt),kt_sum[i]/double(k_max*maxtrace));
    fclose( ft );
    fftw_destroy_plan(hindftt);
  }

  td_out += timer_stop(t_1); // stop measure td_out

  // td_out += t_2 - t_1;
  td_tot = td_init + td_update + td_pois1 + td_pois2 + td_pol + td_bnd + td_out + td_mem;

  printf("|\n| total run time:  %.2f s  =  %.2f min \n", td_tot, td_tot/60.);
 
  td_tot = 100./td_tot;
  
  // total times bottleneck analysis output:
  printf("|\n| absolute run times: ini: %.2f upd: %.2f poi1: %.2f poi2: %.2f pol: %.2f bnd: %.2f out: %.2f mem:%.2f\n",
	 td_init, td_update, td_pois1, td_pois2, td_pol, td_bnd, td_out, td_mem);
  
  // percentual bottleneck analysis output:
 printf("| relative run times: ini: %.2f upd: %.2f poi1: %.2f poi2: %.2f pol: %.2f bnd: %.2f out: %.2f mem:%.2f\n|\n",
	 td_init*td_tot, td_update*td_tot, td_pois1*td_tot, td_pois2*td_tot, td_pol*td_tot, td_bnd*td_tot, td_out*td_tot, td_mem*td_tot);
 
#ifdef _OPENMP 
  fftw_cleanup_threads();  
#endif
  
  printf("| GHW1 end.\n\n");
}
// END OF MAIN ================================================================

// Specify functions:

// ---------------------------------------------------------------------------
// Read input data file "ghw1.inp" and initialise some more parameters
void init_parameters(void)
{ 
  int i,j, ith; char s[80]; FILE *g, *f;
  double rel, para[50];

  for (int i=0;i<=40;i++) printf("_");
  printf("\n| GHW1 start...\n");

  printf("| initialisation ... \n");

  // this sequentially reads all numbers appearing behind any "=" sign in ghw1.inp
  // caution: be aware of changing the numbering in case of additions to ghw1.inp
  g = fopen( "ghw1.inp", "r" );
  i = 0;
  while (!feof(g)) { fscanf(g, "%s ", s); if (strstr(s,"=")) i++; }; 
  rewind(g);
  for (j=1; j<=i; j++) {
    while (!strstr(s,"=")) fscanf(g,"%s ",s); 
    fscanf(g,"%s ",s); para[j]=atof(s); 
  }; 
  fclose( g );

  chat  = para[1];   // nonadiabatic HW coupling parameter (best choice: 0.01 y chat < 10)
                     // ( = "alpha" in Numata et al, and = C in Camargo et al.)
  b_mhw = true;
  if (chat<0.) {chat = fabs(chat); b_mhw = false;};  // ordinary Hasegawa-Wakatani
  
  taui  = para[2];   // ion temperature ratio taui = T_i / ( Z_i * T_e )
  xlwl = 1.; if (taui<0.) { taui = fabs(taui); xlwl = 0.; };
  
  delta = para[3];   // drift scale rho_s / Lperp with Lperp = (d/dx) ln no
                     // (for example: delta = 1/64 = 0.0156)
                     // sets time scale to t*cs/Lperp (as in Camargo et al)
                     // (Numata et al. use t*omega_ci normalisation, so that delta==1)

  diff  = para[4];   // perp. hyperviscosity strength, scales non-trivially with grid resolution
                     // guideline: when doubling resolution, then divide diff by 8.
  
  incon = para[5];   // initial condition: choose between blob, turbulence, vortex, shear flows
                     // (turbulence develops from blob, but faster from "turbulent" bath)
  
  sigma = para[6];   // initial blob or vortex width

  amp   = para[7];   // initial perturbation amplitude (of blob, vortex or "turbulent" bath)

  vorfree = para[8]; // vorticity free blob (for special studies, usually not relevant)

  ly = para[9];       // length of y domain in units of rho_s
  if (incon==5) ly = 1.*TwoPi/sigma;
  
  nx = int(para[10]); // grid points in x
  ny = int(para[11]); // grid points in y

  hy = double(ny)/ly; // finite-differencing factor, assumes square grid
  wsrc = sigma*double(ny)/ly;
  aspect = double(nx)/double(ny); // x-y box aspect ratio (only used in bath initialisation)
  clin = 1.;          // 1 for ExB nonlinearity, 0 for linearisation
  delinv = 1./delta;
  
  dt    = para[12];       // time step
  itstp = int(para[13]);  // inner steps between outputs
  itmax = int(para[14]);  // total steps until end
  ntrace= int(para[15]);  // steps between fast trace output
  maxtrace = int(para[14]/para[15]); // number of recorded trace steps
  if (ntrace==0) maxtrace = 1;

  
  ipade = int(para[16]);  // Pade approximation or Bessel function for gyro operators
  ihype = int(para[17]);  // order and type of hyperviscosity (p=2 is usual choice)
  fflag = FFTW_ESTIMATE;  // FFTW plan estimate or measure
  if (para[18]==1.) fflag = FFTW_MEASURE; 
  if (para[18]==2.) fflag = FFTW_PATIENT; 
  npar  = int(para[19]);  // number of parallel OpenMP threads (depending on available cores)

  b_4th = true; if (para[20] == 0.) b_4th = false; // 2nd order or fourth order Arakawa scheme
  
  // plasma species parameters
  // (change only when more species are added to code)
  taue= -1.;              // fixed electron temperature tau_e = T_e / ( Z_e * T_e ) == -1
  mue =  0.;              // in general mue = m_e/(Z_e*m_i) is very small, can be set to zero
  // mue =- 1./(2.*1836);    // electron mass ratio for deuterium main ion species
  mui = +1.;              // fixed ion mass ratio mui = m_i/(Z_i*m_i) == +1
  aae = -1.;              // fixed electron density ratio
  aai = +1.;              // fixed ion density ratio
  zze = -1.;              // fixed electron charge
  zzi = mui/fabs(mui); if (mui==0.) zzi=0.; // fixed ion charge

  // Initialise Open_MP ...................
#ifdef _OPENMP 
  int np, kp, npmax;
  ith = fftw_init_threads();
  if (ith==0) { printf("| ith=%d : thread init failed... exiting\n",ith); exit(1);  }
  npmax = omp_get_max_threads();
  if (npmax<npar)
    { printf("| npar = %d > %d max. available threads ... exiting\n",npar,npmax); exit(1);  }  
  omp_set_num_threads(npar);
  printf("| parallel threads ");
#pragma omp parallel private(kp)
  { 
    np = omp_get_num_threads();
    kp = omp_get_thread_num();
    printf("%d/",kp);
  } 
  printf(" of %d active.\n",np);
  double dnp=double(np);
#else
  printf("| single core processing\n");
  double dnp=1.;
#pragma omp barrier
#endif
  
  // Initialise parameters ..................

  printtraces = (ntrace>0) ? true : false;

  nxh = nx/2;
  nyh = ny/2; 
  nx1 = nx-1;
  ny1 = ny-1;

  hy2 = hy*hy;
  hysq = hy2/12.;
  g_n = 1.;  // background density: g_n ~ inverse gradient length 1/Lp (= kappa in Numata)
  if (chat==0.) g_n = 0.;

  xyz = 1./double(nx*ny);
  
  if (incon!=2.) {
    f = fopen("eng.dat","w");    fclose(f);
    f = fopen("eng-transfer.dat","w");    fclose(f);
    f = fopen("fne.dat","w");    fclose(f);
    f = fopen("gamma.dat","w");  fclose(f);
    f = fopen("traces.dat","w"); fclose(f);
    f = fopen("time.dat","w");   fclose(f);
    f = fopen("zfx.dat", "w");   fclose(f);
  }

  for (i=1; i<=nx/2; i++) { pkxavg[i] = 0;}
  for (j=1; j<=ny/2; j++) { pkyavg[j] = 0; pkyavgn[j] = 0; pkyavgw[j] = 0;  pkyavge[j] = 0;}

  // coefficients pre-calculation for FFTW Poisson solver
  int ik, jk;
  double kx, ky, dxy = 1./(hy*hy);
  double cnx = 1./double(nx);
  double cny = 1./double(ny);
  double cxy = cnx*cny;
  double kkqq, cinv, bki, bke, taui0;
  double gambese, gambesi;
  double G0i, G0e, G1i, G1e;

  taui0 = (taui==0.) ? 1.e-9 : taui;
  double mue0 = (mue==0.) ? 1.e-9 : mue;

  // for (int i=0; i<=nx1; ++i) for (int j=0; j<=ny1; ++j) {
  for (int i=0; i<=nx1; ++i) for (int j=0; j<=nyh+1-1; ++j) {
      ik = (i>=nxh) ? i - nx : i;
      jk = (j>=nyh) ? j - ny : j;
      
      // k_perp^2:
      kkqq = (TwoPi*hy*ik*cnx)*(TwoPi*hy*ik*cnx) + (TwoPi*hy*jk*cny)*(TwoPi*hy*jk*cny) ;

      bki = taui0*mui*kkqq;
      bke = taue*mue0*kkqq;

      // Gamma operators (gyro-averaging and gyro-screening):
      if (ipade==0) { // Pade approximation
	G0i = 1./(1.+xlwl*bki); G1i = 1./(1.+0.5*bki);
	G0e = 1./(1.+xlwl*bke);	G1e = 1./(1.+0.5*bke);
	cinv = kkqq*( aai*mui*G0i + aae*mue*G0e );
      }
      else if (ipade==1) { // Bessel function
	G0i = gamma0(bki); G1i = sqrt(G0i);
	G0e = gamma0(bke); G1e = sqrt(G0e);
	cinv = (aai/taui0)*(1.-G0i) + (aae/taue)*(1.-G0e) ;
      }


      cvort[i][j] = - cxy/cinv;
      cpoti[i][j] = + cxy*G1i;
      cvfis[i][j] = + cxy/G1i; 

      // Laplace
      clap[i][j] = - cxy*kkqq;
		
      // hyper k^2 = laplace:
      if (ihype==1) chyv[i][j] = +cxy*kkqq; // Navier-Stokes stype viscosity
      // hyper k^4 = normal:
      if (ihype==2) chyv[i][j] = cxy*kkqq*kkqq; // k^4 "usual" hyperviscosity
      // hyper k^8 = enhanced:
      if (ihype==4) chyv[i][j] = cxy*kkqq*kkqq*kkqq*kkqq;
      // hyper k^8 = enhanced:
      if (ihype==8) chyv[i][j] = cxy*kkqq*kkqq*kkqq*kkqq*kkqq*kkqq;

    }
  
  cvort[0][0] = 0.;

  // Test Bessel output: compare with Fig.3, Dorland & Hammett 1993.
  /*
  g = fopen( "bessel.dat", "w" );
  for (int i=0; i<=nxh; ++i) 
    {
      ik = (i>=nxh) ? i - nx : i; 
      kkqq = (TwoPi*hy*ik*cnx)*(TwoPi*hy*ik*cnx);  
      bki = taui*mui*kkqq;
      bke = taue*mue*kkqq;
      gambesi = gamma0(bki);
      gambese = gamma0(bke);
      
      fprintf(g,"%.8e  %.8e  %.8e  %.8e  %.8e  %.8e\n",sqrt(kkqq),
      	      1./(1.+.5*taui*mui*kkqq), sqrt( gambesi ), (kkqq*(aai*mui/(1.+taui*mui*kkqq)) + 1.*(aae*mue/(1.+bke))), ((aai/taui)*(1.-gambesi) + 1.*(aae/taue)*(1.-gambese)), gambesi );
    }
  fclose( g );
  // */
} // end init_parameters


// ---------------------------------------------------------------------------
// initialise single blob perturbation:
void init_add_blob( AXY ne, AXY ni )
{
  int i, j, ixloc;
  double dr = wsrc*wsrc;
    
  printf("| blob initial condition...\n");
  ixloc = nxh;
#pragma omp parallel for private(i,j) 
  for (i=0; i<=nx1; ++i) for (j=0; j<=ny1; ++j) {
      ne[i][j] += amp*( exp(-(i-ixloc)*(i-ixloc)/dr)*exp(-(j-nyh)*(j-nyh)/dr) );
      ni[i][j] = ne[i][j]; 
    }
  if (vorfree==1.) { poisson(ne,cvfis,ni); }  
} // end init_add_blob


// ---------------------------------------------------------------------------
// initialise pseudo-turbulent bath perturbation:
void init_add_turb( AXY ne, AXY ni )
{
  int i, j, ik, jk;
  int nk = 32, mmm = int(aspect);
  double zuf, kk11, turb, ppk, lx = ly*aspect, phs[nk+1][nk+1], cbath[nx][ny];
  fftw_complex wwk[nx][ny], ppx[nx][ny];
  fftw_plan plan_ini_back; 
 
  printf("| bath initial condition...\n");

  // Spectrum for pseudo-turbulent bath initialization:
  double ccc = 1./(lx*ly);    
  for (i=0; i<=nx-1; ++i) for (j=0; j<=ny-1; ++j) {
      ik = (i>=nx/2) ? i - nx : i;
      jk = (j>=ny/2) ? j - ny : j;	  
      kk11 =  TwoPi*TwoPi*(  ( double(ik*ik)/(lx*lx) ) + ( double(jk*jk)/(ly*ly) )  );	  
      cbath[i][j] =  ccc*kk11/(1. + kk11*kk11*kk11 ); 	  
    }    
  
  // use same pseudo "random" seed if reproducibility is desired!
  
#ifdef _OPENMP
  fftw_plan_with_nthreads(npar);
#endif
  plan_ini_back=fftw_plan_dft_2d(nx,ny,&wwk[0][0],&ppx[0][0],FFTW_BACKWARD,FFTW_ESTIMATE);

  // srand(time(NULL) + getpid()); // do not use random for reproducible runs!
 
  // (omp parallel on rand not reproducible)
  for (i=0; i<=nx1; ++i) for (j=0; j<=ny1; ++j) {
      zuf = TwoPi*0.05*(rand() % 20); 
      wwk[i][j][0] = amp*cbath[i][j]*cos(zuf);
      wwk[i][j][1] = amp*cbath[i][j]*sin(zuf);
    }
  fftw_execute(plan_ini_back);
  fftw_destroy_plan(plan_ini_back); 

#pragma omp parallel for private(j)
  for (i=0; i<=nx1; ++i) for (j=0; j<=ny1; ++j) {
      ne[i][j] += ppx[i][j][0]*delinv;
      ni[i][j] += ppx[i][j][0]*delinv;
    }

  if (vorfree==1.) { poisson(ne,cvfis,ni);  }

} // end init_add_turb


// ---------------------------------------------------------------------------
// initialise dual vortex perturbation:
void init_add_dual( AXY ne, AXY ni )
{
  int i, j, ixloc;
  double dr = wsrc*wsrc;
  double elong = 1.0; // e.g. circular for elong=1
  double dry = elong*dr;
  double xxx1, xxx2, vorset;  
  printf("| merger initial condition...\n");

  vorset = 0.01; // to define vorticity from gyrodensity difference
  g_n = 0.; // no background density gradient
  ixloc = nxh;
#pragma omp parallel for private(i,j,pdist,xxx1,xxx2) 
  for (i=0; i<=nx1; ++i) for (j=0; j<=ny1; ++j) {
      pdist = 4.0*wsrc; // set here distance, e.g. 4*wsrc
      xxx1 = (i-ixloc -0.5*pdist);
      xxx2 = (i-ixloc +0.5*pdist);
      ne[i][j] += +amp*( exp(-xxx1*xxx1/dr)*exp(-(j-ny1*1/2)*(j-ny1*1/2)/dry) );
      ne[i][j] += +amp*( exp(-xxx2*xxx2/dr)*exp(-(j-ny1*1/2)*(j-ny1*1/2)/dry) );
    }
  for (i=0; i<=nx1; ++i) for (j=0; j<=ny1; ++j) ni[i][j] = (1.-vorset)*aai*ne[i][j];
  
} // end init_add_dual


// ---------------------------------------------------------------------------
// initialise shear flow perturbation:
void init_add_flow( AXY ne, AXY ni )
{
  int i, j, ixloc = nxh;
  printf("| shear flow initial condition...\n");

  g_n = 0.; chat = 0.; // no background density gradient or HW drive
#pragma omp parallel for private(i,j) 
  for (i=0; i<=nx1; ++i) for (j=0; j<=ny1; ++j) {
      ne[i][j] = amp*cos(TwoPi*1.*double(i-nxh)/double(nx));
      ne[i][j]+= amp*0.001*sin(TwoPi*5.*double(j)/double(ny));
      ni[i][j] = 0.;
    }
} // end init_add_flow


// ---------------------------------------------------------------------------
// initialise linear single ky mode
void init_add_mode( AXY ne, AXY ni )
{
  int i, j;
  printf("| single ky mode initial condition...\n");

  // Here par[6]=sigma is the mode number (rho_s ky).
  // Ly is set to m*TwoPi*sigma to get periodicity (done in initialisation; m=1, ...)!
  clin = 0.;
#pragma omp parallel for private(i,j) 
  for (i=0; i<=nx1; ++i) for (j=0; j<=ny1; ++j) {
      ne[i][j] = amp*sin(sigma*double(j)/hy);
      // ne[i][j]*= sin(TwoPi*0.5*double(i)/double(nx)); // (for kx!=0)
      ni[i][j] = 0.;
    }

  // Write analytical dispersion relation for kx=0 (cf. Camargo 95)
  // for comparison with frequency and growth rate calculated in "diagnose":
  /*
  double ky, a_omr, a_gam, a_sig, a_lam0, a_lam1, a_aaa, a_bbb;
  f = fopen("formula-gamma.dat", "w");
  for (i=1; i<=1000; ++i) {
    ky = 0.01*double(i);
    a_lam0 = chat*(1.+ky*ky)/(ky*ky);
    a_lam1 = 2.*diff*pow(ky,4.);
    a_sig = chat/ky;
    a_aaa = a_lam0/(sqrt(2.)*2.);
    a_bbb = sqrt(1.+16.*a_sig*a_sig/(a_lam0*a_lam0*a_lam0*a_lam0));
    a_omr = a_aaa*sqrt(a_bbb - 1.);
    a_gam = - .5*(a_lam0+a_lam1) + a_aaa*sqrt(a_bbb + 1.);
    fprintf(f,"%.3e  %.3e  %.3e \n",ky,a_omr,a_gam);
  }
  fclose(f);
  */
} // end init_add_mode


// ---------------------------------------------------------------------------
// calculate and write output quantities:
// e.g. energies, transport, k spectra, 2D arrays, profiles, ...
void diagnose( double ttt, double t00, int it, int itmax,
	       AXY n_n, AXY n_e, AXY n_i, AXY p_p, AXY n_g )
{
  int i,j, ik, im,ip,jm,jp;
  double enn, enp, enw, eeb, ezf, fne, fnsol, zfx, dum, rey, freq, grow;
  double etran_grad, etran_adia, etran_visc;
  double lap[nxm][nym]; 
  double avg[nxm], avp[nxm];
  FILE *f, *g, *h, *g1, *g2, *g3, *g4, *g5, *g6, *g7, *g8;

#pragma omp parallel for private(j) shared(avg,avp)
  for (i=0; i<=nx1; i++) // zonal averages of ne and phi
    { avg[i] = 0.; avp[i] = 0.; 
      for (j=0; j<=ny1; j++) {avg[i] += n_e[i][j]; avp[i] += p_p[i][j];}
      avg[i]/=ny; avp[i]/=ny;
    }
  enn = 0.; enp = 0.; enw = 0.; eeb = 0.; ezf = 0.; fne = 0.; 
  etran_grad = 0.; etran_adia= 0.; etran_visc = 0.;
  
  // vorticity
  // laplace(p_p,lap);
  poisson(p_p,clap,lap); 

  // energetic quantities and transport
#pragma omp parallel for private(j,ip,dum) reduction(+:enn,ezf,eeb,fne,enw,etran_grad,etran_adia,etran_visc)
  for (i=0; i<=nx1; i++) {
    ip = (i==nx1) ? 0   : i+1;
      
    for (j=0; j<=ny1; j++) {

      // total thermal free energy:
      enn += n_e[i][j]*n_e[i][j] + taui*n_i[i][j]*n_i[i][j];

      // (mostly zonal) flow energy:
      dum = hy*(avp[ip]-avp[i]);
      ezf += dum*dum;
	  
      // turbulent enstrophy:
      enw += lap[i][j]*lap[i][j];

      // total kinetic energy:
      dum = DFDY(p_p,i,j);
      eeb += dum*dum; 
      dum = DFDX(p_p,i,j);
      eeb += dum*dum; 
      
      // electron ExB particle transport:
      fne += - n_e[i][j]*DFDY(p_p,i,j);

      // gradient drive / transfer:
      etran_grad+= - n_e[i][j]*DFDY(pe,i,j); 
      etran_grad+= - n_i[i][j]*DFDY(pi,i,j)*taui;
   
      // "parallel" resistive dissipation sink / transfer:
      dum = pe[i][j] - n_e[i][j];
      if (b_mhw) dum-= avp[i] - avg[i];
      etran_adia+= -chat*dum*dum;
      
      // perpendicular hyper-viscous sink / transfer:
      dum = hyve[i][j]*(pe[i][j] - n_e[i][j]);
      dum-= hyvi[i][j]*(pi[i][j] + taui*n_i[i][j]);
      etran_visc+= diff*dum;
   
    }
  } 
  
  enn *= .5*xyz; enp *= xyz; enw *= .5*xyz; eeb *= .5*xyz; ezf *= .5*xyz; 
  if (enw==0.) enw=1.e-12; if (ezf==0.) ezf=1.e-12; if (eeb==0.) eeb=1.e-12;

  eeb*=delinv*delinv; ezf*=delinv*delinv; enw*=delinv*delinv; enn*=delinv*delinv;
  fne*= xyz*delinv*delinv;

  etran_grad*= xyz; etran_adia*= xyz; etran_visc*= xyz;
  double etran_tot = etran_grad + etran_adia + etran_visc;
  
  // crash control: stop if nan or inf
  // if ( (isnan(enn)) || (enn<1.e-16) )  
  if (isnan(enn))
    {
      printf("\n|t=%.3f (%.2f): enn=%.3e\n", ttt,(t00+itmax*dt),enn);
      printf("|HGW1 CRASHED!  enn=%.3e,   eno=%.3e\n", enn, eno);
#ifdef _OPENMP 
      fftw_cleanup_threads();
#endif
      exit(1);
    }
  eno = enn;

  // global energy time series output:
  f = fopen("eng.dat","a");
  if (ttt>0.) fprintf(f,"%.8f  %.6e  %.6e  %.6e  %.6e  %.6e\n",
		      ttt, enn, eeb, ezf, (enn+eeb), ezf/eeb);
  fclose(f);
  // global energetic transfer time series output:
  f = fopen("eng-transfer.dat","a");
  if (ttt>0.) fprintf(f,"%.8f  %.6e  %.6e  %.6e  %.6e\n",
		      ttt, etran_grad, etran_adia, etran_visc, etran_tot);
  fclose(f);
  // global transport time series output:
  f = fopen("fne.dat","a");
  if (ttt>0.) fprintf(f,"%.8f  %.6e\n", ttt, fne);
  fclose(f);

  // (non)linear growth rate "grow" and linear frequency estimate "freq":
  if (incon==5.) {
    // compute the (approximate) linear drift wave frequency freq

    double dt_amp = 0.;
    freq = 0.; grow = 0.;
    jsgn_old = jsgn;
    if (p_p[nxh][nyh]==fabs(p_p[nxh][nyh])) {jsgn = 1;} else {jsgn = -1;};
    
    if (jsgn != jsgn_old) {t_old = t_amp; t_amp = ttt;};
    dt_amp = t_amp - t_old;
    if ( (dt_amp !=0.) && (dt_amp != ddtt) && (it>2) )  freq = 0.5*TwoPi/dt_amp;

  }

  if (enwo == 1.e-12) enwo = enw;
  grow = (sqrt(enw)-sqrt(enwo))/(sqrt(enwo)*ddtt);
  enwo = enw;
  f = fopen("gamma.dat","a");
  if (it>1) fprintf(f,"%.3f  %.6e  %.6e\n", ttt, freq, grow);
  fclose(f);
  
  // Plot outputs

  // zonal flow (t,x) 2D plot
  h = fopen( "zfx.dat", "a" );
  for (i=0; i<=nx1; i++)  {
    im = (i==0)   ? nx1 : i-1;
    ip = (i==nx1) ? 0   : i+1;
    for (zfx=0., j=0; j<=ny1; j++)  {
      // zfx += p_p[i][j]/(ny); // zonal potential
      zfx += p_p[ip][j]-p_p[im][j]; // zonal flow <Vy>
    }
    zfx *= .5*hy/(ny); // comment out for pot or vor.
    if (ttt>0.) fprintf(h,"%.3f  %d  %.6e \n", ttt,i,zfx);
  } 
  fprintf(h,"\n");
  fclose(h);

  // x-cuts at y=ny/2:
  g = fopen( "cutx.dat", "w" );
  for (i=0; i<=nx1; i++)
    for (j=nyh; j<=nyh; j++)
      fprintf(g,"%f  %.6e   %.6e   %.6e   %.6e   %.6e\n",
	      double(i)/hy,n_e[i][j],n_i[i][j], p_p[i][j],lap[i][j],n_g[i][j]);
  fclose( g );
  
  g = fopen( "cuty.dat", "w" );
  for (j=0; j<=ny1; j++)
      fprintf(g,"%f  %.6e   %.6e   %.6e   %.6e   %.6e\n",
	      double(j)/hy,n_e[nxh][j],n_i[nxh][j], p_p[nxh][j],lap[nxh][j],n_g[nxh][j]);
  fclose( g );

  // x-profiles averaged over y:
  g = fopen( "xprof.dat", "w" );
  double peprof, neprof, niprof, nsprof, wwprof;
  double pexprof[nxm];
  for (i=0; i<=nx1; i++) {
    peprof = 0.; neprof = 0.; niprof = 0.; nsprof = 0.; wwprof=0.; rey = 0.;
    im = (i==0)   ? nx1 : i-1;
    ip = (i==nx1) ? 0   : i+1;
    for (j=0; j<=ny1; j++) {
      jm = (j==0)   ? ny1 : j-1;
      jp = (j==ny1) ? 0   : j+1;
      // Reynolds stress:
      rey+= .25*hy*hy*(p_p[ip][j]-p_p[im][j])*(p_p[i][jp]-p_p[i][jm]);
      // potential, densities, vorticity
      peprof += p_p[i][j];
      neprof += n_e[i][j];
      niprof += n_i[i][j];
      wwprof += lap[i][j];
    }
    peprof/=ny; neprof/=ny; niprof/=ny; nsprof/=ny; wwprof/=ny; rey/=ny;
    fprintf(g,"%f  %.6e   %.6e   %.6e   %.6e   %.6e\n",
	    double(i)/hy,peprof,neprof,niprof,wwprof,rey);
  }
  fclose( g );

  // 2D (x,y) plots of densities, vorticity, potential
  g1 = fopen( "n2d0.dat", "w" ); g2 = fopen( "w2d0.dat", "w" ); 
  g3 = fopen( "p2d0.dat", "w" ); g4 = fopen( "i2d0.dat", "w" ); 
  for (i=0; i<=nx1; i++) {
    for (j=0; j<=ny1; j++) {
      fprintf(g1,"%d  %d  %.6e\n",i,j,n_e[i][j]); 
      fprintf(g2,"%d  %d  %.6e\n",i,j,lap[i][j]); 
      fprintf(g3,"%d  %d  %.6e\n",i,j,p_p[i][j]); 
      fprintf(g4,"%d  %d  %.6e\n",i,j,(p_p[i][j]-n_e[i][j])); 
    }
    fprintf(g1,"\n"); fprintf(g2,"\n"); fprintf(g3,"\n"); 
    fprintf(g4,"\n"); 
  }
  fclose( g1 ); fclose( g2 ); fclose( g3 ); fclose( g4 );
  rename("n2d0.dat", "n2d.dat"); rename("w2d0.dat", "w2d.dat"); 
  rename("p2d0.dat", "p2d.dat"); rename("i2d0.dat", "i2d.dat"); 

  // Fourier ky spectra
  double py[ny];
  // fftw_complex ky[ny/2], sumky[ny/2];
  fftw_complex ky[ny/2+1], sumky[ny/2+1];
#ifdef _OPENMP 
  fftw_plan_with_nthreads(npar);
#endif
  fftw_plan hindfty;
  hindfty = fftw_plan_dft_r2c_1d(ny, &py[0], &ky[0], FFTW_ESTIMATE);

  // potential phi
  for (j=0; j<=ny/2+1; j++) sumky[j][0] = 0.;
  for (i=0; i<=nx1; i++)
    {
      for (j=0; j<=ny1; j++) py[j] = p_p[i][j];
      fftw_execute(hindfty);
      // for (j=0; j<=ny/2; j++) sumky[j][0] += ky[j][0]*ky[j][0]; 
      for (j=0; j<=ny/2+1; j++) sumky[j][0] += fabs(ky[j][0]); 
    }
  for (j=0; j<=ny/2+1; j++) sumky[j][0] /= double(nx); 
  if (incon==2)
    { for (j=1; j<=ny/2; j++) pkyavg[j] = ( (it-1+1)*pkyavg[j] + sumky[j][0] ) / (it+1); }
  g = fopen( "pky_p.dat", "w" );
  for (j=1; j<ny/2; j++) 
    fprintf(g,"%.6e  %.6e  %.6e\n",TwoPi*double(j)/ly,sumky[j][0]*ny/hy/hy,pkyavg[j]*ny/hy/hy);
  fclose( g );
  
  // density ne
  for (j=0; j<=ny/2+1; j++) sumky[j][0] = 0.;
  for (i=0; i<=nx1; i++)
    {
      for (j=0; j<=ny1; j++) py[j] = ne[i][j];
      fftw_execute(hindfty);
      // for (j=0; j<=ny/2; j++) sumky[j][0] += ky[j][0]*ky[j][0]; 
      for (j=0; j<=ny/2+1; j++) sumky[j][0] += fabs(ky[j][0]); 
    }
  for (j=0; j<=ny/2+1; j++) sumky[j][0] /= double(nx); 
  if (incon==2)
    { for (j=1; j<=ny/2; j++) pkyavgn[j] = ( (it-1+1)*pkyavgn[j] + sumky[j][0] ) / (it+1); }
  g = fopen( "pky_n.dat", "w" );
  for (j=1; j<ny/2; j++) 
    fprintf(g,"%.6e  %.6e  %.6e\n",TwoPi*double(j)/ly,sumky[j][0]*ny/hy/hy,pkyavgn[j]*ny/hy/hy);
  fclose( g );

  // vorticity w
  for (j=0; j<=ny/2+1; j++) sumky[j][0] = 0.;
  for (i=0; i<=nx1; i++)
    {
      for (j=0; j<=ny1; j++) py[j] = lap[i][j];
      fftw_execute(hindfty);
      for (j=0; j<=ny/2+1; j++) sumky[j][0] += ky[j][0]*ky[j][0]; 
    }
  for (j=0; j<=ny/2+1; j++) sumky[j][0] /= double(nx); 
  if (incon==2)
    { for (j=1; j<=ny/2; j++) pkyavgw[j] = ( (it-1+1)*pkyavgw[j] + sumky[j][0] ) / (it+1); }
  g = fopen( "pky_w.dat", "w" );
  for (j=1; j<ny/2; j++) 
    fprintf(g,"%.6e  %.6e  %.6e\n",TwoPi*double(j)/ly,.5*sumky[j][0]*ny/hy/hy,.5*pkyavgw[j]*ny/hy/hy);
  fclose( g );

  // kinetic energy E
  for (j=0; j<=ny/2+1; j++) sumky[j][0] = 0.;
  for (i=0; i<=nx1; i++)  {
    im = (i==0)   ? nx1 : i-1;
    ip = (i==nx1) ? 0   : i+1;
    for (j=0; j<=ny1; j++) {
      jm = (j==0)   ? ny1 : j-1;
      jp = (j==ny1) ? 0   : j+1;
      dum = .5*hy*(p_p[i][jp]-p_p[i][jm]);
      py[j] = dum*dum;
      dum = .5*hy*(p_p[ip][j]-p_p[im][j]);
      py[j]+= dum*dum;
    }
    fftw_execute(hindfty);
    for (j=0; j<=ny/2+1; j++) sumky[j][0] += ky[j][0]*ky[j][0]; 
  }
  for (j=0; j<=ny/2+1; j++) sumky[j][0] /= (nx); 
  if (incon==2)
    { for (j=1; j<=ny/2; j++) pkyavge[j] = ( (it-1+1)*pkyavge[j] + sumky[j][0] ) / (it+1); }
  g = fopen( "pky_e.dat", "w" );
  for (j=1; j<ny/2; j++) 
    fprintf(g,"%.6e  %.6e  %.6e\n",TwoPi*double(j)/ly,.5*sumky[j][0]*ny/hy/hy,.5*pkyavge[j]*ny/hy/hy);
  fclose( g );
  fftw_destroy_plan(hindfty);
  
  // Fourier kx spectra
  double px[nx];
  fftw_complex kx[nx/2+1], sumkx[nx/2+1];
  fftw_plan hindftx;
  hindftx = fftw_plan_dft_r2c_1d(nx, &px[0], &kx[0], FFTW_ESTIMATE);
  for (i=0; i<=nx/2; i++) sumkx[i][0] = 0.;
  for (j=0; j<=ny1; j++)
    {
      for (i=0; i<=nx1; i++) px[i] = pe[i][j];
      fftw_execute(hindftx);
      for (i=0; i<=nx/2+1; i++) sumkx[i][0] += fabs(kx[i][0]); 
    }
  for (i=0; i<=nx/2+1; i++) sumkx[i][0] /= (ny); 
  fftw_destroy_plan(hindftx);
  if (incon==2)
    { for (i=1; i<=nx/2+1; i++) pkxavg[i] = ( (it-1+1)*pkxavg[i] + sumkx[i][0] ) / (it+1); }
  g = fopen( "pkx_p.dat", "w" );
  for (i=1; i<nx/2; i++) 
    fprintf(g,"%.6e  %.6e  %.6e\n",TwoPi*double(i)/ly,sumkx[i][0]*nx/hy/hy,pkxavg[i]*nx/hy/hy);
  fclose( g );

#ifdef _OPENMP 
  fftw_cleanup_threads();
#endif
  // end spectral analysis
  
  // timestamp output
  printf("|t=%.3f (%.3f): enn=%.3e\n", ttt,(t00+itmax*dt),enn);

}

// -----------------------------------------------------------------------------
// 2nd order Arakawa scheme for brackets (see Arakawa 1966, reprint in JCP 1997)
// [u,v] = (du/dx)(dv/dy)-(du/dy)(dv/dx)
void arakawa( AXY uuu, AXY vvv, AXY out )
{ 
  int i00,j00,ip1,jp1,im1,jm1,im2,ip2,jm2,jp2;
  double aj2, afac = .25*hy*hy/3.;
  int i0 = 0, i1 = nx1, j0 = 0, j1 = ny1;

#pragma omp parallel for private(j00,im1,im2,jm1,jm2,ip1,ip2,jp1,jp2,aj2)
  for (i00=0; i00<=nx1; ++i00)  {
    nextp(i00,i0,i1,im1,ip1,im2,ip2); 
 
    for (j00=0; j00<=ny1; ++j00) {
      nextp(j00,j0,j1,jm1,jp1,jm2,jp2);

      // J++ :
      aj2  = (uuu[ip1][j00] - uuu[im1][j00]) * (vvv[i00][jp1] - vvv[i00][jm1]);
      aj2 -= (uuu[i00][jp1] - uuu[i00][jm1]) * (vvv[ip1][j00] - vvv[im1][j00]);
      
      // J+x :
      aj2 += uuu[ip1][j00] * ( vvv[ip1][jp1] - vvv[ip1][jm1] );
      aj2 -= uuu[im1][j00] * ( vvv[im1][jp1] - vvv[im1][jm1] );
      aj2 -= uuu[i00][jp1] * ( vvv[ip1][jp1] - vvv[im1][jp1] );
      aj2 += uuu[i00][jm1] * ( vvv[ip1][jm1] - vvv[im1][jm1] );
      
      // Jx+ :
      aj2 += uuu[ip1][jp1] * ( vvv[i00][jp1] - vvv[ip1][j00] );
      aj2 -= uuu[im1][jm1] * ( vvv[im1][j00] - vvv[i00][jm1] );
      aj2 -= uuu[im1][jp1] * ( vvv[i00][jp1] - vvv[im1][j00] );
      aj2 += uuu[ip1][jm1] * ( vvv[ip1][j00] - vvv[i00][jm1] );
      
      // 2nd order part:
      out[i00][j00] = aj2*afac;      
    };
  }  
} // end function arakawa


// -----------------------------------------------------------------------------
// 4th order Arakawa scheme for brackets (see Arakawa 1966, reprint in JCP 1997)
// [u,v] = (du/dx)(dv/dy)-(du/dy)(dv/dx)
// void arakaw4(double (*uuu)[nym], double (*vvv)[nym], double (*out)[nym])
void arakaw4( AXY uuu, AXY vvv, AXY out )
{ 
  int i00,j00,ip1,jp1,im1,jm1,im2,ip2,jm2,jp2;
  double aj2, aj4, afac = .25*hy*hy/3.;
  int i0 = 0, i1 = nx1, j0 = 0, j1 = ny1;
  
#pragma omp parallel for private(j00,im1,im2,jm1,jm2,ip1,ip2,jp1,jp2,aj2,aj4)
  for (i00=0; i00<=nx1; ++i00)  {
    nextp(i00,i0,i1,im1,ip1,im2,ip2); 
  
    for (j00=0; j00<=ny1; ++j00) {
      nextp(j00,j0,j1,jm1,jp1,jm2,jp2);

      // 2nd order part:
      // J++ :
      aj2  = (uuu[ip1][j00] - uuu[im1][j00]) * (vvv[i00][jp1] - vvv[i00][jm1]);
      aj2 -= (uuu[i00][jp1] - uuu[i00][jm1]) * (vvv[ip1][j00] - vvv[im1][j00]);
      
      // J+x :
      aj2 += uuu[ip1][j00] * ( vvv[ip1][jp1] - vvv[ip1][jm1] );
      aj2 -= uuu[im1][j00] * ( vvv[im1][jp1] - vvv[im1][jm1] );
      aj2 -= uuu[i00][jp1] * ( vvv[ip1][jp1] - vvv[im1][jp1] );
      aj2 += uuu[i00][jm1] * ( vvv[ip1][jm1] - vvv[im1][jm1] );
      
      // Jx+ :
      aj2 += uuu[ip1][jp1] * ( vvv[i00][jp1] - vvv[ip1][j00] );
      aj2 -= uuu[im1][jm1] * ( vvv[im1][j00] - vvv[i00][jm1] );
      aj2 -= uuu[im1][jp1] * ( vvv[i00][jp1] - vvv[im1][j00] );
      aj2 += uuu[ip1][jm1] * ( vvv[ip1][j00] - vvv[i00][jm1] );

      aj2 *= afac;

      // 4th order part:
      // Jxx ;
      aj4  = (uuu[ip1][jp1] - uuu[im1][jm1]) * (vvv[im1][jp1] - vvv[ip1][jm1]);
      aj4 -= (uuu[im1][jp1] - uuu[ip1][jm1]) * (vvv[ip1][jp1] - vvv[im1][jm1]);
      
      // J4x+ :
      aj4 += uuu[ip1][jp1] * ( vvv[i00][jp2] - vvv[ip2][j00] );
      aj4 -= uuu[im1][jm1] * ( vvv[im2][j00] - vvv[i00][jm2] ); // wrong "+" in Arakawa 1966 !
      aj4 -= uuu[im1][jp1] * ( vvv[i00][jp2] - vvv[im2][j00] );
      aj4 += uuu[ip1][jm1] * ( vvv[ip2][j00] - vvv[i00][jm2] );
      
      // J4+x :
      aj4 += uuu[ip2][j00] * ( vvv[ip1][jp1] - vvv[ip1][jm1] );
      aj4 -= uuu[im2][j00] * ( vvv[im1][jp1] - vvv[im1][jm1] );
      aj4 -= uuu[i00][jp2] * ( vvv[ip1][jp1] - vvv[im1][jp1] );
      aj4 += uuu[i00][jm2] * ( vvv[ip1][jm1] - vvv[im1][jm1] );
      
      aj4 *= .5*afac;
      
      // Total (4th order):
      out[i00][j00] = 2.*aj2 - aj4 ;

      // Uncomment for test: 2nd order part only:
      // out[i00][j00] = aj2;
    };
  }  
} // end function arakaw4


// ---------------------------------------------------------------------------
// define 2-neighbourhood of coordinate i,j for 4th order finite differences:
inline void nextp (int i, int i0, int i1, int& im, int& ip, int& im2, int& ip2)
{
  im = (i==i0) ? i1 : i-1;  im2 = (i==i0+1) ? i1 : i-2;  if (i==i0) im2 = i1-1; 
  ip = (i==i1) ? i0 : i+1;  ip2 = (i==i1-1) ? i0 : i+2;  if (i==i1) ip2 = i0+1;
}

// ---------------------------------------------------------------------------
// 4th order y derivative:
inline double DFDY( AXY ain, int i, int j)
{ 
  int jm, jp, jm2, jp2;
  nextp(j,0,ny1,jm,jp,jm2,jp2);
  if (b_4th) { return ( hy*r12*( 8.*( ain[i][jp]-ain[i][jm] ) + ain[i][jm2] - ain[i][jp2] ) ); }
  else {return ( hy*.5*( ain[i][jp]-ain[i][jm] ) ) ; };
}

// ---------------------------------------------------------------------------
// 4th order x derivative:
inline double DFDX( AXY ain, int i, int j)
{ 
  int im, ip, im2, ip2;
  nextp(i,0,nx1,im,ip,im2,ip2); 
  if (b_4th) { return ( hy*r12*( 8.*( ain[ip][j]-ain[im][j] ) + ain[im2][j] - ain[ip2][j] ) ) ; }
  else {  return ( hy*.5*( ain[ip][j]-ain[im][j] ) ) ;  };
}

// ---------------------------------------------------------------------------
// 2D Laplace operator by finite differences: y = del^2 x
void laplace( AXY fi, AXY fo )
{ 
  int i, j, im, ip, jm, jp, im2, ip2, jm2, jp2;
  int i0 = 0, i1 = nx1, j0 = 0, j1 = ny1;  
  double aa, bb, cc, dd;
  
  if (b_4th) { aa = -60.; bb = +16.; cc = -1.; dd = hysq; } // 4th order
  else { aa = -4.; bb = +1.; cc = 0.; dd = hy2; }; // 2nd order

#pragma omp parallel for private(im,ip,im2,ip2,j,jm,jp,jm2,jp2) 
  for (i=i0; i<=i1; i++) {
    nextp(i,i0,i1,im,ip,im2,ip2); 
    
    for (j=0; j<=ny1; j++) {
      nextp(j,j0,j1,jm,jp,jm2,jp2);

      fo[i][j] = aa*fi[i][j] + bb*(fi[ip][j]+fi[i][jp]+fi[im][j]+fi[i][jm])
	+ cc*(fi[ip2][j]+fi[im2][j]+fi[i][jp2]+fi[i][jm2]);
      fo[i][j]*= dd;       
    }
  }
}

// ---------------------------------------------------------------------------
// FFT "Poisson equation" solver: solve in k-space for y in: del^2y = x
// with general arguments defined in cp (not only for the pure Poisson problem).
void poisson( AXY fi, AXY cp, AXY fo )
{
  int i,j,k;
  int nyp = ny/2+1;
  fftw_complex wwk[nx][nyp];
  double wwx[nx][ny];
  fftw_plan hinpsub, herpsub;
  
#ifdef _OPENMP
  omp_set_num_threads(npar);
  fftw_plan_with_nthreads(npar);
#endif

  fftw_import_wisdom_from_string(wisdom_sf); 
  hinpsub=fftw_plan_dft_r2c_2d(nx,ny,*wwx,*wwk,FFTW_ESTIMATE);
  herpsub=fftw_plan_dft_c2r_2d(nx,ny,*wwk,*wwx,FFTW_ESTIMATE);

  int i0 = 0, i1=nx1;

#pragma omp parallel for private(j)
  for (i=i0; i<=i1; ++i) for (j=0; j<=ny1; ++j) { wwx[i][j] = fi[i][j]; }
  fftw_execute(hinpsub);

#pragma omp parallel for private(j,k)
  for (i=i0; i<=i1; ++i) for (j=0; j<=nyp-1; ++j)  {
      for (k=0; k<=1; ++k) wwk[i][j][k] *= cp[i][j];
    }
  fftw_execute(herpsub);

#pragma omp parallel for private(j)
  for (i=i0; i<=i1; ++i) for (j=0; j<=ny1; ++j) fo[i][j] = wwx[i][j];

  fftw_destroy_plan(hinpsub);
  fftw_destroy_plan(herpsub);
}

// ---------------------------------------------------------------------------
inline void f_copy2darray( AXY arrinp, AXY arrout, int iend, int jend )
{
#pragma omp parallel for
  for (int i=0; i<=iend; ++i) for (int j=0; j<=jend; ++j) arrout[i][j] = arrinp[i][j]; 
}


// ---------------------------------------------------------------------------
double gamma0( double bkk)
// Evaluates the gyrofluid Gamma0 operator: Gamma0 = exp(-b)*I0(b):
// based on the modified Bessel function I0 with cosh approximation from:
// Olivares et al, Journal of Physics: Conference series 1043(2018) 012003.
// (Writing cosh in exp form in the Gamma operator avoids the range error!)
{
  return  0.5*(1.+exp(-2.*bkk))*(1.+0.24273*bkk*bkk)
    / ((1.+0.43023*bkk*bkk)*pow((1.+0.25*bkk*bkk),0.25));
}


// ---------------------------------------------------------------------------
inline double timer_start()
{ // measures first omp time
#ifdef _OPENMP 
  t_1 = omp_get_wtime(); 
#endif
  return t_1;
}

// ---------------------------------------------------------------------------
inline double timer_stop(double t_1)
{ // measures second omp time and returns time difference
#ifdef _OPENMP 
  t_2 = omp_get_wtime();
#endif
  return (t_2 - t_1);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

