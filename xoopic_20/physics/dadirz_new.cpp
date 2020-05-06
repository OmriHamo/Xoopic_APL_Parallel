/*
====================================================================

  This is a dynamic ADI solution of the equation
  Del( epsilon * Del(phi) ) = rho
  in z-r coordinates.
  This is an elliptic PDE, which is solved as parabolic PDE 
  (by adding an artificial time-dependent term) and looking at steady state.
  Neumann, Dirichlet, and symmetry boundary conditions are supported.
  Symmetry boundary condition is only applied when r0=0, and the
  neuman flag is specified there with zero slope on the electric field.

  dxdxu + dydyu =  s
  
  The function is based on a Peaceman Rachford Douglas 
  advance of the parabolic equation:  
  
  dtu = dxdxu + dydyu -s
  
  But with the time step dynamically chosen so as to speed convergence 
  to the dtu=0 steady state solution which is what we want.  It is the 
  user's responsiblity to put the initial guess of the solution stored 
  in u when the function is called.  u=0 or u=u(previous time step) 
  are possible choices.  
  
  The function sends back the finishing iteration number, the finishing
  normalized maximum error and location, and the failure code.  
  
  The failure codes are: 
    ifail=0, success actually
	 ifail=1, had to discard too many times
	 ifail=2, used up the iterations
	 
  Send in a negative tol_test if you want to freeze the time step to 
  the initial choice.  Send in a 0 adel_t to let program pick an 
  initial del_t.   

Revision/Programmer/Date
0.?	(Peterm, ??-??-94)	Conversion from C DADI.
0.9	(JohnV Peterm 08-09-95) Modify epsilon in set_coeff for dielectric
		blocks.

====================================================================
*/


#include <math.h>
#include "dadirz.h"
#include "grid.h"
#include <omp.h>
#include <threads_data.h>
#ifndef MAX
#define MAX(x, y)       (((x) > (y)) ? (x) : (y))
#endif

#ifndef MIN
#define MIN(x, y)       (((x) < (y)) ? (x) : (y))
#endif

#ifndef DBL_MIN
#define DBL_MIN         1E-200
#endif

#ifdef _MSC_VER
using std::cout;
using std::cerr;
using std::endl;
#endif
/**********************************************************************
  Single Peaceman Rachford Douglas pass with Direchlet 0 c boundary
  conditions for the equation: 
  
  dtu = dxdxu + dydyu -s, where s is constant in time.  
  
  The Crank-Nicolson finite difference approximation to the 
  above equation leads to the fractional step or 
  ADI equations: 
  
  u*(i,j)-(del_t/2dxdx)[u*(i+1,j)-2u*(i,j)+u*(i-1,j)] 
  = un(i,j)+(del_t/2dydy)[un(i,j+1)-2un(i,j)+un(i,j-1)] - (del_t/2)s(i,j) 
  
  un+1(i,j)-(del_t/2dydy)[un+1(i,j+1)-2un+1(i,j)+un+1(i,j-1)] 
  = u*(i,j)+(del_t/2dxdx)[u*(i+1,j)-2u*(i,j)+u*(i-1,j)] - (del_t/2)s(i,j) 
  
  **********************************************************************/



void Dadirz::adi(Scalar **uadi, Scalar **s, Scalar del_t)
{
  register int i, j;
  Scalar dth;

  dth = .5*del_t;

  /***************************************/
  /* Do z pass.  Set up variables for    */
  /* tridiagonal inversion along z.      */

  for(j=0;j<ng2;j++)
  {
    for(i=0; i<ng1; i++)
    {
      a_x1[i] = -dth*a_x1geom[i][j];
      b_x1[i] = 1 - dth*b_x1geom[i][j];
      c_x1[i] = -dth*c_x1geom[i][j];
		if(j==0)
		  b_x1[i] -= 4*dth*b_x2geom[i][j];
    }
  /*  handle the boundary conditions, neumann and dirichlet */


    for(i=0; i<ng1; i++)
		 if(b_x2geom[i][j]!=0.0)  // if this isn't dirichlet
			 r_x1[i] = uadi[i][j] +  dth*(-s[i][j]
													+((j>0)?(a_x2geom[i][j]*uadi[i][j-1]):0)
													+((j!=0)?b_x2geom[i][j]*uadi[i][j]:0) 
													+((j<nc2)?
													  ((j==0)?4:1)
													  *c_x2geom[i][j]*uadi[i][j+1]
													  :0));
		 else
			 r_x1[i] = uadi[i][j];  // metal sets potential here
									 

    /* Solve tridiagonal system. */
    tridag(a_x1, b_x1, c_x1, r_x1, v_x1, gam_x1, ng1);
    
    
    /* Copy solution into ustar. */
    for(i=0; i<ng1; i++) ustar[i][j] =v_x1[i];

  }
  
  /***************************************/
  /* Do r pass.  Set up variables for    */
  /* tridiagonal inversion along r.      */

  for(i=0;i<ng1;i++) {
	  for(j=0; j<ng2; j++) {
		  
		  a_x2[j] = -dth*a_x2geom[i][j];
		  b_x2[j] = 1 - dth*b_x2geom[i][j];
		  c_x2[j] = -dth*c_x2geom[i][j];
		  //      if(j==0&&N_Y0==2) {
		  if(j==0) {
			  b_x2[j] -= dth*(b_x1geom[i][j]+3*b_x2geom[i][j]);
			  c_x2[j] -= dth*(c_x2geom[i][j]*3);
		  }
	  }


	  for(j=0; j<ng2; j++)
		  if(b_x2geom[i][j]!=0)  //non-dirichlet
			  r_x2[j] = ustar[i][j] +dth*(-s[i][j] 
													+((i>0)?a_x1geom[i][j]*ustar[i-1][j]:0)
													+((j!=0)?b_x1geom[i][j]*ustar[i][j]:0) 
													+((i<nc1)?c_x1geom[i][j]*ustar[i+1][j]:0));
		  else
			  r_x2[j] = ustar[i][j];  // metal sets the potential here.

          
	  /* Solve tridiagonal system. */
	  tridag(a_x2, b_x2, c_x2, r_x2, v_x2, gam_x2, ng2);
    
	  /* Copy solution into ustar. */
	  for(j=0; j<ng2; j++) uadi[i][j] =v_x2[j];
  }
}











void Dadirz::adi_omp(Scalar **uadi, Scalar **s, Scalar del_t)
{
  register int i, j;
  Scalar dth;

  dth = .5*del_t;

	int N_threads_local = N_threads;
	int N_incr=0; //number of times to increment pointer at parallelization start 		

	int npt_x_sweep = (int)ng2 / N_threads_local ; // number of runs per thread		int n_extra_x_sweep = nc2 % N_threads_local; 
	int n_extra_x_sweep = ng2 % N_threads_local; // how many threads do +1 iteration
	int npt_y_sweep = (int)ng1 / N_threads_local ; // number of runs per thread		int n_extra_x_sweep = nc2 % N_threads_local; 
	int n_extra_y_sweep = ng1 % N_threads_local; // how many threads do +1 iteration

  /***************************************/
  /* Do z pass.  Set up variables for    */
  /* tridiagonal inversion along z.      */



#pragma omp parallel 
	{
	//int ID=omp_get_thread_num();
	int ii,jj;
	int  j_start=0,j_min,j_maxp1;
	// nc1 = J , nc2 = K
	//thread_ID[ID]+=1; //thread_ID[3] 
	//if (ID==0) 
	//  num_threads_used=omp_get_num_threads();

	//for(j=0; j<nc2; j++) // subtract 1, the last boundary is in fact the same as the 0'th
	//{	


#pragma omp critical
	{
		for (jj=0; jj<N_incr ; jj++)
			j_start++;

		N_incr++;
	}
		if (j_start<n_extra_x_sweep){		
			j_min = j_start*(npt_x_sweep+1);
			j_maxp1 = j_min + npt_x_sweep + 1;
		}
		else{
			j_min = j_start*npt_x_sweep + n_extra_x_sweep;
			j_maxp1 = j_min + npt_x_sweep;
		}	



	Scalar *a_local,*b_local,*c_local,*r_local,*vx_local,*gam_local ;

	a_local = new Scalar [nc1+1];
	b_local = new Scalar [nc1+1];
	c_local = new Scalar [nc1+1];	
	r_local = new Scalar [nc1+1];	
	vx_local = new Scalar [nc1+1];
	gam_local = new Scalar [nc1+1];

		

	for (jj=j_min; jj<j_maxp1; jj++){ 
  
		for(ii=0; ii<ng1; ii++)
    { // calculate tridiagonal matrix coefficients for ADI method
			a_local[ii] = -dth*a_x1geom[ii][jj]; // for equal spacing and not boundary: - a_ij = [eps(i-1,j-1)+eps(i,j-1)]/[2*dy^2]
			b_local[ii] = 1 - dth*b_x1geom[ii][jj];
			c_local[ii] = -dth*c_x1geom[ii][jj];
		if(jj==0)
		  b_local[ii] -= 4*dth*b_x2geom[ii][jj];
    }
  /*  handle the boundary conditions, neumann and dirichlet */


    for(ii=0; ii<ng1; ii++)
		 if(b_x2geom[ii][jj]!=0.0)  // if this isn't dirichlet
			 r_local[ii] = uadi[ii][jj] +  dth*(-s[ii][jj]
													+((jj>0)?(a_x2geom[ii][jj]*uadi[ii][jj-1]):0)
													+((jj!=0)?b_x2geom[ii][jj]*uadi[ii][jj]:0) 
													+((jj<nc2)?
													  ((jj==0)?4:1)
													  *c_x2geom[ii][jj]*uadi[ii][jj+1]
													  :0));
		 else
			 r_local[ii] = uadi[ii][jj];  // metal sets potential here


				/* Solve tridiagonal system. */
		tridag(a_local, b_local, c_local, r_local, vx_local, gam_local, ng1);
		/* Copy solution into ustar. */
		for(ii=0; ii<ng1; ii++) ustar[ii][jj] =vx_local[ii];  

	} // end for loop ii


	delete [] a_local;
	delete [] b_local;
	delete [] c_local;		
	delete [] r_local;
	delete [] vx_local;		
	delete [] gam_local;

}  // end parallel region



  /***************************************/
  /* Do r pass.  Set up variables for    */
  /* tridiagonal inversion along r.      */

N_incr=0;


#pragma omp parallel 
{


	int ii,jj;
	int  i_start=0;
	int  i_min,i_maxp1;	


#pragma omp critical
	{
	for (ii=0; ii<N_incr ; ii++)
		i_start++;

	N_incr++;
	}

if (i_start<n_extra_y_sweep){		
	i_min = i_start*(npt_y_sweep+1);
	i_maxp1 = i_min + npt_y_sweep + 1;
}
else{
	i_min = i_start*npt_y_sweep + n_extra_y_sweep;
	i_maxp1 = i_min + npt_y_sweep;
}	

Scalar *a_local,*b_local,*c_local,*r1_local,*r2_local,*vx2_local,*gam_local ;

a_local = new Scalar [nc2+1];
b_local = new Scalar [nc2+1];
c_local = new Scalar [nc2+1];	
r1_local = new Scalar [nc2+1];	
r2_local = new Scalar [nc2+1];	
vx2_local = new Scalar [nc2+1];	  
gam_local = new Scalar [nc2+1];


for (ii=i_min; ii<i_maxp1; ii++){   
	  for(jj=0; jj<ng2; jj++) {
		  
		a_local[jj] = -dth*a_x2geom[ii][jj]; // 
		b_local[jj] = 1 - dth*b_x2geom[ii][jj];
		c_local[jj] = -dth*c_x2geom[ii][jj];
		  //      if(jj==0&&N_Y0==2) {
		  if(jj==0) {
			  b_local[jj] -= dth*(b_x1geom[ii][jj]+3*b_x2geom[ii][jj]);
			  c_local[jj] -= dth*(c_x2geom[ii][jj]*3);
		  }
	  }


	  for(jj=0; jj<ng2; jj++)
		  if(b_x2geom[ii][jj]!=0)  //non-dirichlet
			  r2_local[jj] = ustar[ii][jj] +dth*(-s[ii][jj] 
													+((ii>0)?a_x1geom[ii][jj]*ustar[ii-1][jj]:0)
													+((jj!=0)?b_x1geom[ii][jj]*ustar[ii][jj]:0) 
													+((ii<nc1)?c_x1geom[ii][jj]*ustar[ii+1][jj]:0));
		  else
			  r2_local[jj] = ustar[ii][jj];  // metal sets the potential here.

          
	  /* Solve tridiagonal system. */
	  tridag(a_local, b_local, c_local, r2_local, vx2_local, gam_local, ng2);
    
	  /* Copy solution into ustar. */
	  for(jj=0; jj<ng2; j++) uadi[ii][jj] =vx2_local[jj];
  }


delete [] a_local;
delete [] b_local;
delete [] c_local;		
delete [] r1_local;
delete [] r2_local;
delete [] vx2_local;		
delete [] gam_local;

} // end parallel region

}










void Dadirz::init_solve(Grid *grid,Scalar **epsi)
{
  register int i, j;
  //set up freespace default coefficients everywhere.
  //outside boundaries must be set up specially, 
  //but at this point we don't know these yet.
  for(i=0;i<ng1;i++)
	  for(j=0;j<ng2;j++)
		 {
			Boundary *B = grid->GetNodeBoundary()[i][j];
			BCTypes type;
			if(B!=NULL) type = B->getBCType();
			else type = FREESPACE;
			set_coefficient(i,j,type,grid);
		 }

  /************************************************************/
  /* Determine initial step size such that 1=del_t/(dx1*dx1+dx2*dx2).  
     Chosen because this looks like stability condition for 
     an explicit scheme for the above parabolic equation.  
     Also calculate the usual constants.  */
#ifndef sqr
#define sqr(x) ((x)*(x))
#endif
  
  Vector2 **X = grid->getX();
  del_t0 = 0.1*(sqr(X[1][0].e1()-X[0][0].e1()) +sqr(X[0][1].e2()-X[0][0].e2()))/epsi[0][0];

  for(i=1; i<ng1; i++)
    for(j=1;j<ng2;j++)
      del_t0 = MIN(del_t0, 0.1*(sqr(X[i][j].e1()-X[i-1][j].e1()) 
			     +sqr(X[i][j].e2()-X[i][j-1].e2()))/epsi[0][0]);
  
}

/*  set up the coefficient at i,j, which can be one of several types
  epsilon is assumed to be cell-centered so some weighted averaging is
  needed to get it at half-cells. */

void Dadirz::set_coefficient(int i, int j, BCTypes type, Grid *grid) 
{  Vector2 **X=grid->getX();
	int I=grid->getJ();
	int J=grid->getK();
	// dz-
	Scalar dxia = X[i][j].e1()  -X[MAX(0,i-1)][j].e1();
	// dz+
	Scalar dxib = X[MIN(I,i+1)][j].e1()-X[i][j].e1();
	// dz
	Scalar dxave= 0.5*(dxia + dxib);

	Scalar ra = X[i][MAX(0,j-1)].e2();
	Scalar r  = X[i][j].e2();
	Scalar rb = X[i][MIN(J,j+1)].e2();
	Scalar ra2= 0.5*(ra + r);
	Scalar rb2= 0.5*(r +rb);
	Scalar dr2 = rb2*rb2-ra2*ra2;
	Scalar drb = rb -  r;
	Scalar dra = r  - ra;
 

	Scalar da1a,da1b,da2aa,da2ab,da2ba,da2bb;  // Areas of surfaces, for epsilons
	
	Scalar e1a,e1b,e2a,e2b;  /* epsilon (i-1/2, j), (i+1/2,j),(i,j-1/2),(i,j+1/2) */
	
	// the area weightings to use for the epsilons
	da1a = sqr(r) - sqr(ra2);
	da1b = sqr(rb2)-sqr(r);
	da2aa= dxia;
	da2ab= dxib;
	da2ba= dxia;
	da2bb= dxib;
	
	// the epsilons to use for each coefficient below
//	e1a = (epsi[i][j]*da1a   + epsi[i][MIN(J,j+1)]*da1b    ) / ( da1a  + da1b);
//	e1b = (epsi[MIN(I,i+1)][j]*da1a + epsi[MIN(I,i+1)][MIN(j+1,J)]*da1b  ) / ( da1a  + da1b);
//	e2a = (epsi[i][j]*da2aa  + epsi[MIN(I,i+1)][j] * da2ab  ) / ( da2aa + da2ab);
//	e2b = (epsi[i][MIN(J,j+1)]*da2ba  + epsi[MIN(I,i+1)][j] * da2bb  ) / ( da2aa + da2ab);
	int im = MAX(i-1,0);
	int jm = MAX(j-1,0);
	e1a = (epsi[im][jm]*da1a   + epsi[im][MIN(j,J-1)]*da1b ) / ( da1a  + da1b);
	e1b = (epsi[MIN(i,I-1)][jm]*da1a + epsi[MIN(i,I-1)][MIN(j,J-1)]*da1b  ) / ( da1a  + da1b);
	e2a = (epsi[im][jm]*da2aa  + epsi[MIN(i,I-1)][jm] * da2ab  ) / ( da2aa + da2ab);
	e2b = (epsi[im][MIN(j,J-1)]*da2ba  + epsi[MIN(i,I-1)][MIN(j,J-1)] * da2bb  ) / ( da2aa + da2ab);

	// If the dielectric boundary isn't at an edge, treat it like
	// freespace
  if(type==DIELECTRIC_BOUNDARY && (i!=0||i!=I||j!=J||j!=0)) type = FREESPACE;

  
  switch(type) 
	 {
	 case FREESPACE:
		{  
			if(b_x2geom[i][j]==0) break;  // if dirichlet, don't override.
			//Z finite difference
		  if(i) a_x1geom[i][j] = e1a/(dxia*dxave);
		  if(i!=I) c_x1geom[i][j] = e1b/(dxib*dxave);
		  b_x1geom[i][j] = -( a_x1geom[i][j] + c_x1geom[i][j] );

		  //R finite difference

		  if(j) a_x2geom[i][j] = e2a*2*ra2/(dr2 * dra);
		  if(j!=J) c_x2geom[i][j] = e2b*2*rb2/(dr2 * drb);
		  b_x2geom[i][j] = - ( a_x2geom[i][j] + c_x2geom[i][j] );

		  break;
		 }
	 case CONDUCTING_BOUNDARY:
		{
		  a_x1geom[i][j]=0.0;
		  b_x1geom[i][j]=0.0;
		  c_x1geom[i][j]=0.0;
		  a_x2geom[i][j]=0.0;
		  b_x2geom[i][j]=0.0;
		  c_x2geom[i][j]=0.0;
		  break;
		}
	 case DIELECTRIC_BOUNDARY:
		{
			if(i==0&&b_x2geom[i][j]!=0) {  //Don't override a dirichlet with a Neumann
					Scalar xa = X[i][j].e1();
					Scalar xave = 0.5*(2*xa + dxib);

					a_x1geom[i][j] = 0.0;
					c_x1geom[i][j] = e1b/(xave-xa)/dxib;
					b_x1geom[i][j] = -( a_x1geom[i][j] + c_x1geom[i][j] );
					if(j==0) {
						Scalar dria = X[i][1].e2() - X[i][0].e2();
						a_x2geom[i][0] = 0.0;
						b_x2geom[i][0] = -e2b/dria/dria;
						c_x2geom[i][0] = -b_x2geom[i][0];
					}
					if(j==J) { //neumann neumann corner
						c_x2geom[i][j] = 0.0;
						a_x2geom[i][j] = e2a*2*ra2/(dr2*dra);
						b_x2geom[i][j] = -a_x2geom[i][j];
					}

			 }
			if(i==I&&b_x2geom[i][j]!=0) {
					Scalar xb = X[i][j].e1();
					Scalar xave = 0.5*(X[i][j].e1() + X[i-1][j].e1());
					
					a_x1geom[i][j] = e1a/(xb-xave)/dxia;
					c_x1geom[i][j] = 0.0;
					b_x1geom[i][j] = -( a_x1geom[i][j] + c_x1geom[i][j] );
				if(j==0&&b_x2geom[i][j]!=0) { //neumann neumann corner
				  Scalar dria = X[i][j+1].e2() - X[i][j].e2();
				  a_x2geom[i][j] = 0.0;
				  c_x2geom[i][j] = e2b/dria/dria;
				  b_x2geom[i][j] = -c_x2geom[i][j];
				}
				if(j==J&&b_x2geom[i][j]!=0) { //neumann neumann corner
				  c_x2geom[i][j] = 0.0;
				  a_x2geom[i][j] = e2a*2*ra2/(dr2*dra);
				  b_x2geom[i][j] = -a_x2geom[i][j];
				}
			}
			
			if(j==J && i!=0&&b_x2geom[i][j]!=0 ) { 
				Scalar dria = X[i][j].e2()-X[i][j-1].e2();

				a_x2geom[i][nc2] = 2 * e2a * ra2/ ( dria * ( r*r - ra2*ra2));
				b_x2geom[i][nc2] = -a_x2geom[i][nc2];
				c_x2geom[i][nc2] = 0.0;
			}
			break;
		}
	 case CYLINDRICAL_AXIS:
		{ if(i==0||i==I) break;  // don't want cyl_axis in corners
		  Scalar dria = X[i][1].e2() - X[i][0].e2();
		  a_x2geom[i][0] = 0.0;
		  b_x2geom[i][0] = -e2b/dria/dria;
		  c_x2geom[i][0] = -b_x2geom[i][0];

		  a_x1geom[i][j] = e1a/(dxia*dxave);
		  c_x1geom[i][j] = e1b/(dxib*dxave);
		  b_x1geom[i][j] = -( a_x1geom[i][j] + c_x1geom[i][j] );

		  break;
		}
	 default: 
		{
		  cout << "DADI doesn't know about this boundary condition type!\n;" << endl;
		  cout << "(was either PERIODIC_BOUNDARY or SPATIAL_REGION_BOUNDARY)\n" << endl;
		  //e_xit(1);
		}
	 }
 }

BCTypes Dadirz::get_coefficient(int j, int k) {
  return (b_x2geom[j][k]==0)? CONDUCTING_BOUNDARY:FREESPACE;
}

/**********************************************************************/





int Dadirz::solve(Scalar **u_in, Scalar **s, int itermax,Scalar tol_test)
{
	// u_in is del_phi to be calculated 
	register int i, j;
	int iter, ndiscard,N_incr=0;
	static Scalar del_t=0.0;
	Scalar  ratio=0 , del_td=0 , tptop=0, tpbot=0; // del_td - double time step for DADI
	Scalar rnorm=0, rsum=0, res=0;

	rnorm=rsum=0.0;




#pragma omp parallel 
	{

	Scalar rnorm_local=0, rsum_local=0, errchk_local, dxdxutrm_local, dydyutrm_local; // rnorm - normalized residual(?)	
	//int ID=omp_get_thread_num();
	int ii,jj;
	int  i_start=0, i_min,i_maxp1 ;
	int npt = (int)ng1 / N_threads ; // number of runs per thread		int n_extra_x_sweep = nc2 % N_threads_local; 
	int n_extra = ng1 % N_threads; // how many threads do +1 iteration	
	// nc1 = J , nc2 = K
	//thread_ID[ID]+=1; //thread_ID[3] 
	//if (ID==0) 
	//  num_threads_used=omp_get_num_threads();

#pragma omp critical
	{
		for (ii=0; ii<N_incr ; ii++)
			i_start++;

		N_incr++;
	}

	if (i_start<n_extra){		
		i_min = i_start*(npt+1);
		i_maxp1 = i_min + npt + 1;
	}
	else{
		i_min = i_start*npt + n_extra;
		i_maxp1 = i_min + npt;
	}	




	for(ii=i_min; ii< i_maxp1; ii++)
		for(jj=0; jj< ng2; jj++) {


			/* Residual normalization.  */
			// dirichlet conditions don't add to rnorm
			rnorm_local +=( (b_x2geom[ii][jj]==0)? 0:s[ii][jj]*s[ii][jj]) ; // for dirichlet b(i,j)=1


			/*  copy u_in to u for working purposes.  */
			u[ii][jj]=(Scalar)u_in[ii][jj];


			//calculate an initial estimate of the residual
			/* Use the residual as the absolute error and if it is bigger 
			 than the stored maximum absolute error, record new maximum 
			 absolute error and location.  */

			if(ii>0&& jj>0 && ii<nc1 && jj<nc2) { // not on dirichlet boundary, but yes on periodic BC
				dxdxutrm_local=a_x1geom[ii][jj]*u_in[ii-1][jj] +b_x1geom[ii][jj]*u_in[ii][jj] +
					c_x1geom[ii][jj]*u_in[ii+1][jj];
				dydyutrm_local=a_x2geom[ii][jj]*u_in[ii][jj-1] +b_x2geom[ii][jj]*u_in[ii][jj] +
					c_x2geom[ii][jj]*u_in[ii][jj+1];
			}


			/* only include points which are not in structures. */
			errchk_local = ((b_x2geom[ii][jj]==0)? 0:dxdxutrm_local +dydyutrm_local -s[ii][jj]);

			/* Residual sums. */
			rsum_local += errchk_local*errchk_local;
		}


#pragma omp critical
	{
		rsum+=rsum_local;
		rnorm+=rnorm_local;
	}


} // end parallel

// If rnorm is zero, we must deal with it...
if (rnorm == 0.0) {

	// check dirichlet conditions
	for(i=0;i<ng1;i++) for(j=0;j<ng2;j++) {
		// check left
		rnorm += sqr(((i>0&&b_x2geom[i-1][j]!=0)?c_x1geom[i-1][j]*u[i][j]:0));
		// check right
		rnorm += sqr(((i<nc1&&b_x2geom[i+1][j]!=0)?a_x1geom[i+1][j]*u[i][j]:0));
		// check up
		rnorm += sqr(((j>0&&b_x2geom[i][j-1]!=0)?c_x2geom[i][j-1]*u[i][j]:0));
		// check right
		rnorm += sqr(((j<nc2&&b_x2geom[i][j+1]!=0)?c_x2geom[i][j+1]*u[i][j]:0));

	}
	//rnorm = sqrt(rnorm);
	if(rnorm == 0) { //still zero, we don't need to iterate
		for(i=0;i<ng1;i++) for(j=0;j<ng2;j++) u_in[i][j]=0;
		return 0;
	}
}
rnorm=sqrt(rnorm);
res = sqrt(rsum)/rnorm;
#ifdef DADI_DEBUG
printf("dadi: res = %g\n",res);
printf("dadi: rnorm= %lg\n", rnorm);
#endif

if(res<tol_test) return 0;  // no need to iterate

/*************************************************/
if(del_t<1.0e-10*del_t0) del_t=del_t0; else del_t/=4;
del_td= 2.0*del_t;
ndiscard=0;

/********************/
/* Begin iteration. */

for(iter=0; iter<itermax; iter++)
{
	/*************************************************/
	/* Copy u into the work array and storage array. */

	for(i=0; i< ng1; i++)
		for(j=0; j< ng2; j++) uwork[i][j] = ustor[i][j] = u[i][j];

	/************************************/
	/* Two advances of u via ADI at del_t. */
	/* orig code:
	 adi(u, s, del_t);
	 adi(u, s, del_t);

	 /////////////////////////////////////////////////////////
		 // One advance of uwork via ADI at 2*del_t. //

		 adi(uwork, s, del_td);
		 end orig code */



	/* parallel code*/
	adi_omp(u, s, del_t); // u is both input & outpout
	adi_omp(u, s, del_t);

	// One advance of uwork via ADI at 2*del_t.

	adi_omp(uwork, s, del_td);

	/*end parallel code*/







	/*******************************************************/
	/* Calculate test parameter and normalized error.
	 For Dirichlet BCs, no need to worry about boundary
	 points since u,uwork, and ustor should be the same. */

	tptop= tpbot= rsum = 0.0;
	N_incr=0;

#pragma omp parallel 
	{

		Scalar tptop_local=0,tpbot_local=0 , rsum_local=0, errchk_local, dxdxutrm_local, dydyutrm_local; // rnorm - normalized residual(?)	
		//int ID=omp_get_thread_num();
		int ii,jj;
		int  i_start=1; // not includin boundary
		// nc1 = J , nc2 = K
		//thread_ID[ID]+=1; //thread_ID[3] 
		//if (ID==0) 
		//  num_threads_used=omp_get_num_threads();

#pragma omp critical
		{
			for (ii=0; ii<N_incr ; ii++)
				i_start++;

			N_incr++;
		}


		for(ii=i_start; ii< nc1; ii=ii+N_threads) // not including boundary
			for(jj=1; jj< nc2; jj++)
		{
				/* Test paramter sums. */
				tptop_local += (u[ii][jj]-uwork[ii][jj])*(u[ii][jj]-uwork[ii][jj]);
				tpbot_local += (u[ii][jj]-ustor[ii][jj])*(u[ii][jj]-ustor[ii][jj]);

				/* Residual terms. */

				dxdxutrm_local=a_x1geom[ii][jj]*u[ii-1][jj] +b_x1geom[ii][jj]*u[ii][jj] +c_x1geom[ii][jj]*u[ii+1][jj];
				dydyutrm_local=a_x2geom[ii][jj]*u[ii][jj-1] +b_x2geom[ii][jj]*u[ii][jj] +c_x2geom[ii][jj]*u[ii][jj+1];

				/* Use the residual as the absolute error and if it is bigger 
				 than the stored maximum absolute error, record new maximum 
				 absolute error and location.  */
				/* only include points which are not in structures. */
				errchk_local = ((b_x2geom[ii][jj]==0)? 0:dxdxutrm_local +dydyutrm_local -s[ii][jj]);

				/* Residual sums. */
				rsum_local += errchk_local*errchk_local;
			}

#pragma omp critical
		{
			rsum+=rsum_local;
			tptop+=tptop_local;
			tpbot+=tpbot_local;
		}


	}//end parallel


	/* Calculate normalized residual. */
	res = sqrt(rsum)/rnorm;
#ifdef DADI_DEBUG
	printf("dadi: iter= %d, res = %lg del_t=%lf\n", iter, res,del_t);
#endif
		/* If the residual is less than the tolerance, SUCCESS! */
		if ((res < tol_test) && (iter))
	{
#ifdef DADI_DEBUG
			printf("dadi: iter=%d\n", iter);
#endif
				for(i=0;i<ng1;i++)
				for(j=0;j<ng2;j++)
				u_in[i][j]=(Scalar)u[i][j];

			return(0);
		}

	/* Determine ratio used to find the time step change.  If tpbot 
	 is zero but tptop is finite, consider this a case of a large 
	 ratio and act accordingly.  DWH does about the same thing 
	 except he does NOT discard the solution if tpbot=0. */

	if (tpbot > 0.0) ratio=tptop/tpbot;
	if (tpbot== 0.0) ratio=1.0;

	/* Get next time step. */
	if (ratio < 0.02) del_t *= 8.0;
	if (ratio >=0.02 && ratio < 0.05) del_t *= 4.0;
	if (ratio >=0.05 && ratio < 0.10) del_t *= 2.0;
	if (ratio >=0.10 && ratio < 0.30) del_t *= 0.80;
	if (ratio >=0.30 && ratio < 0.40) del_t *= 0.50;
	if (ratio >=0.40 && ratio < 0.60) del_t *= 0.25;
	for(i=0;i<ng1;i++)
		for(j=0;j<ng2;j++)
		u_in[i][j]=(Scalar)u[i][j];


	/* Ratio is too large. */
	if (ratio >=0.60)
	{ 
		ndiscard++;
#ifdef DADI_DEBUG
		printf("ndiscard= %d, iter=%d step=%lf\n", ndiscard, iter,del_t); 
#endif      
			/* Check if too many discards. */
			if (ndiscard > 20)
		{
				for(i=0;i<ng1;i++)
					for(j=0;j<ng2;j++)
					u_in[i][j]=(Scalar)u[i][j];
				del_t = del_t0;
				printf("Poisson solve FAILURE: dadi: iter= %d, ndiscard>10\n", iter);
				return 1;
			}	
		/* Discard by replacing u with what we started with. */
		for(i=0; i< ng1; i++)
			for(j=0; j< ng2; j++) u[i][j] = ustor[i][j];

		/* Reduce del_t. */
		del_t /= 8.0;
		//del_t = del_t0;
	}
	del_td=2*del_t;
}
/* Fail if used up the maximum iterations. */
printf("Poisson solve FAILURE:dadi:  iter>= %d\n",itermax);
for(i=0;i<ng1;i++)
for(j=0;j<ng2;j++)
u_in[i][j]=(Scalar)u[i][j];

return(2);
}





Dadirz::~Dadirz()
{

}

