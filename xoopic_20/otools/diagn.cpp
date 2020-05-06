/*
====================================================================

diagn.cpp

This file contains all the arrays and functions necessary for maintaining 
and initializing the diagnostics' data.


1.1  (KC 01-6-97) Not calculation the KE of the particles here anymore,
     moved it to the push.  Add KE by species. 
2.01 (Bruhwiler 10-08-99) initialization of synchRadiationFlag added
2.02 (Bruhwiler 10-05-2000) added history plots for RMS beam quantities

====================================================================
*/

#ifdef UNIX
#include <config.h>
#endif


#include <unistd.h>
#include <cmath>

#include <omp.h>

#include <stdlib.h>
#include <stdio.h>

#include "fields.h"
#include "grid.h"

#include "ovector.h"
#include <sptlrgn.h>
#include <ptclgrp.h>
#include "globalvars.h"
#include "oopiclist.h"
#include "diagn.h"
#include "gpdist.h"
#include "history.h"
#include "newdiag.h"


// hdf5 includes
#include "h5diagg.h"



// NPARTICLES specifies that, when the number of particles
//  of a given species exceeds this number, then XGrafix
//  plots only every 2nd particle.  It does not affect the
//  physics simulation.
#define NPARTICLES 800000
//#define HISTMAX 4096

Scalar AllDiagnostics = 1;
Scalar EnerPlots = 1;
Scalar PhaseSpacePlots = 1;

oopicList <Diagnostics> *theDiagnostics;

// omri write to file init
#define iter_to_avg 10
#define iter_to_save 100
char* theInputFile;
string out_dir;


extern int N_threads;
#ifdef MPI_VERSION
#include <mpi.h>
extern MPI_Comm XOOPIC_COMM;
extern int MPI_RANK;
#define SRB_ANNOUNCE_TAG 100
#define SRB_LINK_TAG 200
#define N_DIAGNOSTIC_TAG 300


typedef struct {
  int index;
  int linkedP;
  char *name;	
} SRBdat;
#endif


//  These are the diagnostics arrays:  they're global in scope but used
//  only here and in Maintian_Diagnostics and in xgrafix.

Diagnostics::Diagnostics(SpatialRegion *thisSpace) {
  number_of_species = 1;
  Bz0=0.0;
  divderrormag = 0;
  Show_loaded_densityFlag = FALSE;
  theSpace=thisSpace;
  RegionName = theSpace->getName();
  const int jm = theSpace->getJ();
  Jm = jm;
  const int km = theSpace->getK();

  number_of_species=theSpace->get_nSpecies();
  int j,k,isp;
  B = theSpace->getBNode();
  E = theSpace->getENode();
#ifdef DEBUG_FIELDS
  intBdS = theSpace->getIntBdS();
  intEdl = theSpace->getIntEdl();
#endif
  blist = theSpace->getBoundaryList();
  BoundaryIhist = new oopicList<Ihistdiag>;
  BoundaryPFhist = new oopicList<PFhistdiag>;
  BDynamic = theSpace->getBNodeDynamic();
  I = theSpace->getI();
  rho = theSpace->getRho();

  /** Set the local pointer to oopicList<NGD>* to point
   * to the corrensponding list on the grid.
   * dad: 01/24/2001
   */
  ptrNGDList = theSpace->getPtrNGDList();

  boltzmannFlag = theSpace->getBoltzmannFlag();
  electrostaticFlag = theSpace->getElectrostaticFlag();
  synchRadiationFlag = theSpace->getSynchRadiationFlag();
  if (boltzmannFlag){
		totalRho = theSpace->getTotalRho();
	}
  loaded_density = theSpace->getloaddensity();
  Show_loaded_densityFlag = theSpace->getShowInitialDensityFlag();
  phi = theSpace->getPhi();
  divderror = theSpace->getDivDerror();
  rho_species = theSpace->get_rho_species();
  theSpecies= new SpeciesDiag[number_of_species];
  CellVolumes=theSpace->get_halfCellVolumes();

  // get memory for the grid arrays
  x1_array = new Scalar[jm+2];
  x2_array = new Scalar[km+2];

  //get memory for the diagnostics
  S_array= new Vector3 *[jm+2];
  Ue= new Scalar * [jm+2];
  Ub= new Scalar * [jm+2];
		
#ifdef DEBUG_PHI
		phi_err = new Scalar *[jm+2];
#endif

  // Set up the user-defined diagnostics
  oopicListIter<Diag> DiagIter(*thisSpace->getDiagList());
  for(DiagIter.restart();!DiagIter.Done();DiagIter++) 
	 DiagIter.current()->setDiagnostics(this);
  
  for(j=0;j<=jm;j++)
	 {
		S_array[j]= new Vector3 [km+2];
		for(k=0;k<=km;k++)
			S_array[j][k] = Vector3(0,0,0);

		Ue[j]=new Scalar[km+2];
		memset(Ue[j],0,(km+2)*sizeof(Scalar));

		Ub[j]=new Scalar[km+2];
		memset(Ub[j],0,(km+2)*sizeof(Scalar));
		
#ifdef DEBUG_PHI
		phi_err[j] = new Scalar[km+2];
#endif

	 }

  //  Get memory for the arrays that contain particle positions for xgrafix

  for(isp=0;isp<number_of_species;isp++) {
	 // search for the name corresponding to this species
	 oopicListIter<Species> siter(*theSpace->getSpeciesList());
	 for(siter.restart();!siter.Done();siter++) {
		if(isp==siter.current()->getID()) {
		  theSpecies[isp].name=siter.current()->get_name_ptr();
		  theSpecies[isp].spec = siter.current();
		  theSpecies[isp].ID = siter.current()->getID(); 
		}
	 }

	 theSpecies[isp].memory_size=NPARTICLES;
	 theSpecies[isp].nparticles=0;
	 theSpecies[isp].nparticles_plot = 0;
	 theSpecies[isp].x1positions= new Scalar[NPARTICLES];
	 theSpecies[isp].x2positions= new Scalar[NPARTICLES];
	 theSpecies[isp].x1velocities=new Scalar[NPARTICLES];
	 theSpecies[isp].x2velocities=new Scalar[NPARTICLES];
	 theSpecies[isp].x3velocities=new Scalar[NPARTICLES];
	 memset(theSpecies[isp].x1positions,0,NPARTICLES*sizeof(Scalar));
	 memset(theSpecies[isp].x2positions,0,NPARTICLES*sizeof(Scalar));
	 memset(theSpecies[isp].x1velocities,0,NPARTICLES*sizeof(Scalar));
	 memset(theSpecies[isp].x2velocities,0,NPARTICLES*sizeof(Scalar));
	 memset(theSpecies[isp].x3velocities,0,NPARTICLES*sizeof(Scalar));
  }


  number = new Scalar_History* [number_of_species];
  ke_species = new Scalar_History* [number_of_species];
  ngroups = new Scalar_History* [number_of_species];
  total_density = new Scalar_History* [number_of_species];
  Ave_KE = new Scalar_History* [number_of_species];
  Ave_KE_gridprob = new Scalar [number_of_species];	//IgalK
  number_gridprob = new Scalar [number_of_species];	//IgalK
  total_density_gridprob = new Scalar [number_of_species];		//IgalK
  Ir_gridprob = new Scalar [64];	


	 

  HISTMAX = theSpace->get_histmax(); // retrieve time array length

  for(isp=0;isp<number_of_species;isp++) {
    number[isp]= new Scalar_History(HISTMAX,2);
    ke_species[isp] = new Scalar_History(HISTMAX,2);
    ngroups[isp] = new Scalar_History(HISTMAX,2);
    total_density[isp] = new Scalar_History(HISTMAX,2); 
    Ave_KE[isp] = new Scalar_History(HISTMAX,2);
  }
  
  /* **************************************************************
   * Need to collect scalar history of RMS beam quantities
   * for beam physics applications:  velocities, beam sizes,
   * energy spread, emittance, etc.
   * (Bruhwiler, revisions started on 10/04/2000)
   * (Bruhwiler/Dimitrov,  started on 10/09/2000)
   ****************************************************************/
  int sh;
  aveBeamSize = new Scalar_History* [2];
  rmsBeamSize = new Scalar_History* [2];
  for ( sh = 0; sh < 2; sh++ ) {
    aveBeamSize[sh]  = new Scalar_History(HISTMAX,2); 
    rmsBeamSize[sh]  = new Scalar_History(HISTMAX,2); 
  }

  // get memory for the velocity (three components in this case)
  aveVelocity = new Scalar_History* [3];
  rmsVelocity = new Scalar_History* [3];
  for ( sh = 0; sh < 3; sh++ ) {
    aveVelocity[sh] = new Scalar_History(HISTMAX,2);  
    rmsVelocity[sh] = new Scalar_History(HISTMAX,2);
  } 
  // get memory for the emittance (two components in this case)
  rmsEmittance = new Scalar_History* [2];
  for ( sh = 0; sh < 2; sh++ ) {
    rmsEmittance[sh] = new Scalar_History(HISTMAX,2);
  }  
  // get memory for the aveEnergy_eV and rmsEnergy_eV
  ///aveEnergy_eV = new Scalar_History(HISTMAX,2);
  ///rmsEnergy_eV = new Scalar_History(HISTMAX,2);

  //  Set up the grid arrays
  for(j=0;j<=jm;j++)
	 x1_array[j]=theSpace->getMKS(j,0).e1();
  for(k=0;k<=km;k++)
	 x2_array[k]=theSpace->getMKS(0,k).e2();

  //  get the size of the grid
  x1max = theSpace->getMKS(jm,km).e1();
  x2max = theSpace->getMKS(jm,km).e2();
  x1min = theSpace->getMKS(0,0).e1();
  x2min = theSpace->getMKS(0,0).e2();

  /**
   * allocate memory and initialize the 
   * x1_arrayNGD and x2_arrayNGD
   */
  jmNGD = jm;
  kmNGD = km;
  Vector2 x;  
  x1_arrayNGD = new Scalar[jm+1];
  x2_arrayNGD = new Scalar[km+1];
  for(j=0;j<jm;j++) {
    x.set_e1( static_cast<Scalar>(j) + 0.5 );
    x.set_e2( 0.5 );
    x1_arrayNGD[j] = (theSpace->getMKS(x)).e1();
  }
  for(k=0;k<km;k++) {
    x.set_e1( 0.5 );
    x.set_e2( static_cast<Scalar>(k) + 0.5 );
    x2_arrayNGD[k] = (theSpace->getMKS(x)).e2();
  }

  //set up diagnsotic history arrays for boundaries which
  //have particle diagnostics on
  oopicListIter<Boundary> nextb(*blist);
  for(nextb.restart();!nextb.Done();nextb++) { //if it's taking diagnostics of particles
     if(nextb.current()->get_particle_diag()) {
        Ihistdiag *Ihistory= new Ihistdiag(number_of_species);
        int len = nextb.current()->get_Ihist_len();
        int avg = nextb.current()->get_Ihist_avg();
        Ihistory->Ihist = new Scalar_Ave_History(len, avg);
        Ihistory->p_dist = nextb.current()->get_particle_diag();
        strcpy(Ihistory->name,nextb.current()->getBoundaryName().c_str());
        for (int i=0; i<number_of_species; i++) {
           Ihistory->Ihist_sp[i] = new Scalar_Ave_History(len,avg);
        }
        BoundaryIhist->add(Ihistory);
     }
  }
  for(nextb.restart();!nextb.Done();nextb++){ 
     if(nextb.current()->getPF()) {
        PFhistdiag *PFhistory= new PFhistdiag;
        PFhistory->PFhist = new Scalar_History(HISTMAX,2);
        PFhistory->PFlocal = new Scalar_Local_History(HISTMAX,MAX(1,HISTMAX/16));
        PFhistory->p_flux = nextb.current()->getPF();
        strcpy(PFhistory->name,nextb.current()->getBoundaryName().c_str());
        BoundaryPFhist->add(PFhistory);
     }
  }

// find anode and cathode indices
	int count_b=0; 
	oopicListIter<Ihistdiag> nextdiag(*BoundaryIhist);
		for(nextdiag.restart();!nextdiag.Done();nextdiag++) { 
			if (nextdiag.current()->p_dist->j1==0 && nextdiag.current()->p_dist->j2==0)
				theSpace->set_anode_index (count_b); //nextdiag.current()->p_dist->k2
			if (nextdiag.current()->p_dist->j1==theSpace->getJ() && nextdiag.current()->p_dist->j2==theSpace->getJ())
				theSpace->set_cathode_index (count_b); //nextdiag.current()->p_dist->k2

			count_b++;
		}



//  omri write to file - create diagnostics arrays

	  t_avg = simulation_time + (iter_to_avg + 1)/2.0*dt; // initial value of average time. does not need to actually averge (average can be calculated analytically)
	  KE_e_vec_save = new Scalar[iter_to_save];
	  KE_i_vec_save = new Scalar[iter_to_save];
	  KE_n_vec_save = new Scalar[iter_to_save];
	  ne_avg_vec_save = new Scalar[iter_to_save];
	  ni_avg_vec_save = new Scalar[iter_to_save];	
	  nn_avg_vec_save = new Scalar[iter_to_save];
	  Ie_FS_vec_save = new Scalar[iter_to_save];
	  Ii_FS_vec_save = new Scalar[iter_to_save];
	  Ie_anode_vec_save = new Scalar[iter_to_save];
	  Ii_anode_vec_save = new Scalar[iter_to_save];
	  In_FS_vec_save = new Scalar[iter_to_save];	
	  In_anode_vec_save = new Scalar[iter_to_save];		
	  Ie_FS_beam_vec_save = new Scalar[iter_to_save];
	  Ii_anode_beam_vec_save = new Scalar[iter_to_save];	
	phi_ij_vec_save = new Scalar[iter_to_save];

	num_bound_Ir = theSpace->get_num_bound_Ir();

	Ir_spc_avg = new Scalar **[iter_to_save];
	for(int i=0;i<iter_to_save;i++) {
		Ir_spc_avg[i] = new Scalar *[num_bound_Ir];
		for(j=0;j<num_bound_Ir;j++) {
			Ir_spc_avg[i][j]=new Scalar[number_of_species];
			//zero the memory
			memset(Ir_spc_avg[i][j],0,(number_of_species)*sizeof(Scalar)); // zero array
		}
	}



	Ir_spc_sum = new Scalar *[num_bound_Ir];
	for(j=0;j<num_bound_Ir;j++) {
		Ir_spc_sum[j]=new Scalar[number_of_species];
		//zero the memory
		memset(Ir_spc_sum[j],0,(number_of_species)*sizeof(Scalar)); // zero array
	}

	
	
	  KE_e_sum = KE_i_sum = KE_n_sum = ne_avg_sum = ni_avg_sum = nn_avg_sum = Ie_FS_sum = Ii_FS_sum = Ie_anode_sum = Ii_anode_sum = In_FS_sum = In_anode_sum = Ii_anode_beam_sum = Ie_FS_beam_sum =  phi_ij_sum = 0.0;
      i_time_step = i_avg = 0;



	// set save file name############################3

	out_dir = theInputFile;

    k = out_dir.find_last_of("/");
	
	out_dir =  out_dir.substr (0,k+1)   +   "diagnostics_history.txt";
	//char* c_out_dir = &out_dir[0];


	/// omri write to file init end
}

Diagnostics::~Diagnostics() {

	BoundaryIhist->deleteAll();
	delete BoundaryIhist;
	BoundaryIhist=0;
	BoundaryPFhist->deleteAll(); 
	delete BoundaryPFhist;
	BoundaryPFhist = 0 ;

  DiagList *dlist = theSpace->getDiagList();
  dlist->deleteAll();

  int sh;
  for ( sh = 0; sh < 2; sh++ ) {
    delete rmsEmittance[sh];
  } 
  delete [] rmsEmittance;

  for (sh = 0; sh < 3; sh++ ) {
    delete aveVelocity[sh];
    delete rmsVelocity[sh]; 
  }
  delete [] aveVelocity;
  delete [] rmsVelocity;

  for (sh = 0; sh < 2; sh++ ) {
    delete aveBeamSize[sh];
    delete rmsBeamSize[sh];
  }
  delete [] aveBeamSize;
  delete [] rmsBeamSize; 

  int isp;
  for(isp=0;isp<number_of_species;isp++) {
      delete number[isp];
      delete ke_species[isp];
      delete ngroups[isp];
      delete total_density[isp]; 
      delete Ave_KE[isp];
    }
    delete [] number; 
    delete [] ke_species;
    delete [] ngroups;
    delete [] total_density;
    delete [] Ave_KE;
    delete [] Ave_KE_gridprob;	//IgalK change
    delete [] number_gridprob;	//IgalK change
    delete [] total_density_gridprob;	//IgalK change
	delete [] Ir_gridprob; //IgalK change

  for(isp=0;isp<number_of_species;isp++) {
    delete [] theSpecies[isp].x1positions;
    delete [] theSpecies[isp].x2positions;
    delete [] theSpecies[isp].x1velocities;
    delete [] theSpecies[isp].x2velocities;
    delete [] theSpecies[isp].x3velocities;
  }
  delete [] theSpecies;

  int j;
  for(j=0;j<=Jm;j++) {
    delete [] S_array[j];
    delete [] Ue[j];
    delete [] Ub[j];		
#ifdef DEBUG_PHI
    delete [] phi_err[j];
#endif
  }
  delete [] S_array;
  delete [] Ue;
  delete [] Ub;
#ifdef DEBUG_PHI
  delete [] phi_err;
#endif  
 
  delete [] x1_array;
  delete [] x2_array;

  delete [] x1_arrayNGD;
  delete [] x2_arrayNGD;


  delete [] KE_e_vec_save;
  delete [] KE_i_vec_save;	
  delete [] KE_n_vec_save;
  delete [] ne_avg_vec_save;	
  delete [] ni_avg_vec_save;		
  delete [] nn_avg_vec_save;	
  delete []	Ie_FS_vec_save;
  delete []	Ii_FS_vec_save;
  delete []	Ie_anode_vec_save;	
  delete []	Ii_anode_vec_save;
  delete []	In_FS_vec_save;	
  delete []	In_anode_vec_save;		
  delete [] Ie_FS_beam_vec_save;
  delete []	Ii_anode_beam_vec_save;	
  delete []	phi_ij_vec_save;	


	for(int i=0; i<iter_to_save; i++)
	{
		for(j=0; j<num_bound_Ir; j++) 
			delete [] Ir_spc_avg[i][j];
		delete [] Ir_spc_avg[i];
	}
	delete [] Ir_spc_avg;



	for(j=0; j<num_bound_Ir; j++) 
		delete [] Ir_spc_sum[j];

	delete [] Ir_spc_sum;
	  	
}

/**********************************************************************
 * a helper function to calculate averages and standard deviations
 * for an array of "random" variables distributed according to
 * the array of weight elements _Pq. It uses a two-pass algorithm.
 * DAD: Thu Oct 19 11:23:14 MDT 2000
 **********************************************************************/
void  Diagnostics::Get_Statistical_Moments(int N, Scalar* W, Scalar* Data,
					   Scalar& Ave, Scalar& Var)
{
  // first calculate the average
  Ave = 0.0;
  int i;
  for ( i = 0; i < N; i++ )
    Ave += W[i] * Data[i];
  // now calculate the variance on the second pass.
  Scalar tmp;
  Scalar sum  = 0.0;
  Scalar sum2 = 0.0;
  for ( i = 0; i < N; i++ ) {
    tmp   = Data[i] - Ave;
    sum  += W[i] * tmp;
    sum2 += W[i] * tmp * tmp;
  }

  Var = sum2 - (sum * sum);
  // check if it is positive
  if ( Var < 0.0 ) {
    if ( Var > -Scalar_EPSILON ) { // a precision problem set Var to zero 
      Var= 0.0;
    } else {
      printf(" ");
      printf("WARNING: Negative variance found in ");
	    printf("Diagnostics::Get_Statistical_Moments(...)");
	    printf("This should NEVER happen!!!");
	    printf("Setting the value of the variance to ZERO!");
    }
  }
  
  return;
}
/**********************************************************************
 * A helper function for the calculation of the averages and rms' of
 * positions, velocities, emittances, and energy for a specific 
 * species. DAD: Mon Oct 23 11:05:04 MDT 2000
 **********************************************************************/
void Diagnostics::Set_Diagnostics_Moments(int Geometry_Flag, 
			       oopicListIter<ParticleGroup>& PGIterator, 
			       Grid* RMSGrid, Scalar* _Pq, int NumWeightElements, 
			       Scalar& AveX1, Scalar& AveX2, 
			       Scalar& RMSX1, Scalar& RMSX2, 
			       Scalar* _aveU, Scalar* _rmsU, Scalar* _rmsEmit,
			       Scalar& _aveE_eV, Scalar& _rmsE_eV)
{
  int numParticles;
  int counter;
  ParticleGroup* RMSParticles = 0;  

  Vector2* position;  // pointer to particle position array

#if (!defined(_MSC_VER) && (  defined(__DECCXX) || defined(_CRAYT3E) || defined(__KCC) || defined(__xlC__) || defined(__ICC)))
  Scalar* X1 = new Scalar(NumWeightElements); // position arrays
  Scalar* X2 = new Scalar(NumWeightElements);
  Scalar* U1 = new Scalar(NumWeightElements); // velocity arrays
  Scalar* U2 = new Scalar(NumWeightElements);
  Scalar* U3 = new Scalar(NumWeightElements);
#else 
#if defined(__BCPLUSPLUS__) || defined(_MSC_VER) 
  Scalar* X1 = new Scalar[NumWeightElements]; // position arrays
  Scalar* X2 = new Scalar[NumWeightElements];
  Scalar* U1 = new Scalar[NumWeightElements]; // velocity arrays
  Scalar* U2 = new Scalar[NumWeightElements];
  Scalar* U3 = new Scalar[NumWeightElements];
#else 
  Scalar X1[NumWeightElements]; // position arrays
  Scalar X2[NumWeightElements];
  Scalar U1[NumWeightElements]; // velocity arrays
  Scalar U2[NumWeightElements];
  Scalar U3[NumWeightElements];
#endif  // __BCPLUSPLUS__
#endif  

  ///Scalar E_eV[NumWeightElements]; // energy array for the particles in a group in eV

  // initialize X1, X2, velocity, and energy arrays.
  Vector3* U;  // pointer to particle velocity array
  counter = 0;
  for(PGIterator.restart(); !PGIterator.Done(); PGIterator++) {
    RMSParticles = PGIterator.current();
    numParticles = RMSParticles->get_n();
    position     = RMSParticles->get_x();
    U            = RMSParticles->get_u();
    for (int j=0; j<numParticles; j++) {
      X1[counter]   = RMSGrid->getMKS( position[j] ).e1();
      X2[counter]   = RMSGrid->getMKS( position[j] ).e2();
      U1[counter]   = U[j].e1();
      U2[counter]   = U[j].e2();
      U3[counter]   = U[j].e3();
      ///E_eV[counter] = RMSParticles->energy_eV(j);
      counter++;
    }
  }
  
  // get the mean and standard deviation values for the positions
  Scalar VarX1, VarX2;
  Get_Statistical_Moments(NumWeightElements, _Pq, X1, AveX1, VarX1);
  Get_Statistical_Moments(NumWeightElements, _Pq, X2, AveX2, VarX2);
  RMSX1 = sqrt( VarX1 );
  if ( Geometry_Flag == 0 ) { // cylindrical geometry according to "grid.h"   
    // average of X2 == r must be zero since it is the average
    // of the vector r. (Bruhwiler/Dimitrov, 10/19/2000)
    VarX2 += AveX2*AveX2; 
    RMSX2 = sqrt( VarX2 );
    AveX2 = 0.; 
  } else if ( Geometry_Flag == 1 ) { // xy geometry
    RMSX2 = sqrt( VarX2 );
  }
  
  // get the mean and standard deviation values for the velocities
  Scalar VarU[3];
  Get_Statistical_Moments(NumWeightElements, _Pq, U1, _aveU[0], VarU[0]);
  Get_Statistical_Moments(NumWeightElements, _Pq, U2, _aveU[1], VarU[1]);
  Get_Statistical_Moments(NumWeightElements, _Pq, U3, _aveU[2], VarU[2]);
  for ( int j = 0; j < 3; j++ ) 
    _rmsU[j] = sqrt(VarU[j]);

  // calculate the rms Emittance for X1 and X2;  
  if ( Geometry_Flag == 0 ) { // cylindrical geometry according to "grid.h"   
    // by default for cylindrical geometry the beam axis is X1
    // this will be calculated only if fabs(_aveU[1]/_aveU[0]) < 1.0e-3
    // and fabs(_aveU[2]/_aveU[0]) < 1.0e-3 which will be use as a condition 
    // for a beam, otherwise the emittance will be set to the negative
    // value of 0.0
    if ( (fabs(_aveU[1]/_aveU[0]) < 1.0e-3) && (fabs(_aveU[2]/_aveU[0]) < 1.0e-3) ) {
      Scalar ave_XU[2];
      int i;
      for ( i = 0; i < 2; i++ ) {
	ave_XU[i] = 0.0;
      }
      for ( i = 0; i < NumWeightElements; i++ ) {
	ave_XU[0] += _Pq[i] * (X1[i] - AveX1) * (U1[i] - _aveU[0]);
	ave_XU[1] += _Pq[i] * (X2[i] - AveX2) * (U2[i] - _aveU[1]);
      }
      
      _rmsEmit[0] = VarX1 * VarU[0] - ave_XU[0] * ave_XU[0];
      _rmsEmit[1] = VarX2 * VarU[1] - ave_XU[1] * ave_XU[1];
      
      for ( i = 0; i < 2; i++ ) {
	if ( _rmsEmit[i] < 0.0 ) {
	  if ( _rmsEmit[i] > -Scalar_EPSILON ) { /* a precision problem,  
						    set _rmsEmit[i] to zero */
	    _rmsEmit[i] = 0.0;
	  } else {
          std::cout << "WARNING: Negative value for the square of the " 
	       << "Emittance E[" << i+1 << "] found in " << endl 
	       << "Diagnostics::Set_RMS_Cylindrical(...)" << endl
	       << "This should NEVER happen!!! Now setting the value of" << endl
	       << "the emittance to ZERO!" << endl;
	  }
	} else {
	  _rmsEmit[i] = sqrt(_rmsEmit[i]) / fabs( _aveU[0]);
	}
      }
    } else {
      for ( int i = 0; i < 2; i++ ) 
	_rmsEmit[i] = 0.0;
    }
  } else if ( Geometry_Flag == 1 ) { // xy geometry
    // for cartesian geometry we have to determine the 
    // direction of the beam, if we have a beam at all. 
    static int BeamDirectionIndex = -1; // assume initially there is no beam
    int AveUmaxIndex;
    int ib = 0;
    ib = (fabs(_aveU[0]) > fabs(_aveU[1]) ) ? 0 : 1;
    AveUmaxIndex = (fabs(_aveU[ib]) > fabs(_aveU[2])) ? ib : 2;
    // is it a beam?
    if ( AveUmaxIndex != ib ) {
      if ( fabs(_aveU[ib]/_aveU[AveUmaxIndex]) < 1.0e-3 )
	BeamDirectionIndex = AveUmaxIndex;
    } else {
      if ( AveUmaxIndex == 0 )
	ib = (fabs(_aveU[1]) > fabs(_aveU[2]) ) ? 1 : 2;
      else if ( AveUmaxIndex == 1 ) 
	ib = (fabs(_aveU[0]) > fabs(_aveU[2]) ) ? 0 : 2;
      if ( fabs(_aveU[ib]/_aveU[AveUmaxIndex]) < 1.0e-3 )
	BeamDirectionIndex = AveUmaxIndex;
    } 
    
    int i;
    if ( BeamDirectionIndex == -1 ) { // no beam
      for ( i = 0; i < 2; i++ ) 
	_rmsEmit[i] = 0.0; // set it ZERO for no beam
    } else { 
      // there is a beam
      // calculate the rms Emittance for X1 and X2;
      Scalar ave_XU[2];
      for ( i = 0; i < 2; i++ ) {
	ave_XU[i] = 0.0;
      }
      for ( i = 0; i < NumWeightElements; i++ ) {
	ave_XU[0] += _Pq[i] * (X1[i] - AveX1) * (U1[i] - _aveU[0]);
	ave_XU[1] += _Pq[i] * (X2[i] - AveX2) * (U2[i] - _aveU[1]);
      }
      
      _rmsEmit[0] = VarX1 * VarU[0] - ave_XU[0] * ave_XU[0];
      _rmsEmit[1] = VarX2 * VarU[1] - ave_XU[1] * ave_XU[1];
      for ( i = 0; i < 2; i++ ) {
	if ( _rmsEmit[i] < 0.0 ) {
	  if ( _rmsEmit[i] > -Scalar_EPSILON ) { /* a precision problem,  
						    set _rmsEmit[i] to zero */
	    _rmsEmit[i] = 0.0;
	  } else {
	  std::cout << "WARNING: Negative value for the square of the " 
	       << "Emittance E[" << i+1 << "] found in " << endl 
	       << "Diagnostics::Set_RMS_Cylindrical(...)" << endl
	       << "This should NEVER happen!!! Now setting the value of" << endl
	       << "the emittance to ZERO!" << endl;
	  }
	} else {
	  _rmsEmit[i] = sqrt(_rmsEmit[i]) / fabs( _aveU[0]);
	}
      }
    }
  }

  // get the mean and standard deviation values for the Energy
  ///Scalar varE;
  ///Get_Statistical_Moments(NumWeightElements, _Pq, E_eV, _aveE_eV, varE);
  ///_rmsE_eV = sqrt(varE);

  // Clean up the locally allocated memory:
#if defined(__DECCXX) || defined(__BCPLUSPLUS__) || defined(_CRAYT3E) || defined(__KCC) || defined(_MSC_VER) || defined(__xlC__)
  delete [] X1;
  delete [] X2;
  delete [] U1;
  delete [] U2;
  delete [] U3;
#endif
  
  return;
} 

/************************************************************************/
/* time history accumulator; calculates and stores all history values,  */
/* and performs combing on history values when low on memory            */
/************************************************************************/
void Diagnostics::history()
{
	register int isp, k, j;

	/*****************************************/
	/*  Fixing the diagnostic arrays         */

	simulation_time=theSpace->getTime();
	int jm=theSpace->getJ();
	int km=theSpace->getK();
	oopicListIter<Ihistdiag> nextdiag(*BoundaryIhist);
	oopicListIter<PFhistdiag> nextdiag2(*BoundaryPFhist);

	DiagList *dlist = theSpace->getDiagList();
	oopicListIter<Diag> nextd(*dlist);

	// ------Updating newDiag-----//
	for(nextd.restart();!nextd.Done(); nextd++) 
		nextd.current()->MaintainDiag((Scalar)simulation_time);
	//--------------------------//

	// update the history objects

	Uktot = 0;
	for(isp=0;isp<number_of_species;isp++) {
		number[isp]->add(theSpecies[isp].nparticles,simulation_time);
                number_gridprob[isp]=theSpecies[isp].nparticles; //IgalK
		ngroups[isp]->add(theSpecies[isp].ngroups,simulation_time);
		ke_species[isp]->add(theSpecies[isp].spec->getKineticEnergy(),simulation_time);
		Uktot += theSpecies[isp].spec->getKineticEnergy();
	}
	energy_e.add(Uetot,simulation_time);
	energy_b.add(Ubtot,simulation_time);
	energy_k.add(Uktot,simulation_time);
	energy_all.add(Uetot+Ubtot+Uktot,simulation_time);
	divderrorhis.add(divderrormag,simulation_time);

	for(nextdiag2.restart();!nextdiag2.Done();nextdiag2++)
		{ //if it's taking diagnostics of Poyting Flux
			Scalar p_flux = *nextdiag2.current()->p_flux;
			nextdiag2.current()->PFhist->add(p_flux,simulation_time);
			nextdiag2.current()->PFlocal->add(p_flux,simulation_time);
		}


	int no_collec_curr = 0; //IgalK number of boundaries collecting currents
	Scalar Ir = 0;

	SpeciesList* speciesList=theSpace->getSpeciesList();
	oopicListIter<Species> spiter(*speciesList);
	
	for(nextdiag.restart();!nextdiag.Done();nextdiag++) { 
	  //if it's taking diagnostics of particles
	  Ir = nextdiag.current()->p_dist->getdQ()/dt;
	  nextdiag.current()->Ihist->add(Ir,simulation_time);
	  for (spiter.restart (); !spiter.Done () ; spiter++) {
		if (spiter.current()->cycleNow()){
			int isp = spiter.current()->getID();
			Ir = nextdiag.current()->p_dist->getdQ(isp)/dt/spiter.current()->get_subcycleIndex();

			nextdiag.current()->Ihist_sp[isp]->add(Ir,simulation_time);

			/*if (abs(Ir)>0 && theSpace->getTime()>4.5e-9)
			  Scalar ddd=1;

			if (abs(Ir)>0 && isp>0)
				Scalar ddd2 = 1; // theSpace.cathode_index*/

			theSpace->set_Ir_spc(no_collec_curr,isp,Ir);

			//Scalar ttt=1; //theSpace->get_i_index()			
		}
	  }
	  no_collec_curr++;  
	}



	
	
	
	// Calculation of the total density
	Scalar total_dens = 0.0;
	Scalar icharge = 0.0;
	Scalar total_phys_par = 0.0;
	for(isp=0;isp<number_of_species;isp++) {
		total_dens = 0.0;
		total_phys_par = 0.0;
		icharge = 1/theSpecies[isp].spec->get_q(); 
		for (j=0; j<=jm; j++){
			for (k=0; k<=km; k++){
				// total_dens += CellVolumes[j][k]*rho_species[isp][j][k]; // omri change
				total_dens +=rho_species[isp][j][k]; // omri change
				total_phys_par += CellVolumes[j][k]*rho_species[isp][j][k];
			}
		}
		//total_dens *= icharge; // omri change
		total_dens = total_dens*icharge/(km+1)/(jm+1); // omri change - eventually displays the average number density for each specie
		total_phys_par*=icharge;
		theSpace->set_n_avg_species(isp,total_dens);
		
		total_density_gridprob[isp] = total_dens; //IgalK change
		total_density[isp]->add(total_dens,simulation_time);
		if (total_phys_par) {
			Ave_KE[isp]->add(theSpecies[isp].spec->getKineticEnergy()/(fabs(ELECTRON_CHARGE)*total_phys_par), simulation_time);
			Ave_KE_gridprob[isp] = theSpecies[isp].spec->getKineticEnergy()/(fabs(ELECTRON_CHARGE)*total_phys_par); //IgalK change
		}                                           
		else {
			Ave_KE[isp]->add(0.,simulation_time);
			Ave_KE_gridprob[isp] = 0;
		}
		theSpace->set_KE_species(isp,Ave_KE_gridprob[isp]);
	}

/* add omri print diagnostics to file @@@@@@@@@@@@@@@@22    // write to file */	

// calculation of the omri diagnostics	
		i_time_step++;
		ne_avg_sum+=theSpace->get_n_avg_species(theSpace->get_e_index());
		ni_avg_sum+=theSpace->get_n_avg_species(theSpace->get_i_index());
		nn_avg_sum+=theSpace->get_n_avg_species(theSpace->get_n_index());	
		KE_e_sum+=theSpace->get_KE_species(theSpace->get_e_index());
		KE_i_sum+=theSpace->get_KE_species(theSpace->get_i_index());
		KE_n_sum+=theSpace->get_KE_species(theSpace->get_n_index());
	    Ie_FS_sum += theSpace->getIr_spc(theSpace->get_cathode_index(),theSpace->get_e_index());
		Ii_FS_sum += theSpace->getIr_spc(theSpace->get_cathode_index(),theSpace->get_i_index()) ;
		Ie_anode_sum += theSpace->getIr_spc(theSpace->get_anode_index(),theSpace->get_e_index());
		Ii_anode_sum += theSpace->getIr_spc(theSpace->get_anode_index(),theSpace->get_i_index());
		In_FS_sum += theSpace->getIr_spc(theSpace->get_cathode_index(),theSpace->get_n_index());
		In_anode_sum += theSpace->getIr_spc(theSpace->get_anode_index(),theSpace->get_n_index());	
		Ie_FS_beam_sum += theSpace->getIbeam_spc(theSpace->get_cathode_index(),theSpace->get_e_index());
		Ii_anode_beam_sum += theSpace->getIbeam_spc(theSpace->get_anode_index(),theSpace->get_i_index());	
	 phi_ij_sum += phi[0][0];


	 for(int i_bound=0; i_bound<num_bound_Ir; i_bound++){  //if it's taking diagnostics of particles
		 for(isp=0;isp<number_of_species;isp++) {
			 Ir_spc_sum[i_bound][isp]+= theSpace->getIr_spc(i_bound,isp);
		 }
	 }
	 

	 
	
		if (i_time_step%iter_to_avg == 0) {
			ne_avg_vec_save[i_avg] = ne_avg_sum/iter_to_avg;
			ne_avg_sum=0;
			ni_avg_vec_save[i_avg] = ni_avg_sum/iter_to_avg;
			ni_avg_sum=0;
			nn_avg_vec_save[i_avg] = nn_avg_sum/iter_to_avg;
			nn_avg_sum=0;
			
			KE_e_vec_save[i_avg] = KE_e_sum/iter_to_avg;
			KE_e_sum=0;
			KE_i_vec_save[i_avg] = KE_i_sum/iter_to_avg;
			KE_i_sum=0;			
			KE_n_vec_save[i_avg] = KE_n_sum/iter_to_avg;
			KE_n_sum=0;

			Ie_FS_vec_save[i_avg] = Ie_FS_sum/iter_to_avg;
			Ie_FS_sum=0;			
			Ii_FS_vec_save[i_avg] = Ii_FS_sum/iter_to_avg;
			Ii_FS_sum=0;	
			Ie_anode_vec_save[i_avg] = Ie_anode_sum/iter_to_avg;
			Ie_anode_sum=0;	
			Ii_anode_vec_save[i_avg] = Ii_anode_sum/iter_to_avg;
			Ii_anode_sum=0;
			In_FS_vec_save[i_avg] = In_FS_sum/iter_to_avg;
			In_FS_sum=0;	
			In_anode_vec_save[i_avg] = In_anode_sum/iter_to_avg;
			In_anode_sum=0;				
			Ii_anode_beam_vec_save[i_avg] = Ii_anode_beam_sum/iter_to_avg;
			Ii_anode_beam_sum=0;
			Ie_FS_beam_vec_save[i_avg] = Ie_FS_beam_sum/iter_to_avg;
			Ie_FS_beam_sum=0;			
			phi_ij_vec_save[i_avg] = phi_ij_sum/iter_to_avg;
			phi_ij_sum=0;	



			for(j=0;j<num_bound_Ir;j++) {
				for(isp=0;isp<number_of_species;isp++) {
					Ir_spc_avg[i_avg][j][isp]=Ir_spc_sum[j][isp] / iter_to_avg;
					//zero the memory
					Ir_spc_sum[j][isp]=0;
				}
			}
			
			
			i_avg++;
		}


		if (i_time_step%iter_to_save == 0) {

		FILE *myfile_con;
		myfile_con = fopen(out_dir.c_str(),"a");		

		int size = ftell(myfile_con); // see if file is empty (size 0)

			if(size==0) {
				//fprintf(myfile_con,"   t_avg[sec]      elc_avg_n[m^-3]              ion_avg_n[m^-3]             ntr_avg_n[m^-3]     	      elc_KE[ev]           ion_KE[ev]            ntr_KE[ev]		    Ie_FS[A]        Ii_FS[A]       Ie_anode[A]        Ii_anode[A]     In_FS[A]     In_anode[A]    Ie_FS_beam[A]        Ii_anode_beam[A]       phi_ij[V]       I_boun_spc_arr[A] , dt=%E \n", dt); // headline
				fprintf(myfile_con,"   t_avg[sec]      elc_avg_n[m^-3]              ion_avg_n[m^-3]             ntr_avg_n[m^-3]     	      elc_KE[ev]           ion_KE[ev]            ntr_KE[ev]		    Ie_FS[A]        Ii_FS[A]       Ie_anode[A]        Ii_anode[A]     In_FS[A]     In_anode[A]    Ie_FS_beam[A]        Ii_anode_beam[A]"); // headline
				//fprintf(myfile_con,"   t_avg[sec]      elc_avg_n[m^-3]              ion_avg_n[m^-3]             ntr_avg_n[m^-3]     	      elc_KE[ev]           ion_KE[ev]            ntr_KE[ev]");
				oopicListIter<Boundary> nextb(*blist);
				for(nextb.restart();!nextb.Done();nextb++) { //if it's taking diagnostics of particles
					if(nextb.current()->get_particle_diag()) {
						for(isp=0;isp<number_of_species;isp++) 
							fprintf(myfile_con,"         I_%s_%s [A]",theSpecies[isp].name,nextb.current()->getBoundaryName().c_str());
					}
				}


				fprintf(myfile_con, "       phi_ij[V]  		   dt=%E \n", dt); // headline
			}


			
			for (int i_save=0; i_save<i_avg; i_save++){
				
				fprintf (myfile_con, "%.10E         ", t_avg) ;				
				fprintf (myfile_con, "%.3E        ",ne_avg_vec_save[i_save]);			
				fprintf (myfile_con, "%.3E        ",ni_avg_vec_save[i_save]);			
				fprintf (myfile_con, "%.10E        ",nn_avg_vec_save[i_save]);				
				fprintf (myfile_con, "%.3E         ", KE_e_vec_save[i_save]) ;
				fprintf (myfile_con, "%.3E         ", KE_i_vec_save[i_save]) ;
				fprintf (myfile_con, "%.3E         ", KE_n_vec_save[i_save]) ;
				fprintf (myfile_con, "%.3E         ", Ie_FS_vec_save[i_save]) ;
				fprintf (myfile_con, "%.3E         ", Ii_FS_vec_save[i_save]) ;
				fprintf (myfile_con, "%.3E         ", Ie_anode_vec_save[i_save]) ;				
				fprintf (myfile_con, "%.3E         ", Ii_anode_vec_save[i_save]) ;
				fprintf (myfile_con, "%.3E         ", In_FS_vec_save[i_save]) ;				
				fprintf (myfile_con, "%.3E         ", In_anode_vec_save[i_save]) ;								
				fprintf (myfile_con, "%.3E         ", Ie_FS_beam_vec_save[i_save]) ;
				fprintf (myfile_con, "%.3E         ", Ii_anode_beam_vec_save[i_save]) ;		
				for(int i_bound=num_bound_Ir-1; i_bound>-1; i_bound--) { //for some reason the code reverse the boundary order , and we need to reverse back on printing
					for(isp=0;isp<number_of_species;isp++) 
						fprintf (myfile_con, "%.3E         ", Ir_spc_avg[i_save][i_bound][isp]) ;

				}
				fprintf (myfile_con, "%.3E         \n", phi_ij_vec_save[i_save]) ;	

				
							
				t_avg+=iter_to_avg*dt;
			}
		fclose (myfile_con);
		i_avg=0;

		}

	/* add omri print end@@@@@@@@@@@2*/
	

	/***************************begin********************************
	 * Calculation of the rms beam size (presumably for a beam)
	 * (Bruhwiler/Dimitrov, 10/09/2000)
	 ****************************************************************/
	// This is done for species number iSpecesRMS, assuming that the
	// rmsBeamSizeFlag has been set to 1 (default value 0) for at
	// least one species.  If the flag has been set for two different
	// species, only the first one will get plotted.
	// We set the initial value to -1, so that nothing is plotted,
	// if the rmsBeamSizeFlag has not been set to 1 for any species.
	int iSpeciesRMS = -1;
	
	// We will loop over all species in the simulation.
	// For each species, we will check whether the
	// rmsDiagnosticsFlag has been set to 1 (default value is 0).
	oopicList<ParticleGroup> *pgList;
	oopicListIter<ParticleGroup> pgIterator;
	ParticleGroup* tmpGroup;
	Species* tmpSpecies;
	int tmpFlag;
	BOOL Species_Flag = false;
	isp = 0; 
	while ( (isp < number_of_species ) && (!Species_Flag) ) {
	  oopicListIter<ParticleGroup> 
	    pgscan(*theSpace->getParticleGroupListArray()[isp]);
	  pgscan.restart();
	  if ( pgscan.isEmpty() ) {
	    tmpGroup = pgscan.current();
	    tmpSpecies = tmpGroup -> get_species();
	    tmpFlag = tmpSpecies -> get_rmsDiagnosticsFlag();
	    
	    if ( tmpFlag ) {
	      iSpeciesRMS = isp;
	      pgList = theSpace->getParticleGroupListArray()[iSpeciesRMS];
	      pgIterator = oopicListIter<ParticleGroup>(*pgList);
	      Species_Flag = true;
	    }
	  }
	  isp++;
	}

	if (iSpeciesRMS != -1) {
	  // initialize a Scalar array with the 
	  // normalized the weight elements (distribution density)
	  // 
	  // (i) count the number of weight elements
	  ParticleGroup* RMSParticles = 0;
	  int numParticles; 
	  int numWeightElements = 0;
	  for(pgIterator.restart(); !pgIterator.Done(); pgIterator++) {
	    RMSParticles = pgIterator.current();
	    numWeightElements += RMSParticles->get_n();
	  }
	  // (ii) allocate memory for the weight elements
	  Scalar* Pq = new Scalar[numWeightElements];
	  // (iii) initialize Pq and the normalization constant
	  numWeightElements = 0;
	  Scalar sum_q   = 0.;
	  for(pgIterator.restart(); !pgIterator.Done(); pgIterator++) {
	    RMSParticles = pgIterator.current();
	    numParticles = RMSParticles->get_n();
	    for (int j=0; j<numParticles; j++) {
	      Pq[numWeightElements] = fabs( RMSParticles->get_q(j) );
	      sum_q   += Pq[numWeightElements];
	      numWeightElements++;
	    }
	  }
	  // (iv) normalize Pq
	  for (int j=0; j< numWeightElements; j++) {
	      Pq[j]  /= sum_q;
	  }

	  // the Grid knows what geometry it is defined for via the variable
	  // "int geometryFlag" which is a data member of the Grid class. 
	  // By default this variable is set to 0 == cylindrical geometry.
	  // It can be set to cartesian geometry via adding the line
	  // "Geometry = 1" in the input file when the Grid parameters are 
	  // specified. DAD::Wed Oct 11 14:59:01 MDT 2000
	  Grid* rmsGrid = theSpace->get_grid();
	  int tmp_Geometry_Flag = rmsGrid->query_geometry();
	  
	  // Next, we use helper functions to calculate the Ave & RMS values.
	  Scalar aveX1;
	  Scalar aveX2;
	  Scalar rmsX1;
	  Scalar rmsX2;
	  Scalar aveU[3];
	  Scalar rmsU[3];
	  Scalar rmsEmit[2]; // for the Emittance
	  Scalar aveE_eV;    // for the average and rms Energy in eV
	  Scalar rmsE_eV;
	  
	  if ( (tmp_Geometry_Flag != 0) && (tmp_Geometry_Flag != 1) ) { 
	    printf("WARNING: Undefined geometry flag detected in diagn.cpp");
	  } else {
	    Set_Diagnostics_Moments(tmp_Geometry_Flag, pgIterator, rmsGrid, Pq, 
				    numWeightElements, aveX1, aveX2, rmsX1, rmsX2, 
				    aveU, rmsU, rmsEmit, aveE_eV, rmsE_eV);   
	  }
          
	  aveBeamSize[0]->add(aveX1,simulation_time);
	  aveBeamSize[1]->add(aveX2,simulation_time);	  
	  rmsBeamSize[0]->add(rmsX1,simulation_time);
	  rmsBeamSize[1]->add(rmsX2,simulation_time);

	  aveVelocity[0]->add(aveU[0],simulation_time);
	  aveVelocity[1]->add(aveU[1],simulation_time);
	  aveVelocity[2]->add(aveU[2],simulation_time);
	  
	  rmsVelocity[0]->add(rmsU[0],simulation_time);
	  rmsVelocity[1]->add(rmsU[1],simulation_time);
	  rmsVelocity[2]->add(rmsU[2],simulation_time);

	  rmsEmittance[0]->add(rmsEmit[0],simulation_time);
	  rmsEmittance[1]->add(rmsEmit[1],simulation_time);

	  ///aveEnergy_eV->add(aveE_eV,simulation_time);
	  ///rmsEnergy_eV->add(rmsE_eV,simulation_time);
	  
	  delete[] Pq;
	}
	/****************************************************************
	 * (end of changes by Bruhwiler/Dimitrov, 10/09/2000)
	 ****************************************************************/

	return;
}

/* another function for testing. */
void Diagnostics::compute_div_derror_magnitude(Scalar **divderror,int jm,int km) {
	int j,k;
	Scalar rhomagave;
	divderrormag=0;
	rhomagave = 0;
	for(j=1;j<jm-1;j++)
		for(k=1;k<km-1;k++) {
			rhomagave+=fabs(rho[j][k]);
			divderrormag+=fabs(divderror[j][k]);
		}
	if(rhomagave!=0) divderrormag/=rhomagave;
}


/*  UpdatePreDiagnostics does processing on raw physics data such
	which should be done before a physics iteration so that the user
	gets fields + particle positions which are reasonably synchronized.
 */

void Diagnostics::UpdatePreDiagnostics()
{  
	int i,j;

#ifndef PHI_TEST

/* parallel data*/
		double t00=omp_get_wtime();	  
//		int num_threads_used=-1;
//		int num_threads_count=0;	

		//int N_threads_local = 4 ;
		//omp_set_num_threads(N_threads_local);

/*parallel data end*/
	
	if(AllDiagnostics) {  //collect diagnostics at all?
		 
		if(PhaseSpacePlots!=0.0) {
			long int nparticles_thread[N_threads];
			long int nparticles_thread_cum[N_threads+1];
			
#pragma omp parallel 
			{	
			int My_thread_ID = omp_get_thread_num();
			register int jj;
			int flag_break; // indicating weather to continue running on particle groups
			int n_local;

			register int i_local;
			Grid *grid_local = theSpace->get_grid();
			register Vector2 *xarr;
			register Vector3 *varr;
			register Vector2 x ;
			register Vector3 v ;
			//register SpeciesDiag *thisSpc_local= thisSpc;
			register Vector2 ** X_local = 	theSpace->get_X_local_thread(My_thread_ID);				
			int isp;
			int skip;
			long int n_groups_local  ;
			long int n_particles_local  ;
			//  maintain the phase space diagnostics.
			
			for(isp = 0; isp < number_of_species; isp++) {
				
				register SpeciesDiag *thisSpc_local= &theSpecies[isp];
#pragma omp single
				{
				thisSpc_local->nparticles=0;
				thisSpc_local->ngroups=0;
			}
				n_groups_local = 0 ;
				n_particles_local = 0 ;
				register oopicListIter<ParticleGroup>  pgscan_local(*theSpace->getParticleGroupListArray()[isp]);
				//i=0;  //place in phase space diagnostic
				
				// First count the particles, store in thisSpc=theSpecies[isp]
				flag_break=0;
				pgscan_local.restart();

					if (pgscan_local.Done() )  // case 0 particles
						flag_break=1;
		  

						for (jj=0; jj<My_thread_ID ; jj++){
							pgscan_local++;
							if (pgscan_local.Done() ){
								flag_break=1;
								break;
								}
						}
				

  					while (flag_break==0){ 

					n_particles_local+=(int)pgscan_local.current()->get_n();
					n_groups_local++;
						  
						for (jj=0; jj<N_threads ; jj++){
							pgscan_local++;
							if (pgscan_local.Done()){
								flag_break=1;
								break;
							}
						}
					}// close while
				

#pragma omp atomic 
				thisSpc_local->ngroups+=n_groups_local;
#pragma omp atomic 
				thisSpc_local->nparticles+=n_particles_local;
				
				nparticles_thread[My_thread_ID] = n_particles_local;

#pragma omp barrier
				// Skip particles in the physics so we don't overflow the
				// particle array in the diagnostics
				skip = thisSpc_local->nparticles/thisSpc_local->memory_size +1;
#pragma omp single 
				{
					//theSpecies[isp].spec->set_ngroups(thisSpc_local->ngroups); 
					thisSpc_local->nparticles_plot = thisSpc_local->nparticles/skip;
					nparticles_thread_cum[0]=0;					
					for (jj=0; jj<N_threads ; jj++){
					nparticles_thread_cum[jj+1] = nparticles_thread_cum[jj] + nparticles_thread[jj];
					
					}

				}
#ifdef MPI_VERSION
				Scalar l_nparticles = thisSpc_local->nparticles;
				Scalar t_nparticles = 0;
				MPI_Reduce ( &l_nparticles, &t_nparticles,1,MPI_Scalar,MPI_SUM,0,XOOPIC_COMM);
				if(MPI_RANK==0)
				  thisSpc->nparticles = (int)t_nparticles;
#endif /* MPI_VERSION */

				/* orig code start@@@@@@@@@
				
				for(pgscan.restart();!pgscan.Done();pgscan++) {
					int npart =(int) pgscan.current()->get_n();
					for(j=0;j<npart&&i<NPARTICLES;j+=skip,i++) {
							Vector2 *xarr=pgscan.current()->get_x();
							Vector3 *varr=pgscan.current()->get_u();
							Vector2 x = theSpace->getMKS(xarr[j]);
							Vector3 v = varr[j];
						
							thisSpc->x1positions[i]=x.e1();
							thisSpc->x2positions[i]=x.e2();
							thisSpc->x1velocities[i] = v.e1();
							thisSpc->x2velocities[i] = v.e2();
							thisSpc->x3velocities[i] = v.e3();
						}
				}

				 orig code end@@@@@@@@@*/
				

/* omri parallel start @@@@@@@@@@@*/



				//oopicListIter<ParticleGroup> pgscan_local(*theSpace->getParticleGroupListArray()[isp]);
				flag_break=0;
				i_local = nparticles_thread_cum[My_thread_ID];

				pgscan_local.restart();

				if (pgscan_local.Done() )  // case 0 particles
					flag_break=1;


				for (jj=0; jj<My_thread_ID ; jj++){
					pgscan_local++;
					if (pgscan_local.Done() ){
						flag_break=1;
						break;
					}
				}


				
  					while (flag_break==0){ 
						//collect_charge_to_grid(pgscan_local(),rho_local[ID]);
						n_local=pgscan_local.current()->get_n();				  

						for (jj=0; jj<n_local && i_local<NPARTICLES; jj+=skip){

							xarr=pgscan_local.current()->get_x();
							varr=pgscan_local.current()->get_u();
							x = grid_local->getMKS_local(xarr[jj], X_local);
							v = varr[jj];
						
							thisSpc_local->x1positions[i_local]=x.e1();
							thisSpc_local->x2positions[i_local]=x.e2();
							thisSpc_local->x1velocities[i_local] = v.e1();
							thisSpc_local->x2velocities[i_local] = v.e2();
							thisSpc_local->x3velocities[i_local] = v.e3();

							i_local++;
						}

						  
						for (jj=0; jj<N_threads ; jj++){
							pgscan_local++;
							if (pgscan_local.Done()){
								flag_break=1;
								break;
							}
						}
					}// close while
					
				} // close specie loop

				
 /*omri parallel end @@@@@@@@@@@*/

				


				
			} // close #pragma omp
		} // end if PhaseSpacePlots
		
	} // end if 2
#endif


//double t1 = omp_get_wtime()-t00;
//static double t1=0;
//t1 += omp_get_wtime()-t00;
//double speed_up_f=0.365/t1;
//printf (" The speed up factor is %f with %d threads \n",speed_up_f,N_threads_local); //

//printf (" execution time is %f with %d threads \n",t1,N_threads); //	
}

void Diagnostics::UpdatePostDiagnostics() {
  int j,k;
  int jm = theSpace->getJ();
  int km = theSpace->getK();

	if (boltzmannFlag)
		theSpace->CopyBoltzmannRho();

  //  initialize time history varaibles to zero
  Uetot=0; Ubtot=0;
  if(AllDiagnostics) {  //collect diagnostics at all?
		Grid* grid = theSpace->get_grid();
		theSpace->updateDivDerror();

		if(EnerPlots){
			for( j=0;j<=jm;j++)
				for( k=0;k<=km;k++){

					/* get the E energy density */
                                  Ue[j][k] = 0.5*(theSpace->get_fields()->get_epsi(MIN(j,jm-1),MIN(k,km-1))+
                                                  theSpace->get_fields()->get_epsi(MAX(j-1,0),MIN(k,km-1))+
                                                  theSpace->get_fields()->get_epsi(MIN(j,jm-1),MAX(k-1,0))+
                                                  theSpace->get_fields()->get_epsi(MAX(j-1,0),MAX(k-1,0)))/4.0
                                    *(E[j][k]*E[j][k]);
					/* get the B energy density */
					Ub[j][k]=sqr(BDynamic[j][k].e1());
					Ub[j][k]+=sqr(BDynamic[j][k].e2());
					Ub[j][k]+=sqr(BDynamic[j][k].e3());
					Ub[j][k]*=0.5*iMU0;

#ifdef DEBUG_PHI
					if(phi[j][k] > 1e-20)
						//absolute error (uncomment one)
					phi_err[j][k] = fabs((phi[j][k] - (sin( M_PI * grid->getX()[j][k].e1()) 
						* sin(M_PI*grid->getX()[j][k].e2()))));

						//relative error
//					phi_err[j][k] = fabs((phi[j][k] - (sin( M_PI * grid->getX()[j][k].e1()) 
//						* sin(M_PI*grid->getX()[j][k].e2())))/phi[j][k]);
					else phi_err[j][k] = 0;
#endif
                                          //collect totals of energy for volumetric data
					Uetot+=Ue[j][k]*CellVolumes[j][k];
					Ubtot+=Ub[j][k]*CellVolumes[j][k];
				}
			
			// This section makes the poyting vector on the interior points. 
			// To get the edge points it is a little harder.
			if (!electrostaticFlag){
				for( j=0;j<=(jm-1);j++)
					for( k=1;k<=(km-1);k++)
						S_array[j][k] = ((E[j][k].cross(BDynamic[j][k]))).
							jvMult(grid->dS(j,k))*iMU0;
			}
		}

		history();
		
  }
} 
#ifdef BENCHMARK
//  This function writes a 'trace', 
//  the first component of the electric field
//  so that you can check the output of a simulation
//  which has no GUI.

void write_validation() {
  FILE *trace_file;
  int J = theSpace->getJ();
  int K = theSpace->getK();

  if((trace_file=fopen("trace.dat","w"))==NULL) {
    stringstream ss (stringstream::in | stringstream::out);
    ss<< "diagn::write_validation: Error: \n"<<
      "can not open trace.dat \n"<<endl;
    
    std::string msg;
    ss >> msg;
    Oops oops(msg);
    throw oops;    // exit()  not called

  }
  for(int j=0;j<J;j++) 
    for(int k=0;k<K;k++) 
			fprintf(trace_file,"%10.4g\n",E[j][k].e1());

  fclose(trace_file);
}
#endif


#ifdef HAVE_HDF5
void Diagnostics::dumpGridH5(dumpHDF5 *H5DumpObj)  
{
  
  //cerr << "Dumping grid information to diagnostic dump file.\n";
  int rank = 2;
  int *size;
  int J,K;
  
// dump grid information
  J = theSpace->getJ();
  K = theSpace->getK();
  rank = 1;
  size = new int[rank];
  
  size[0] = J+1;
  H5DumpObj->writeSimple("x1array",x1_array,rank,size);
  
  size[0] = K+1;
  H5DumpObj->writeSimple("x2array",x2_array,rank,size);
  
  delete[] size;
  return;
}

void Diagnostics::dumpAllDiagsH5(oopicList<dumpHDF5> *H5DiagDumpObjs) 
{

// create iterator over dump objects, increment at end of loop
  oopicListIter<dumpHDF5> nObj(*H5DiagDumpObjs);
// for now assume all diagnostics output to same file, so only dump grid info once
  nObj.restart();
  dumpGridH5(nObj.current() );
 
    for (nObj.restart(); !nObj.Done(); nObj++){
//      cerr << "datasetName "<<(nObj.current())->datasetName<<"\n";
      dumpDiagH5(nObj.current());
    }

  return;
}


void Diagnostics::dumpDiagH5(dumpHDF5 *H5DumpObj) 
{

//  hid_t fileId, groupId, datasetId,dataspaceId;
//  herr_t status;
      
  string diagName = H5DumpObj->datasetName;
// get rid of null character in diagName if present
  if(diagName.find('\0') != string::npos){
    diagName = diagName.substr(0, diagName.length()-1);
  }
  //  cerr << "Dumping "<<diagName<<" diagnostic.\n";

    Scalar *data; 
    int rank = 2;
    int *size;
    int J,K;

    //  dumpGridH5(H5DumpObj);

    J = theSpace->getJ();
    K = theSpace->getK();
  
    if(diagName=="E1" || diagName=="E2" || diagName=="E3" || 
       diagName=="B1" || diagName=="B2" || diagName=="B3" || 
       diagName=="I1" || diagName=="I2" || diagName=="I3" || diagName=="ESPotential"){
      rank = 3;
      size = new int[rank];
      
      size[0] = J+1; size[1] = K+1; size[2] = 1;
      data = new Scalar[(J+1)*(K+1)*1];
      
    
      
      
      if(diagName=="E1"){
	for (int j=0; j<=J;j++)
	  for (int k=0; k<=K;k++){
	    data[j*(K+1) + k] = E[j][k].e1();
	  }
      }else
	if(diagName=="E2"){
	  for (int j=0; j<=J;j++)
	    for (int k=0; k<=K;k++){
	      data[j*(K+1) + k] = E[j][k].e2();
	    }
	  
      }else  
      if(diagName=="E3"){
	for (int j=0; j<=J;j++)
	  for (int k=0; k<=K;k++){
	    data[j*(K+1) + k] = E[j][k].e3();
	  }

      }else
      if(diagName=="B1"){
	for (int j=0; j<=J;j++)
	  for (int k=0; k<=K;k++){
	    data[j*(K+1) + k] = B[j][k].e1();
	  }

      }else
	if(diagName=="B2"){
	  for (int j=0; j<=J;j++)
	    for (int k=0; k<=K;k++){
	    data[j*(K+1) + k] = B[j][k].e2();
	    }
	  
	}else
	  if(diagName=="B3"){
	    for (int j=0; j<=J;j++)
	      for (int k=0; k<=K;k++){
		data[j*(K+1) + k] = B[j][k].e3();
	      }
	  }else
	if(diagName=="ESPotential"){
	  for (int j=0; j<=J;j++)
	    for (int k=0; k<=K;k++){
	      data[j*(K+1) + k] = phi[j][k];
	    }
	}else
      if(diagName=="I1"){
	for (int j=0; j<=J;j++)
	  for (int k=0; k<=K;k++){
	    data[j*(K+1) + k] = I[j][k].e1();
	  }

      }else
	if(diagName=="I2"){
	  for (int j=0; j<=J;j++)
	    for (int k=0; k<=K;k++){
	    data[j*(K+1) + k] = I[j][k].e2();
	    }
	  
	}else
	  if(diagName=="I3"){
	    for (int j=0; j<=J;j++)
	      for (int k=0; k<=K;k++){
		data[j*(K+1) + k] = I[j][k].e3();
	      }
	  }    
    
      // H5DumpObj->createFile();
      H5DumpObj->writeSimple(diagName,data,rank,size);
      
      delete[] size;
      delete[] data;
      return;
   
    }else
      if(diagName=="NGD"){
	rank = 3;
	size = new int[rank];
	
	size[0] = J; size[1] = K; size[2] = 1;
	data = new Scalar[(size[0])*(size[1])*(size[2])];
	
	
	if(diagName=="NGD"){
	  
	  // for now take the first NDG data pointer
	  oopicListIter<NGD> NGDIter(*ptrNGDList);
	  NGDIter.restart();  
	  
	  for (int j=0; j<J;j++)
	    for (int k=0; k<K;k++){
	      data[j*(K) + k] = ((NGDIter.current())->getNGDdata())[j][k];
	    }
	}
	  
	  
	H5DumpObj->writeSimple(diagName,data,rank,size);
	
	  
	// write the grid points for NGD
	string datasetName;
	datasetName ="/NGDx1array";
	rank = 1;
	size[0] = J;
	H5DumpObj->writeSimple(datasetName,x1_arrayNGD,rank,size);
	  
	datasetName = "/NGDx2array";
	size[0] = K;
	H5DumpObj->writeSimple(datasetName,x2_arrayNGD,rank,size);	
	
	
	
	delete[] size;
	delete[] data;
	return;
      }else
	for(int isp=0;isp<number_of_species;isp++) {
// loop over all species, dumping number density for each
// if it exists 
	  cerr << "diagName : "<<diagName<<"\n";
	  string nm = theSpecies[isp].name;
	  nm = nm + "_numberDensity";	
	  if(diagName == nm){
	    rank = 3;
	    size = new int[rank];
	    size[0] = J+1; size[1] = K+1; size[2] = 1;
	    data = new Scalar[(J+1)*(K+1)*1];
	    for (int j=0; j<=J;j++)
	      for (int k=0; k<=K;k++){
		data[j*(K+1) + k] = rho_species[isp][j][k];
	      }

	    cerr << "data: "<< data[size[0]-1]<<"\n";

	    H5DumpObj->writeSimple(nm,data,rank,size);
	   
	  
	  delete[] size;
	  delete[] data;
	  return;
	  }
	}
    
    
    
    DiagList *dlist = theSpace->getDiagList();
    oopicListIter<Diag> nextd(*dlist);

    int iii =0;
    for(nextd.restart();!nextd.Done(); nextd++){ 
      
      cerr << "names of diag "<<iii<<": "<< nextd.current()->getVarName().str()<<"\n";
      iii++;
}

    for(nextd.restart();!nextd.Done(); nextd++){ 
      cerr << "next names of diag "<<iii<<": "<< nextd.current()->getVarName().str()<<"\n";
      if( nextd.current()->getVarName().str()==diagName){
	nextd.current()->dumpDiagH5(H5DumpObj);
      }
    }
    return;
    
    
    cerr << "Input file calls for dump of " << diagName
	 << " diagnostic, which does not exist.\n\n";
    return;
}



void Diagnostics::dumpH5(dumpHDF5 &dumpObj) 
{
  //   cerr << "Entered Diagnostics::dumpH5(dumpObj) .\n";
 
  //   hid_t fileId, groupId, datasetId,dataspaceId;
//     herr_t status;
//     hid_t scalarType = dumpObj.getHDF5ScalarType();


//     fileId= H5Fopen(dumpObj.getDumpFileName().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
//     groupId = H5Gcreate(fileId,"Diagnostics",0);
//     int J = theSpace->getJ();
//     int K = theSpace->getK();
//     hsize_t fieldDim[3];
//     fieldDim[0] = J+1; fieldDim[1] = K+1; fieldDim[2] = 1;
   
  
//     hid_t EspaceId = H5Screate_simple(3, fieldDim, NULL);
//     hid_t EsetId = H5Dcreate(groupId, "E1", scalarType, EspaceId, H5P_DEFAULT);

//     Scalar *Edata = new Scalar[(J+1)*(K+1)*1];
   
//     for (int j=0; j<=J;j++)
//       for (int k=0; k<=K;k++){
//         Edata[j*(K+1) + k] = E[j][k].e1();
//         //   Edata[j*(K+1)*3 + k*3 + 1] =  E[j][k].e2();
//         //  Edata[j*(K+1)*3 + k*3 + 2] = E[j][k].e3();
//       }
//     status = H5Dwrite(EsetId, scalarType, H5S_ALL, EspaceId, H5P_DEFAULT, Edata);
//     status = H5Dclose(EsetId);
   
//     EsetId = H5Dcreate(groupId, "E2", scalarType, EspaceId, H5P_DEFAULT);
//     for (int j=0; j<=J;j++)
//       for (int k=0; k<=K;k++){
//         Edata[j*(K+1) + k] =  E[j][k].e2();
//        //  Edata[j*(K+1)*3 + k*3 + 2] = E[j][k].e3();
//  	    }
//     status = H5Dwrite(EsetId, scalarType, H5S_ALL, EspaceId, H5P_DEFAULT, Edata);
//     status = H5Dclose(EsetId);
   
//     EsetId = H5Dcreate(groupId, "E3", scalarType, EspaceId, H5P_DEFAULT);
//     for (int j=0; j<=J;j++)
//       for (int k=0; k<=K;k++){
//         Edata[j*(K+1) + k] = E[j][k].e3();
//       }
//     status = H5Dwrite(EsetId, scalarType, H5S_ALL, EspaceId, H5P_DEFAULT, Edata);
//     status = H5Sclose(EspaceId); 
//     status = H5Dclose(EsetId);
   
//     status = H5Gclose(groupId);
//     status = H5Fclose(fileId);

  return;
}










void Diagnostics::dumpHistH5(dumpHDF5 &H5DumpObj) // by Matteo
{
	//create the diagnostic group. if you don't number it and there is more than one region per instance it will be broken

	 hid_t fileId, groupId, histgroupId, datasetId, dataspaceId, bgroupId;
	  herr_t status;
	    hid_t scalarType = H5DumpObj.getHDF5ScalarType();
	  hsize_t RANK  = 3;

	  fileId= H5Fopen(H5DumpObj.getDumpFileName().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	     groupId = H5Gcreate(fileId,"History",0);

	     	 int length = *number[0]->get_hi();

	 		hsize_t histLength;//length;
	 		histLength = length;

	 		string datasetName;


	     //time vector
	 		datasetName = "time";
	 		dataspaceId = H5Screate_simple(1, &histLength, NULL);

	 		datasetId = H5Dcreate(groupId, datasetName.c_str(), scalarType, dataspaceId,
	 		 					H5P_DEFAULT);

	 		Scalar* data = new Scalar[histLength];

	 		for(int i=0; i<histLength; i++)
	 			data[i] = number[0]->get_time_array()[i];

	 		status = H5Dwrite(datasetId, scalarType, H5S_ALL, H5S_ALL,
	 		 				    H5P_DEFAULT, data); //time vector made longer by Matteo
	 		delete[] data;
	 		H5Sclose(dataspaceId);
	 		H5Dclose(datasetId);

	//create fields for different diagnostics. histories are 2D matrices, first index specie second index time vector


	 		//ke_species

	 		histgroupId = H5Gcreate(groupId,"KE_tot",0);

	 		for(int isp=0;isp<number_of_species;isp++) {
	 		// loop over all species, dumping for each
	 			  string nm = theSpecies[isp].name;

	 		 		datasetName = nm + "_KE_tot";
	 		 		dataspaceId = H5Screate_simple(1, &histLength, NULL);

	 		 		datasetId = H5Dcreate(histgroupId, datasetName.c_str(), scalarType, dataspaceId,
	 		 		 					H5P_DEFAULT);

	 		 		Scalar* data = new Scalar[histLength];

	 		 		for(int i=0; i<histLength; i++)
	 		 			data[i] = ke_species[isp]->get_data()[i];

	 		 		status = H5Dwrite(datasetId, scalarType, H5S_ALL, H5S_ALL,
	 		 		 				    H5P_DEFAULT, data);
	 		 		delete[] data;
	 		 		H5Sclose(dataspaceId);
	 		 		H5Dclose(datasetId);

	 			  }

	 		 status = H5Sclose(histgroupId);

	//total_density

		 		histgroupId = H5Gcreate(groupId,"n_tot",0);

		 		for(int isp=0;isp<number_of_species;isp++) {
		 		// loop over all species, dumping for each
		 			  string nm = theSpecies[isp].name;

		 		 		datasetName = nm + "_n_tot";
		 		 		dataspaceId = H5Screate_simple(1, &histLength, NULL);

		 		 		datasetId = H5Dcreate(histgroupId, datasetName.c_str(), scalarType, dataspaceId,
		 		 		 					H5P_DEFAULT);

		 		 		Scalar* data = new Scalar[histLength];

		 		 		for(int i=0; i<histLength; i++)
		 		 			data[i] = total_density[isp]->get_data()[i];

		 		 		status = H5Dwrite(datasetId, scalarType, H5S_ALL, H5S_ALL,
		 		 		 				    H5P_DEFAULT, data);
		 		 		delete[] data;
		 		 		H5Sclose(dataspaceId);
		 		 		H5Dclose(datasetId);

		 			  }

		 		 status = H5Sclose(histgroupId);

	//Ave_KE

			 		histgroupId = H5Gcreate(groupId,"KE_ave",0);

			 		for(int isp=0;isp<number_of_species;isp++) {
			 		// loop over all species, dumping for each
			 			  string nm = theSpecies[isp].name;

			 		 		datasetName = nm + "_KE_ave";
			 		 		dataspaceId = H5Screate_simple(1, &histLength, NULL);

			 		 		datasetId = H5Dcreate(histgroupId, datasetName.c_str(), scalarType, dataspaceId,
			 		 		 					H5P_DEFAULT);

			 		 		Scalar* data = new Scalar[histLength];

			 		 		for(int i=0; i<histLength; i++)
			 		 			data[i] = Ave_KE[isp]->get_data()[i];

			 		 		status = H5Dwrite(datasetId, scalarType, H5S_ALL, H5S_ALL,
			 		 		 				    H5P_DEFAULT, data);
			 		 		delete[] data;
			 		 		H5Sclose(dataspaceId);
			 		 		H5Dclose(datasetId);

			 			  }

			 		 status = H5Sclose(histgroupId);

	//create fields for boundary diagnostics

				 		histgroupId = H5Gcreate(groupId,"Iboundary",0);

			 		 	oopicListIter<Ihistdiag> nextdiag(*BoundaryIhist);

			 			if(!BoundaryIhist->isEmpty())
			 			{
			 				nextdiag.restart();

			 				int length = *nextdiag.current()->Ihist->get_hi();

			 				hsize_t histLength;//length;
			 				histLength = length;

			 				datasetName = "boundary_time";
			 				dataspaceId = H5Screate_simple(1, &histLength, NULL);

			 				datasetId = H5Dcreate(histgroupId, datasetName.c_str(), scalarType, dataspaceId, H5P_DEFAULT);

			 				Scalar* data = new Scalar[histLength];

			 				for(int i=0; i<histLength; i++)
			 					data[i] = nextdiag.current()->Ihist->get_time_array()[i];

			 			    status = H5Dwrite(datasetId, scalarType, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
			 			    delete[] data;
			 				H5Sclose(dataspaceId);
			 				H5Dclose(datasetId);
			 			}

			 			for(nextdiag.restart();!nextdiag.Done();nextdiag++) {

					 		bgroupId = H5Gcreate(histgroupId,nextdiag.current()->name,0);

					 		datasetName = "Ihist";
					 		dataspaceId = H5Screate_simple(1, &histLength, NULL);

					 		datasetId = H5Dcreate(bgroupId, datasetName.c_str(), scalarType, dataspaceId, H5P_DEFAULT);

					 		Scalar* data = new Scalar[histLength];

					 		for(int i=0; i<histLength; i++)
					 			data[i] = nextdiag.current()->Ihist->get_data()[i];

					 		status = H5Dwrite(datasetId, scalarType, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
					 		delete[] data;
					 		H5Sclose(dataspaceId);
					 		H5Dclose(datasetId);

			 				  for (int isp=0; isp<number_of_species; isp++) {
					 			  string nm = theSpecies[isp].name;

					 			 datasetName = "Ihist_" + nm;
					 			 dataspaceId = H5Screate_simple(1, &histLength, NULL);

					 			 datasetId = H5Dcreate(bgroupId, datasetName.c_str(), scalarType, dataspaceId, H5P_DEFAULT);

					 			 Scalar* data = new Scalar[histLength];

			 					 for(int i=0; i<histLength; i++)
			 					 	data[i] = nextdiag.current()->Ihist_sp[isp]->get_data()[i];

			 					status = H5Dwrite(datasetId, scalarType, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
			 					delete[] data;
			 					H5Sclose(dataspaceId);
			 					H5Dclose(datasetId);
			 				  }
						 		 status = H5Sclose(bgroupId);
			 				}

				 		 status = H5Sclose(histgroupId);

	//create fields for custom diagnostics

				 		histgroupId = H5Gcreate(groupId,"Custom",0);

				 		DiagList *dlist = theSpace->getDiagList();
				 			oopicListIter<Diag> nextd(*dlist);

				 			if(!dlist->isEmpty())
				 			{
				 				nextd.restart();

				 				int length = *nextd.current()->getHistory()->get_hi();

				 				hsize_t histLength;//length;
				 				histLength = length;

				 				datasetName = "custom_time";
				 				dataspaceId = H5Screate_simple(1, &histLength, NULL);

				 				datasetId = H5Dcreate(histgroupId, datasetName.c_str(), scalarType, dataspaceId, H5P_DEFAULT);

				 				Scalar* data = new Scalar[histLength];

				 				for(int i=0; i<histLength; i++)
				 					data[i] = nextd.current()->getHistory()->get_time_array()[i];

				 			    status = H5Dwrite(datasetId, scalarType, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
				 			    delete[] data;
				 				H5Sclose(dataspaceId);
				 				H5Dclose(datasetId);
				 			}
				 			// ------Updating newDiag-----//
				 			for(nextd.restart();!nextd.Done(); nextd++)
				 			{
				 				datasetName = nextd.current()->getTitle();
				 				dataspaceId = H5Screate_simple(1, &histLength, NULL);

				 				datasetId = H5Dcreate(histgroupId, datasetName.c_str(), scalarType, dataspaceId, H5P_DEFAULT);

				 				Scalar* data = new Scalar[histLength];

				 				for(int i=0; i<histLength; i++)
				 					data[i] = nextd.current()->getHistory()->get_data()[i];

				 			    status = H5Dwrite(datasetId, scalarType, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
				 			    delete[] data;
				 				H5Sclose(dataspaceId);
				 				H5Dclose(datasetId);

				 			}

					 		 status = H5Sclose(histgroupId);

	//close everything


	 status = H5Sclose(groupId);
	 status = H5Dclose(fileId);

	return;
}

#endif