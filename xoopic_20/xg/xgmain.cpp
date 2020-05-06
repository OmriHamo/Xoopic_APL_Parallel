extern "C++" void Physics(void);
extern "C++" void Maintain_Diagnostics();
extern "C" void XGMainLoop();
extern "C" void Quit(void);
#include <iostream>
#include "ovector.h"
#include "oopiclist.h"
#include "diagn.h"
#include "oops.h"
extern oopicList<Diagnostics> *theDiagnostics;

#include <omp.h>

#ifdef TIMING
extern "C" {
#include <sys/time.h>
}
double Time0=-1;
long int TimeIndex = 0;
#endif

void XGMainLoop()
{

#ifdef TIMING
  struct timeval TheTime;
  double TimeUsed;

  gettimeofday(&TheTime,0);

  // initialize the first time
  if(Time0==-1)  Time0 = 1.0e-6*(double)TheTime.tv_usec + (double)TheTime.tv_sec;

  TimeUsed =  1.0e-6*(double)TheTime.tv_usec + (double)TheTime.tv_sec - Time0;
  
  fprintf(stderr,"\n%ld %f",TimeIndex,(Scalar)TimeUsed);
  TimeIndex++;

#endif //TIMING

  // We have to catch exceptions here, because XGrafix can't do it

//double t0 = omp_get_wtime();
	
  try {
   Physics();  // fields then particles
	  {
	  oopicListIter <Diagnostics> Diter(*theDiagnostics);
	  for(Diter.restart();!Diter.Done();Diter++)
		 Diter.current()->UpdatePreDiagnostics();
	}
	  {
	  oopicListIter <Diagnostics> Diter(*theDiagnostics);
	  for(Diter.restart();!Diter.Done();Diter++)
		 Diter.current()->UpdatePostDiagnostics();
	}
  }
  catch(Oops& oopsPhys){
    std::cout << "\nAn error occurred while processing XGMainLoop(). "
		<< "Back trace follows.\n"
		<< oopsPhys.getMessage() << endl;
    Quit();
    exit(1);
  }

//double t1 = omp_get_wtime()-t0;
//double speed_up_f=1.595/t1;
//printf (" The time step duration is %f [sec] with %d threads \n",t1,1); //
//printf (" The speed up factor is %f with %d threads \n",speed_up_f,1); //

}

void taskYield()
{}
