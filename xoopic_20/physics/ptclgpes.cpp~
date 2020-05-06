/*
====================================================================

PARTICLEGROUPES.CPP

Electrostatic version of ParticleGroup.

Revision History
0.99	(JohnV 01-23-95)	Original version.
0.991	(JohnV 05-28-95)	Complete and integrate.

Should optimize for B=0.

====================================================================
*/

#include "ptclgpes.h"
#include "fields.h"
#include	"particle.h"

void ParticleGroupES::advance(Fields& fields, Scalar dt)
{
	Scalar	f = 0.5*dt*q_over_m;
	Scalar	q_dt = q/dt;
	Vector3	uPrime;
	Vector3 uOld;
	Vector3	a;
	Vector3	t;
	Vector3	s;
	Vector2* x = ParticleGroup::x;
	Vector3* u = ParticleGroup::u;
	Vector3	dxMKS;
	Boundary*	bPtr;
	Grid&	grid = *fields.get_grid();
	Vector3	B;
//	register Scalar temp;
	double localEnergy = 0; // overflows if float

	for (int i=0; i<n; i++, x++, u++)
	{
		uOld = *u;
		bPtr = 0;
		a = f*fields.E(*x);
		*u += a;									//	half acceleration
		B = fields.B(*x);   // speed-ups make an array for f*E and f*B
		                    // have a and t returned with one call
		                    // if test B.isNonZero outside this loop
		                    // then f*B doesn't have to be returned and
								  // an if test and += are saved per particle
		if (B.isNonZero())
		{
			t = f*B;								// reuse B!
			uPrime = *u + u->cross(t);
			s = 2*t/(1 + t*t);
			*u += uPrime.cross(s);			//	rotation
		}
		*u += a;									//	half acceleration
		dxMKS = grid.differentialMove(fields.getMKS(*x), *u, *u, dt);
		while (fabs(dxMKS.e1()) + fabs(dxMKS.e2()) > 1E-25 && !bPtr)
			bPtr = grid.translate(*x, dxMKS);
		if (bPtr)
		{
		  //	send this particle to boundary
		  bPtr->collect(*(new Particle(*x, *u, species, get_np2c(i),
					       (BOOL)(qArray!=0))), dxMKS);
		  n--; // decrement number of particles
		  //	Move last particle into this slot and advance it.
		  if (i == n) break; //	last particle in array?
		  *x = ParticleGroup::x[n]; //	move last particle into open slot
		  *u = ParticleGroup::u[n];
		  if (qArray) qArray[i] = qArray[n];
		  i--; // Set index to do i again.
		  x--; // Setup to redo the array element
		  u--; // it gets inc in for loop
		}
		else{
		  if (qArray)
		    localEnergy += qArray[i]**u*uOld;
		  else
		    localEnergy += *u*uOld;
			/*Scalar uOld_x = uOld.e1();
			Scalar uOld_y = uOld.e2();
			Scalar uOld_z = uOld.e3();
			Scalar uOld_norm = uOld.magnitude (); 

			Vector3 unew = *u;
			Scalar unew_x = unew.e1();
			Scalar unew_y = unew.e2();
			Scalar unew_z = unew.e3();
			Scalar unew_norm = unew.magnitude (); 

			Scalar diff_percent = (unew_norm - uOld_norm) / uOld_norm *100;*/
		}
#ifdef DEBUG
		if(strncmp(HOSTTYPE, "alpha", 5) != 0)
		  {
		    if(((*x).e1()>grid.getJ() || (*x).e1()<0 || (*x).e2()>grid.getK()
			|| (*x).e2()<0)&&!bPtr) {
		      printf("Particle escaped.\n");
		      int *crash = 0;	
		      *crash = 10;
		    }
		  }
		
#endif //DEBUG
		
	}

	if (qArray)
	  species->addKineticEnergy(0.5*localEnergy/q_over_m);
	else
	  species->addKineticEnergy(0.5*localEnergy*get_m());
}


	


void omp_advance(ParticleGroup *group, Vector3** ENode_local, Vector3** BNode_local , Grid&	grid , Vector2** X_local  , Scalar dt , Boundary*** bc1_local   , Boundary*** bc2_local, double * grp_KE_qTrue, double * grp_KE_qFalse)
{
	register Scalar  q = group->get_q();  //charge * np2c
	register Scalar* qArray = group->get_qArray();
	register Scalar	f = 0.5*dt*group->get_q_over_m ();  
	register Scalar	q_dt = q/dt;
	register Vector3	uPrime;
	register Vector3 uOld;
	register Vector3	a;
	register Vector3	t;
	register Vector3	s;
	register Vector2* x = group->get_x();
	register Vector3* u = group->get_u();
	register Vector3	dxMKS;
	Boundary*	bPtr;
//	Grid&	grid = *fields.get_grid();
	Vector3	B;
//	register Scalar temp;
	double localEnergy = 0; // overflows if float
	int n=group->get_n ();
	
	for (int i=0; i<n; i++, x++, u++)
	{
		uOld = *u;
		bPtr = 0;
		//a = f * interpolateBilinear (ENode_local,*x);
		a = f*grid.interpolateBilinear (ENode_local,*x);
		//a = f*fields.E(*x);
		*u += a;									//	half acceleration
		//B = interpolateBilinear (BNode_local,*x);
		B = grid.interpolateBilinear (BNode_local,*x);
		//B = fields.B(*x);   // speed-ups make an array for f*E and f*B
		                    // have a and t returned with one call
		                    // if test B.isNonZero outside this loop
		                    // then f*B doesn't have to be returned and
								  // an if test and += are saved per particle
		if (B.isNonZero())
		{
			t = f*B;								// reuse B!
			uPrime = *u + u->cross(t);
			s = 2*t/(1 + t*t);
			*u += uPrime.cross(s);			//	rotation
		}
		*u += a;									//	half acceleration

		dxMKS = grid.differentialMove(grid.getMKS_local(*x, X_local), *u, *u, dt);
		//dxMKS = Vector3(u->e1()*dt, u->e2()*dt, u->e3()*dt); // increases speed by 25% if replaced !!!!

		
		while (fabs(dxMKS.e1()) + fabs(dxMKS.e2()) > 1E-25 && !bPtr) //translate incremently until done or crossed a boundary
			bPtr = grid.translate_omp(*x, dxMKS , bc1_local, bc2_local);
			//bPtr = grid.translate(*x, dxMKS);			
		if (bPtr) // crossed a boundary
		{
		//	send this particle to boundary
		#pragma omp critical  
			{
				bPtr->collect(*(new Particle(*x, *u, group->get_species(), group->get_np2c(i),(BOOL)(qArray!=0))), dxMKS);
			} // bPtr.normal
		  group->set_n(n-1); // decrement number of particles
		  n--;
		  //	Move last particle into this slot and advance it.
		  if (i == n) break; //	last particle in array?
		  *x = group->get_x(n); //	move last particle into open slot
		  *u = group->get_u(n);
		  if (qArray) qArray[i] = qArray[n];
		  i--; // Set index to do i again.
		  x--; // Setup to redo the array element
		  u--; // it gets inc in for loop
		}
		else{ // didn't cross a boundary
		  if (qArray)
		    localEnergy += qArray[i]**u*uOld;
		  else
		    localEnergy += *u*uOld;
			/*Scalar uOld_x = uOld.e1();
			Scalar uOld_y = uOld.e2();
			Scalar uOld_z = uOld.e3();
			Scalar uOld_norm = uOld.magnitude (); 

			Vector3 unew = *u;
			Scalar unew_x = unew.e1();
			Scalar unew_y = unew.e2();
			Scalar unew_z = unew.e3();
			Scalar unew_norm = unew.magnitude (); 

			Scalar diff_percent = (unew_norm - uOld_norm) / uOld_norm *100;*/
		}
 #ifdef DEBUG
		if(strncmp(HOSTTYPE, "alpha", 5) != 0)
		  {
		    if(((*x).e1()>grid.getJ() || (*x).e1()<0 || (*x).e2()>grid.getK()
			|| (*x).e2()<0)&&!bPtr) {
		      printf("Particle escaped.\n");
		      int *crash = 0;	
		      *crash = 10;
		    }
		  }
		
#endif //DEBUG
		
	} // end loop on particles 

		
		if (qArray){
			//group->get_species()->addKineticEnergy(0.5*localEnergy/group->get_q_over_m());
			*grp_KE_qTrue+=0.5*localEnergy/group->get_q_over_m();
		}
		else{
			//group->get_species()->addKineticEnergy(0.5*localEnergy*group->get_m());
			*grp_KE_qFalse+=0.5*localEnergy*group->get_m();
		}
	
	}
		
