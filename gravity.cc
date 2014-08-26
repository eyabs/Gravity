/*
Gravity Simulator
gravity.cc
member function definition file

Author: Eric Y.  Date: 25 Aug 2014
*/

#include <cmath>
#include <cstdlib>
#include <iostream>
#include "gravity.h"

using namespace std;

///Default constructor
/**ParticleSystem()
 **uses predefined variables when constructing object
 **/
ParticleSystem::ParticleSystem(){
    srand48( time(NULL) );
    setNumDims(2);
    setNumParticles(500);
    setdT(0.1);
    setTmax(100.0);
    
    double lowerLims[2] = {-10,-10};
    double upperLims[2] = {10,10};
    
    iLLims = new double[iNumDims];
    iULims = new double[iNumDims];
    setLLims( lowerLims );
    setULims( upperLims );
    
    createBins();
    
    clearMatrix( iVelocities );
    // cout << "\n" << iNumParticles << "\n" << iNumDims << endl;
    initPositions();
    setInit(1);
}

///Destructor
/**~ParticleSystem()
 **Frees dynamically allocated memory
 **/
ParticleSystem::~ParticleSystem(){
    delete [] iLLims;
    delete [] iULims;

    for( unsigned long particle = 0; particle < iNumParticles; ++particle){
        delete [] iGravities[particle];
        delete [] iVelocities[particle];
        delete [] iPositions[particle];
    }

    delete [] iGravities; 
    delete [] iVelocities;
    delete [] iPositions;
}

///Member Function createBins()
/**Dynamically allocated matricies to hole program data
 **i.e. gravity, acceleration, velocities and positions
 **/
void ParticleSystem::createBins(){
    iGravities = new double*[iNumParticles];
        for( unsigned long particle = 0; particle < iNumParticles; particle++)
            iGravities[particle] = new double[iNumDims]; 
    iVelocities = new double*[iNumParticles];
        for( unsigned long particle = 0; particle < iNumParticles; particle++)
            iVelocities[particle] = new double[iNumDims]; 
    iPositions = new double*[iNumParticles];
    	for( unsigned long particle = 0; particle < iNumParticles; ++particle)
        	iPositions[particle] = new double[iNumDims]; 
} 
  
///Member Function clearMatrix( double **)
/**Clears contents of input array
 **/
void ParticleSystem::clearMatrix( double **aMatrix ){
	for( unsigned long rowLoop = 0; rowLoop < iNumParticles; rowLoop++ ){
		for( unsigned long colLoop = 0; colLoop < iNumDims; colLoop++ )
			aMatrix[rowLoop][colLoop] = 0.0;
	}
}

///Member Function setLLims( douuble *aLLims)
/**Fills the internal array with the lower bounds of plotting from input array
 **/
void ParticleSystem::setLLims( double *aLLims ){
    for( unsigned long dim = 0; dim < iNumDims; ++dim )
        iLLims[dim] = aLLims[dim];
}

///Member Function setULims( douuble *aLLims)
/**Fills the internal array with the upper bounds of plotting from input array
 **/
void ParticleSystem::setULims( double *aULims ){
    for( unsigned long dim = 0; dim < iNumDims; ++dim )
        iULims[dim] = aULims[dim];
}

///Member Function setLLims( douuble *aLLims)
/**Fills the internal array with the lower bounds of plotting from input array
 **/
void ParticleSystem::initPositions(){
    for(unsigned long particle = 0; particle < iNumParticles; ++particle){
        for(unsigned long dim = 0; dim < iNumDims; ++dim)             
            iPositions[particle][dim] = ( iULims[dim]-iLLims[dim] ) * drand48() + iLLims[dim];  
    }
}

///Member Function dumpState()
/**Prints all the particle locations for the current state
 **/
void ParticleSystem::dumpState(){
     for(unsigned long particle = 0; particle < iNumParticles; ++particle){
        for(unsigned long dim = 0; dim < iNumDims; ++dim)             
             cout << iPositions[particle][dim] << " ";
        cout << endl;
    }
}

///Member Function gravForce(double *, double *)
/**Calculates a gravittional force component between to particles w/ input position arrays
 ** F = G * m * m / r ^ (3) * <displacement vector component>
 **/
double ParticleSystem::gravForce( unsigned long aDim, /*unsigned long numDims,*/ double *thisParticle, double *otherParticle ){
    double gForce;
    double distance = 0;
    
	for(unsigned long dim = 0; dim < iNumDims; dim++)
        distance += pow( (otherParticle[dim] - thisParticle[dim]), 2);
    		
	distance = pow(distance, 0.5);
	
	//
	//Avoiding a divide by zero error
	if( pow(distance, -2.0) > 1E2 ){
		if(otherParticle[aDim] > thisParticle[aDim]) return 1;
		else return -1;
	
	}
    gForce = gConst / ( pow(distance, 3) ) * (otherParticle[aDim] - thisParticle[aDim]);

    return gForce;
}

///Member Function findGravity()
/**using the gravForce function this function fills the matrix of 
 ** Grav. Acceleration components for each particle
 **/
void ParticleSystem::findGravity(){
    
    clearMatrix( iGravities );
    
    //for each particle
    for(unsigned long thisParticle = 0; thisParticle < iNumParticles; thisParticle++){
   		
        //find the gravity from each other particle and make a running sum  
        for(unsigned long otherParticle = 0; otherParticle < iNumParticles; otherParticle++){      	
        	
        	for(unsigned long thisDim = 0; thisDim < iNumDims; thisDim++){	
        			
        		if( otherParticle != thisParticle ){  	        			
        				iGravities[thisParticle][thisDim] +=
        					gravForce(thisDim, iPositions[thisParticle], iPositions[otherParticle]);
        		}
        		
        	}		
        }
    }
}

///Member Function stepVelocities()
/**for each particle in each dimension
 ** Change the velocities by dv = g * dt	
 **/
void ParticleSystem::stepVelocities(){
	//for each particle
    for(unsigned long particle = 0; particle < iNumParticles; particle++){
   	    //in each dimension
   	    for(unsigned long dim = 0; dim < iNumDims; dim++){
   	    	//Change the velocities by dv = g * dt	
        	iVelocities[particle][dim] += iGravities[particle][dim] * idT;        	        		
        }		
    }
}

///Member Function stepVelocities()
/**for each particle in each dimension
 ** Change the position by dx = v * dt	
 **/
void ParticleSystem::stepPositions(){
	//for each particle
    for(unsigned long particle = 0; particle < iNumParticles; particle++){
   	    //in each dimension
   	    for(unsigned long dim = 0; dim < iNumDims; dim++){
   	    	//Change the position by dx = v * dt	
        	 iPositions[particle][dim] += iVelocities[particle][dim] * idT; 		
        }
    }
}



























