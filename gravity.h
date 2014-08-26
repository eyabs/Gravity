/*
Gravity Simulator
gravity.h - object definition file

Author: Eric Y.  Date: 25 Aug 2014
*/

#include <iostream>
#include <cstdlib>
#include <cmath>
#define gConst 1.0E-2

using namespace std;

class ParticleSystem{
  //
  //private internal variables
  private:
    bool iIsInitialized;
    unsigned long iNumDims;
    unsigned long iNumParticles;
    double idT;
    double iTmax;
    double *iLLims;
    double *iULims;
    double **iGravities;
    double **iVelocities;
    double **iPositions;

  public:
    
    //
    //Constructors
    ParticleSystem();
    ParticleSystem( unsigned long aNumDims, unsigned long aNumParticles,
                    double adT, double *aLLims, double *aULims );
    //
    //

    //
    //public setter functions
    void setInit( bool aIsInitialized ) { iIsInitialized = aIsInitialized; }
    void setNumDims( unsigned long aNumDims ) { iNumDims = aNumDims; }
    void setNumParticles( unsigned long aNumParticles ) { iNumParticles = aNumParticles; }
    void setdT( double adT ) { idT = adT; }
    void setTmax( double aTmax ) { iTmax = aTmax; }
    void setLLims( double *aLLims );
    void setULims( double *aULims );
    void setGravs( double **aGravities );
    void setVels( double **aVels );
    void setPosits( double **aPositions );

    //
    //public getter functions
    bool initialized() { return iIsInitialized; }
    unsigned long numDims() { return iNumDims; }
    unsigned long numParticles() { return iNumParticles; }
    double dT() { return idT; }
    double lowerLimit( unsigned long aDim ) { return iLLims[aDim]; }
    double upperLimit( unsigned long aDim ) { return iULims[aDim]; }
    double gravity( unsigned long aParticleID, unsigned long aDim )
        { return iGravities[aParticleID][aDim]; }
    double velocity( unsigned long aParticleID, unsigned long aDim )
        { return iVelocities[aParticleID][aDim]; }
    double position( unsigned long aParticleID, unsigned long aDim )
        { return iPositions[aParticleID][aDim]; }

    //
    //Member Functions
    void initPositions();
    void createBins();
	void clearMatrix( double ** );
	void dumpState(); 
	double gravForce(unsigned long, double* , double*);
	void findGravity();
	void stepVelocities();
	void stepPositions();
    
	//
    //Destructors
    ~ParticleSystem();
};
























