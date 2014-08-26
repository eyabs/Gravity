/*
Gravity Simulator
gravityMain.cc - main program file wich uses the ParticleSystem Class definitions
defined in gravity.h and gravity.cc

Author: Eric Y.  Date: 25 Aug 2014
*/

#include <cmath>
#include <cstdlib>
#include <iostream>
#include "gravity.h"

using namespace std;


int main(void){

    ParticleSystem MySystem;		//defining the system
    
    //
    //step the system 10000 times
    for( int time = 0; time < 10000; time++){
    	if(time%50==0)					//
    		MySystem.dumpState();		//print every 50th state
    	MySystem.findGravity();
    	MySystem.stepVelocities();
    	MySystem.stepPositions();
    	cout << "\n\n";
    }
       
    return 0;
}
