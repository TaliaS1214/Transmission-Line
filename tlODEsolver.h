//
//  tlODEsolver.h
//  TransmissionLine
//
//  Created by Sean Talia on 5/28/14.
//  Copyright (c) 2014 Sean Talia. All rights reserved.
//

#ifndef __TransmissionLine__tlODEsolver__
#define __TransmissionLine__tlODEsolver__

#include <iostream>
#include <vector>

using namespace std;

void tl_ODE_solver(void (*g)(double, int, double*, double*),    // ODE function
                   double h,                                    // Step size
                   double T,                                    // Final Time
                   double *y0,                                  // Initial Data
                   int D);

#endif /* defined(__TransmissionLine__tlODEsolver__) */
