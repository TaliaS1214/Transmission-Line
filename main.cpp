//
//  main.cpp
//  TransmissionLine
//
//  Created by Sean Talia on 5/28/14.
//  Copyright (c) 2014 Sean Talia. All rights reserved.
//

#include <math.h>
#include <iostream>
#include "transLine.h"
#include "tlODEsolver.h"

using namespace std;

int main(int argc, const char * argv[])
{
    // Function to Use
    void (*f)(double, int, double *, double *) = transLine;

    // Step Size
    double h = pow(10, -6);

    // Final Time
    double T = 80.0;

    // Number of segments, size of system
    int N = 15;
    int M = 2*N+1;

    // Initial Data (must be an array of length M)
    double y0[M];
    for (int i = 0; i < M; i++){
        y0[i] = 0.0;
    }

    tl_ODE_solver(f, h, T, y0, M);

}
