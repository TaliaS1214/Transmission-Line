//
//  transLine.cpp
//  TransmissionLine
//
//  Created by Sean Talia on 5/28/14.
//  Copyright (c) 2014 Sean Talia. All rights reserved.
//

#include <math.h>
#include <vector>
#include "transLine.h"

using namespace std;

static const double A_0 = 5.0;            // Voltage source amplitude
static const double omega = 0.5;          // Voltage source frequency
static const double R = 0.2;              // Total line resistance
static const double L = 2.5;              // Total line inductance
static const double C = 4.5;              // Total line capacitance
static const double R_load = 10.0;        // Load resistance
static const double L_load = 15.0;        // Load inductance

double V(double t){

    return A_0 * cos(omega * t);

}

// Note that N here is NOT the size of the array y, but rather floor(array_size / 2)
// So if, for example, y is an array of length 11, then N should be 5.

void transLine(double t, int N, double *y, double *f){

    // Setting V'_1 through V'_N
    for(int i=0; i < N; i++){
        f[i] = (N / C) * (y[N+i] - y[N+i+1]);
    }

    // Setting I'_1
    f[N] = (N / L) * (V(t) - y[0] - (R / N) * y[N]);

    // Setting I'_2 though I'_N
    for(int i=0; i < N-1; i++){
        f[N+i+1] = (N / L) * (y[i] - y[i+1] - (R / N) * y[N+i+1]);
    }

    // Setting I'_load
    f[2*N] = (1 / L_load) * (y[N-1] - R_load * y[2*N]);

}
