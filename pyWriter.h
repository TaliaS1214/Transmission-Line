//
//  pyWriter.h
//  TransmissionLine
//
//  Created by Sean Talia on 8/27/14.
//  Copyright (c) 2014 Sean Talia. All rights reserved.
//

#ifndef __TransmissionLine__pyWriter__
#define __TransmissionLine__pyWriter__

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

void writeArray(double * cppArray, string pyArray, int n, ofstream &pyFile);

void pyWriter(double *x, double **y, int N, int pointCount);


#endif /* defined(__TransmissionLine__pyWriter__) */
