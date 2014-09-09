//
//  tlODEsolver.cpp
//  TransmissionLine
//
//  Created by Sean Talia on 5/28/14.
//  Copyright (c) 2014 Sean Talia. All rights reserved.
//

#include <iostream>
#include <vector>
#include <math.h>
#include "tlODEsolver.h"
#include "transLine.h"
#include "pyWriter.h"

using namespace std;

static void init(double *k, int D);
static void displayResults(double *v, int D);
static void quickShow(double *v, int D);

void tl_ODE_solver(
    void    (*g)(double, int, double*, double*),                              // ODE function
    double  h,                                                                // Step size
    double  T,                                                                // Final Time
    double  *y0,                                                              // Initial Data
    int     M                                                                 // Dimension of System
)
{
// ------------------- Setting up our Necessary Variables------------------------------

    void (*f)(double, int, double*, double*) = g;                             // The function to use

    double t = 0;                                                             // Initial time
    int steps = ceil(T / h);                                                  // Number of time steps

    int N = floor(M / 2);                                                     // This is the number of transmission line
                                                                              // segments, which we use as input for the function
                                                                              // pointer

    double *y = y0;                                                           // Initializing y

    double *k1 = new double[M]; double *k1_tmp = new double[M];               // Vectors needed for RK4 method
    double *k2 = new double[M]; double *k2_tmp = new double[M];
    double *k3 = new double[M]; double *k3_tmp = new double[M];
    double *k4 = new double[M];

    init(k1, M); init(k1_tmp, M);                                             // Initializing the temp vectors
    init(k2, M); init(k2_tmp, M);
    init(k3, M); init(k3_tmp, M);
    init(k4, M);

// ------------------- Variables for Output to NumPy Array --------------------

    // In this section, we'll need to create k arrays, where k is the number of
    // points we want to plot. Each array should be of length M.
    // As we move forward in time via the numerical method, we'll need to periodically
    // record the values stored in *y, so that we can track the values of
    // the current and voltages through each segment in time.

    int point_count = 351;                                                  // Number of points to be plotted in Python graph
    double x_py[point_count];                                               // Set of x values to output to Python array
    double *y_py[N];                                                        // Creating as many arrays as points we plan to plot.
                                                                            // Each array will contain M values (which we record from y*
                                                                            // while the numerical method executes)

    for (int i = 0; i < N; ++i) {                                           // Setting each entry of the pointer array to be a
        y_py[i] = new double[point_count];                                  // static array of length however many points we plot.
    }

    init(x_py, point_count);

    for(int i=0; i < N; i++){
        init(y_py[i], point_count);
    }

    double spacing = T / (point_count - 1);                                 // Spacing the points equally (pointCount - 1 because we want
                                                                            // to include the endpoints 0 and T as x-values

    for(int j = 0; j < point_count; j++){
        x_py[j] = double(j * spacing);
    }

    int inner_count = 0; int outer_count = 0;                               // Here we're seting our initial counters to 0
    int inner_count_max = spacing / h;                                      // Specifying the threshold we'll reach before we reset
                                                                            // the inner counter (inner count goes up by one after every
                                                                            // step in time, outer count goes up by one every time we
                                                                            // reach a time at which we want to record data)

    for (int i = 0; i < N; ++i){                                            // Inserting our initial data as the first vector in the matrix
        y_py[0][i] = y[N+i];
    }

    outer_count += 1;


// ---------------------Executing the Numerical Method ---------------------

    for(int i = 0; i <= steps; ++i){                                        // Time step loop

        f(t, N, y, k1);                                                     // Setting k1 = f(t, y)

        for(int j = 0; j < M; ++j) {
            k1_tmp[j] = y[j] + (h / 2.0) * k1[j];
        }

        f(t + h/2.0, N, k1_tmp, k2);                                        // Setting k2 = f(t + h/2, y + (h/2) * k1)

        for(int j = 0; j < M; ++j) {
            k2_tmp[j] = y[j] + (h / 2.0) * k2[j];
        }

        f(t + h/2.0, N, k2_tmp, k3);                                        // Setting k3 = f(t + h/2, y + (h/2) * k2)

        for(int j = 0; j < M; ++j) {
            k3_tmp[j] = y[j] + h * k3[j];
        }

        f(t + h, N, k3_tmp, k4);                                            // Setting k4 = f(t + h, y + h * k3)

        for(int j = 0; j < M; ++j) {
            y[j] += (h / 6) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);      // Updating the value of y
        }

        // Here we construct our Python matrix of y-values to be plotted

        if (inner_count >= inner_count_max && outer_count < point_count) {  // If we've reached a critical x value,
                                                                            // record data and reset the inner count
            for (int i = 0; i < N; ++i){
                y_py[i][outer_count] = y[N+i];
            }

            outer_count += 1;
            inner_count = 0;
        }
        else {
            inner_count += 1;                                               // If we haven't, don't do anything
        }

        t += h;
    }

// --------------------- Clean Up and Display Results -------------------------------

    delete[] k1; delete[] k1_tmp;                                           // Removing the temp vectors from memory
    delete[] k2; delete[] k2_tmp;
    delete[] k3; delete[] k3_tmp;
    delete[] k4;

    //displayResults(y, M);                                                   // Display f(T) in Terminal window

    // Here we'll need to send the data to our Python function, which should plot
    // the vectors. We'll send it our multidimensional array, which Python is to
    // parse through so that it may plot the currents.

    for (int i=0; i < N; i++) {
        y_py[i][point_count-1] = y[N+i];
        //quickShow(y_py[i], point_count);
    }

    pyWriter(x_py, y_py, N, point_count);
}

void displayResults(double *v, int D) {
    cout << "The approximate solution is y(T) = " << '\n';
    for(int i=0; i < D; i++){
        if (i != (D-1)) {
            if (i == 0) {
                cout << "[(Voltages through segments)" << '\n';
            }
            if (i == floor(D/2)) {
                cout << "(Currents through segments)" << '\n';
            }
            cout << v[i] << ", " << '\n';
        }
        else {
            cout << "(Current through load)" << '\n';
            cout << v[i] << "] " << '\n';
        }
    }
}

void quickShow(double *v, int D){
    for(int i=0; i < D; i++) {
        if (i < (D-1)) {
            if (i == 0) { cout << "["; }
            cout << "(" << i << ")" << v[i] << ", ";

        }
        else { cout << "(" << i << ")" << v[i] << "]" << '\n'; }
    }
    cout << "---------------------------------" << '\n';
}

void init(double *k, int D) {
    for(int i = 0; i < D; i++)
        k[i] = 0.0;
}
