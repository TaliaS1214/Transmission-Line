//
//  pyWriter.cpp
//  TransmissionLine
//
//  Created by Sean Talia on 8/27/14.
//  Copyright (c) 2014 Sean Talia. All rights reserved.
//

#include <iostream>
#include <fstream>
#include "pyWriter.h"

using namespace std;


void writeArray(double * cppArray,
                string   pyArray,
                int      n,
                ofstream &pyFile);



void pyWriter(double *x_vals,
              double **y_vals,
              int          N,
              int    point_count) {
    
    ofstream pyFile;
    
    string indent = "    ";
    
    pyFile.open("/Users/SeanTalia/Coding_Files/Thesis_Files/2_TransLine/transLinePlot.py");
    pyFile << "import numpy as np"                                                                                       << endl;
    pyFile << "import matplotlib.pyplot as plt"                                                                          << endl;
    pyFile << "def main():"                                                                                              << endl;
    pyFile                                                                                                               << endl;
    pyFile << indent << "current_matrix = np.zeros((" << N << ", " << point_count << "))"                                << endl;

    writeArray(x_vals, "x_vals", point_count, pyFile);
    
    string current_vector_name;
    
    for (int i=0; i < N; i++) {
        string k = to_string(i);
        current_vector_name = "current_matrix[" + k + "]";
        writeArray(y_vals[i], current_vector_name, point_count, pyFile);
    }
    
    pyFile                                                                                                               << endl;
    pyFile << indent << "fig = plt.figure()"                                                                             << endl;
    pyFile << indent << "ax = fig.add_subplot(1, 1, 1)"                                                                  << endl;
    pyFile                                                                                                               << endl;
    pyFile << indent << "# Labeling Axes"                                                                                << endl;
    pyFile << indent << "ax.set_xlabel('Time')"                                                                          << endl;
    pyFile << indent << "ax.set_ylabel('Current')"                                                                       << endl;
    pyFile                                                                                                               << endl;
    pyFile << indent << "# Curve Colors"                                                                                 << endl;
    pyFile << indent << "colors = ['#ff1493', '#ff0000', '#ff4500', '#ff6347', '#ff7f50',"                               << endl;
    pyFile << indent << indent << "'#ffa500', '#ffb600', '#ffd700', '#fff000', '#ffff00',"                               << endl;
    pyFile << indent << indent << "'#9acd32', '#adff2f', '#00ff00', '#00cc44', '#00aa22',"                               << endl;
    pyFile << indent << indent << "'#228b22', '#7fff7f', '#4cff4c', '#00ffff', '#00fa9a',"                               << endl;
    pyFile << indent << indent << "'#00bfff', '#4169e1', '#0000cd', '#000080', '#191970',"                               << endl;
    pyFile << indent << indent << "'#5518ab', '#4c177d', '#44146f', '#3b1261', '#330f53',"                               << endl;
    pyFile << indent << indent << "'#2a0d45', '#220a37', '#190729', '#11051b', '#82020d']"                               << endl;
    pyFile                                                                                                               << endl;
    pyFile << indent << "#Plotting our curves"                                                                           << endl;
    pyFile << indent << "for i in range(" << N << "):"                                                                   << endl;
    pyFile << indent << indent << "ax.plot(x_vals, current_matrix[i], color = colors[i])"                                << endl;
    pyFile                                                                                                               << endl;
    pyFile << indent << "# Show the Plot!"                                                                               << endl;
    pyFile << indent << "plt.show()"                                                                                     << endl;
    pyFile                                                                                                               << endl;
    pyFile << "main()"                                                                                                   << endl;
    pyFile                                                                                                               << endl;
    pyFile.close();
}

void writeArray(double *cppArray,
                string   pyArray,
                int      n,
                ofstream &pyFile) {
    
    string indent = "    ";
    
    pyFile << indent << pyArray << " = np.array(\\" << endl;
    pyFile << indent << "[";
    for(int i=0; i < n; i++){
        if(i < n-1) {
            pyFile << cppArray[i] << ", ";
        }
        else {
            pyFile << cppArray[i] << "])" << endl;
        }
    }
}