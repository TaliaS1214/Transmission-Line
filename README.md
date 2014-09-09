Transmission Line ODE Solver
==========================================
The files included here are a combination of C++ and Python files. The ODE solver can actually be used to solve any system of 1st order
ordinary differential equations, but I only used it to model a lumped-element transmission line divided into N segments. The function(s) accept 11 different parameters and solve for the currents flowing through each segment. All of the heavy lifting is 
done by the C++ files, which 1) solve a differential equation via a numerical method and 2) write a Python file (transLinePlot.py) that plots the solutions. In order to use this software, just follow these instructions:

1. Open the transLine.cpp file and manipulate the parameters A_0, omega, R, L, C, R_load, and L_load.
2. Open the main.cpp file and manipulate the parameters h T, N, and y0.
3. In pyWriter.cpp, change line 32 to "pyFile.open("[PATHNAME]/transLinePlot.py")",
   where PATHNAME is the path to the current directory .
4. In your terminal, cd to the current directory and run the command "make run".

That's it! 


