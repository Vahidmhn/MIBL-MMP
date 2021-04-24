# MIBL-MMP
This is a C++ package to solve Mixed-Integer bi-linear Maximum Multiplicative Problems. This class of optimization problems arises in many applications such as finding Nash bargaining solution (Nash social welfare optimization), capacity allocation market, reliability optimization, etc. Users are wellcome to use it for academic purposes. The following publication can be cited for reference:

....


## Requirements
This package uses CPLEX 12.10 to solve needed subproblems. Therefore, it is necessary to have this solver properly installed and setup on your machine. Remember to compile the code you need to link the following libraries in the _makefile_:
1)ibm/ILOG/CPLEX_Studio1210/cplex/lib/x86-64_linux/static_pic/libilocplex.a 
2)ibm/ILOG/CPLEX_Studio1210/cplex/lib/x86-64_linux/static_pic/libcplex.a 
3)ibm/ILOG/CPLEX_Studio1210/concert/lib/x86-64_linux/static_pic/libconcert.a 
4)-lm 
5)-lpthread 
6)-ldl

For users convenience, a _makefile_ is provided in the package. To use the default _makefile_, please remember to replace [...] by the installation directory of CPLEX in your machine in "makefiles/Makefile-Release.mk" and "makefiles/Makefile-Debug.mk". 

(** The aforementioned explanation is for a linux OS. Please find the proper settings for linking CPLEX in other operating systems otherwise.)

## Setup and how to use
After compling the code, an executive file named "MILBMMP" will be generated.
