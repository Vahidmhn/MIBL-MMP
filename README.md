# MIBL-MMP
This is a C++ package to solve Mixed-Integer bi-linear Maximum Multiplicative Problems. This class of optimization problems arises in many applications such as finding Nash bargaining solution (Nash social welfare optimization), capacity allocation market, reliability optimization, etc. Users are wellcome to use it for academic purposes. The following publication can be cited for reference:

....


## Requirements
This package uses CPLEX 12.10 to solve needed subproblems. Therefore, it is necessary to have this solver properly installed and setup on your machine. Remember to compile the code you need to link the following libraries in the _makefile_:<br/>
1)ibm/ILOG/CPLEX_Studio1210/cplex/lib/x86-64_linux/static_pic/libilocplex.a <br/>
2)ibm/ILOG/CPLEX_Studio1210/cplex/lib/x86-64_linux/static_pic/libcplex.a <br/>
3)ibm/ILOG/CPLEX_Studio1210/concert/lib/x86-64_linux/static_pic/libconcert.a <br/>
4)-lm <br/>
5)-lpthread <br/>
6)-ldl

For users convenience, a _makefile_ is provided in the package. To use the default _makefile_, please remember to replace [...] by the installation directory of CPLEX in your machine in "makefiles/Makefile-Release.mk" and "makefiles/Makefile-Debug.mk". 

(** The aforementioned explanation is for a linux OS. Please find the proper settings for linking CPLEX in other operating systems otherwise.)

## Setup and how to use
After compling the code, an executive file named "MILBMMP" will be generated. To solve your instance you need to run this executive program three following inputs:
1) formatted lp file of the instance
2) number of pieces for McCormick dual finder operation (_k*D_ in the paper)
3) number of pieces for McCormick cut (_k^P_ in the paper)

The only requirement for the lp files for this program is that two first constraints must indicate the variables representing the linear terms. For example: <br/>
"Minimize <br/>
 obj: <br/>
Subject To <br/>
 c1:    x1  = 0<br/>
 c2:    x2  = 0<br/>
 c3:    ... "
 
 Where there has to be a constraint that equates _x1_ to the first linear term and _x2_ to the second linear term of the multiplicative objective function.
 
### Example
An example lp file "1.lp" has been provided in the packege. By runnin the follwoing command in the terminal the program will solve this instance and show the following result:<br/>
"#_of_constraints #_of_variables instance_name<br/>
Lower_bound Upper_bound     'Time = ' solution_time  'y1= ' value_of_first_term  'y2= ' value_of_second_term"<br/>
for example:<br/>
./MIBLMMP "1.lp" 0 0 <br/>
1000  1000  1.lp<br/>
51.617 51.6171    Time = 4.74645   y1= 7.88948   y2= 6.54251"<br/>

In addition, a text file named "Output.txt" will be generated that includes some more detail information of algorithm's solution:<br/>
1) Instance: name of the instance
2) nCons: number of the constraints
3) nVar: number of the variables 
4) GLB: global lower bound
5) GUB: global upper bound
6) NodeNumber: number of explored nodes
7) Time: run time
8) McDualImp: The number of times McCormick dual finder improves the solution
9) McCutImp: The number of times McCormick cut improves the solution 
10) McDualPieceNum: K^D 
11) McCutPieceNum: K^P
12) Iteration: number of iteration
13) McDualCallNum: number of times any McCormick operation is used
14) CplexDualNodeNum: Number of nodes CPLEX generates during the primal operation
15) CplexPrimalNodeNum: Number of nodes CPLEX generates during the dual operation
