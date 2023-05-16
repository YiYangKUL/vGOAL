# vGOAL
Get the results of the vGOAL approach:  
Put AMR_Spec.py, Interpreter.py, MG_01.txt, MG_02.txt, MG_03.txt, MG_04.txt, MG_05.txt, and MG_012.txt into one folder.  
Compile and run AMR_Spec.py, Interpreter.py in Pycharm or simply in the terminal.  
To test an error-free run for the given task, give inputs: 1. 
To test a non-fatal error run, we use ‘E1’ as an example, give inputs: 2. 
To test a non-fatal error run, we use ‘E4’ as an example, give inputs: 3.  
Example:  
1: an error-free run test  
2: a non-fatal error test  
3: a fatal error run test  
Which test do you want to pereform? (Please enter 1, 2, or 3.)3.

Output explanation:  
Cycle1 
A1 move3(6,3)  
A2 move3(8,4)\
The decision-making generation takes: 0.8210129737854004  
Explanation: For the first reasoning cycle, the interpreter generates two safe decisions for A1 and A2, and the decision-making generation took 0.8210129737854004s.

Cycle2  
Waiting for new information!  
Explanation: For the first reasoning cycle, the interpreter cannot generate any new decisions, because A1 and A2 are executing the durative actions, and A3 is waiting for the location permission.  


Cycle30   
A3 changes the goal due to the error information: E1.  
A3 move5(3,5)  
Explanation:  In the fatal-error run, A3 encounters a non-fatal error, E1, i.e., docking error, A3 abandons the current delivery goal, and it goes to the waiting point.   


Cycle28. 
A3 crashed!  
A1  received new goals.  
Waiting for new information!  
Explanation:  In the fatal-error run, A3 encounters a fatal error, and its remaining goals are redistributed to A1.   


Cycle41  
No active goals!  
Explanation: There is no active goals in the autonomous system.  
