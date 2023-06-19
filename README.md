# *Using the  Conjugate gradient method to solve the Poisson's Equation in OpenMP*
---
## 1. Introduction
Conjugate gradient method is an algorithum method to solve a large system of _linear equations_.
In this project, we use this method to solve a simple Poisson's equation problem. 
To increase the efficiency of the process, we use the OpenMP to calculate some difficult steps in parallel.
## 2. Technologies
* C++
* OPENMP
Project is created with:
* T110E5, jsp (Ning Chen): Create the first version and second version of CG method.
* Kuo-Jui: Do the final edition and add the parallel parts.
## 3. Setup
Download the document and run the program below:

`$ make`

`./a.out`
## 4. Result
We choose N as 64, 80, 96, 112, and 128 to calculate how well the parallelization does. 
The number of threads was chosen from 1 to 8. All of the settings were changing on hands (They are not in the script.)
It is lucky for us to have a nice CG method performance. Nonetheless, from the plot of efficiency to the number of threads, we know the former result was accidentally correct. The issue is that the iteration will change after doing the same process for several times. It may mean that our code or parallelization wasn't good enough. There may be some mistakes when making the **_A_** matrix. To make our code better, we should try to sparse the A matrix first because there are lots of zeros inside it.
## 5. Reference
* https://github.com/blackcata/Poisson_Equation/tree/master
* https://sites.cs.ucsb.edu/~gilbert/cs240a/old/cs240aSpr2011/hw2
* An Introduction to the Conjugate Gradient Method Without the Agonizing Pain, Jonathan Richard Shewchuk
August 4, 1994
____
## Ps.
Our SOR method didn't put at here. We use the old version on Python, which is the same as the homework before. However, it was still a disaster.

**_Readme is written by Kuo-Jui. If there is any question, please mail me._**
