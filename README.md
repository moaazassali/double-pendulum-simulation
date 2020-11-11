# double-pendulum-simulation
A double pendulum simulation in MATLAB which demonstrates the chaotic nature of the system


### Parameters to change the initial conditions of the double pendulum system:

#### Global variables: 
- L1: Length of first ball
- L2: Length of second ball
- m1: Mass of first ball
- m2: Mass of second ball
- g: Acceleration due to gravity
  

#### Initial conditions of the system:
- th1(1): Initial angle of the first ball from the vertical in radians
- th2(2): Initial angle of the second ball from the vertical in radians
- w1(1): Initial angular velocity of the first ball
- w2(1): Initial angular velocity of the second ball
  
  
#### Properties of the simulation runtime:
- time: How long the simulation will run in seconds (represents the system's simulation accurately for "time" seconds; doesn't actually run in real-time "time" seconds)
- timeStep: Step size of time increment for the ODE solvers
  

#### Choice of ODE solver:
- Euler's Method: A first-order ode solver method which is less accurate and not desirable for chaotic systems given the large impact small changes lead to
- 4th-order Runge-Kutta (RK4): A fourth-order ode solver method that gives much more accurate results with minimal truncation errors

The preferred ode-solver method can be chosen by commenting out the other method in the main file. By default, "Euler's Method" is commented out and RK4 is used.
  
