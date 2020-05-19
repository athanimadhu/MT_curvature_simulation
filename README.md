# MT_curvature_simulation
Modeling the tip of an MT as an ABP and treating the motor binding events as a poisson process.
——————————————————————————————————————————————————————————————————————————————————————————————————————————————---------------
Main
——————————————————————————————————————————————————————————————————————————————————————————————————————————————---------------
INPUT PARAMETERS IN THE CODE:
rate - Rate of binding events
num_events - Number of binding events
vel - Velocity of the ABP (MT)
dt - Time interval of the ABP
Drot - Rotational diffusion of the ABP
Dtrans - Translational diffusion of the ABP
theta - Initial theta value
theta_0 - Maximum allowed absolute value of orientation angle at each ABP step. Models the bending stiffness of the MT.
guess - Needed to numerically find roots of a transcendental equation.

OUTPUT OF THE CODE: 2 files “Curvature.txt” and “Position_and_orientation.txt” will be gererated.
1. Curvature.txt - Has 4 columns (Event number, Curvature, Theta_final_of_ABP, delta_theta_of_each_ABP). We need to find curvature. The histogram of the 2nd column can be plotted using ‘histogram.py’ that is attached.
2. Position_and_orientation.txt - Contains the position and orientation of the walker at all times. 

The function 'Time_step' is called first to generate the time intervals between each motor binding event. This also gives the number of steps of an ABP that will be called between each binding event. The function 'ABP' simulates the motion of the tip of the MT between each motor binding event. This takes into account the initial orientating on the ABP. After each ABP, the final position, the displacement to distance ratio and the initial orientation angle will uniquely determine the arc of motion of the ABP. This arc gives a curvature. It will also intern set the initial orientation angle of the next ABP which takes place in the next time interval between the binding events. The functions 'To_get_theta_final' determines the unique arc, and intern also gives the final angle of the ABP which serves as the initial angle for the next ABP. The function 'Curvature' computes the curvature of this unique arc. 
——————————————————————————————————————————————————————————————————————————————————————————————————————————————---------------
Functions
——————————————————————————————————————————————————————————————————————————————————————————————————————————————---------------
There are 4 parts to this code defined as functions, namely 'Time_step', 'ABP', 'To_get_theta_final' and 'Curvature'.

1. Time_step(rate, num_events(int), dt): This gives the time interval between each binding event.
Input: 
rate - Rate of the event occurring (Motor binding event).
num_events - Total number of motor binding events to simulate (this is an integer number)
dt - To calculate the number of time steps an ABP needs to be run in each time interval of time.
Output:
inter_arrival_time - Time interval between each binding event.
event_time - This is the cumulative progress of time.
n_steps = int(inter_arrival_time/dt). Integer value of the number of time steps. (this will be an input function - ABP)

2. ABP(nsteps(int), theta, theta_0, dt, vel, Drot, Dtrans):
Input:
nsteps - Number of steps the ABP should take. This is an output from the previous function.
theta - Initial orientation angle of the walker
theta_0 - Maximum allowed absolute value of orientation angle at each ABP step. Models the bending stiffness of the MT.
dt - Time interval of each ABP step. This parameter should be same as the dt in the previous function.
Drot - Rotational diffusion.
Dtrans - Translational diffusion.

3. To_get_theta_final(theta_initial, guess, d_by_l): Under the assumption that the MT is like an arc of a circle between each binding events, the final position of the ABP and the distance travelled by the ABP can give a unique final orientation angle and the curvature of the arc. This arc can be analytically found to be a single variable linear transcendental equation which needs to be solved to find the final theta at each binging event. This function uses \fsolve (an inbuilt function to find roots).
Input:
theta_initial - Initial orientation of the ABP.
guess - An initial guess value needs to be provided to solve the equation.
d_by_l - This is the ratio of the displacement to the distance travelled by the ABP.
Output:
theta_final - The final orientation of the ABP. This will serve as the theta_initial for the next ABP in the next poisson distributed time interval.
Note: The transcendental equation that this function solves is given by another function
‘func(x,dis_len)’ --> dis_len*x - np.sin(x). Where dis_len = d_by_l and x is the theta variable.

4. Curvature(theta, theta_final, length): This gives the curvature of the ABP between each binding event. There is one curvature value associates with every ABP (between each binding event).
Input:
theta - Initial orientation of the ABP.
theta_final - Find orientation of the ABP determined from solving the transcendental equation using the previous function. 
length - distance travelled by the ABP.
Output:
curvature - The curvature of the arc that approximates the movement of the MT between two binding events.
