# MT_curvature_simulation

Modeling the tip of an MT as an ABP and treating the motor binding events as a Poisson process.

---

## Main

### Input Parameters in the Code:
- **rate**: Rate of binding events
- **num_events**: Number of binding events
- **vel**: Velocity of the ABP (MT)
- **dt**: Time interval of the ABP
- **Drot**: Rotational diffusion of the ABP
- **Dtrans**: Translational diffusion of the ABP
- **theta**: Initial theta value
- **theta_0**: Maximum allowed absolute value of orientation angle at each ABP step. Models the bending stiffness of the MT.
- **guess**: Needed to numerically find roots of a transcendental equation.

### Output of the Code:
Two files, `Curvature.txt` and `Position_and_orientation.txt`, will be generated:

1. **Curvature.txt**
   - Contains 4 columns: (Event number, Curvature, Theta_final_of_ABP, Delta_theta_of_each_ABP).
   - The curvature can be plotted as a histogram using the `histogram.py` script.

2. **Position_and_orientation.txt**
   - Contains the position and orientation of the walker at all times.

### Workflow:
- The function `Time_step` is called first to generate the time intervals between each motor binding event. This also gives the number of steps of an ABP that will be called between each binding event.
- The function `ABP` simulates the motion of the tip of the MT between each motor binding event, considering the initial orientation of the ABP.
- After each ABP, the final position, displacement-to-distance ratio, and initial orientation angle uniquely determine the arc of motion of the ABP. This arc gives the curvature.
- The `To_get_theta_final` function determines the unique arc and provides the final angle of the ABP, which serves as the initial angle for the next ABP.
- The `Curvature` function computes the curvature of this unique arc.

---

## Functions

### 1. `Time_step(rate, num_events, dt)`
Generates the time interval between each binding event.

#### Input:
- **rate**: Rate of the motor binding event.
- **num_events**: Total number of motor binding events to simulate (integer).
- **dt**: Time step to calculate the number of steps an ABP needs to run in each interval.

#### Output:
- **inter_arrival_time**: Time interval between each binding event.
- **event_time**: Cumulative progress of time.
- **n_steps**: Integer value of the number of steps (calculated as `int(inter_arrival_time / dt)`).

### 2. `ABP(n_steps, theta, theta_0, dt, vel, Drot, Dtrans)`
Simulates the motion of the ABP.

#### Input:
- **n_steps**: Number of steps the ABP should take (output from `Time_step`).
- **theta**: Initial orientation angle of the walker.
- **theta_0**: Maximum allowed absolute value of orientation angle at each ABP step.
- **dt**: Time interval of each ABP step (same as in `Time_step`).
- **Drot**: Rotational diffusion.
- **Dtrans**: Translational diffusion.

### 3. `To_get_theta_final(theta_initial, guess, d_by_l)`
Determines the unique final orientation angle and curvature of the arc based on the ABP's position and displacement.

#### Input:
- **theta_initial**: Initial orientation of the ABP.
- **guess**: Initial guess value for solving the equation.
- **d_by_l**: Ratio of displacement to the distance traveled by the ABP.

#### Output:
- **theta_final**: Final orientation of the ABP, used as the initial angle for the next ABP.

#### Notes:
This function solves a transcendental equation using `fsolve`:
```
func(x, dis_len) = dis_len * x - np.sin(x)
```
Where `dis_len = d_by_l` and `x` is the theta variable.

### 4. `Curvature(theta, theta_final, length)`
Computes the curvature of the ABP between binding events.

#### Input:
- **theta**: Initial orientation angle.
- **theta_final**: Final orientation angle.
- **length**: Arc length of the ABP.

#### Output:
- **curvature**: Curvature value for each ABP between binding events.
