#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 21:29:43 2020

@author: athani_m

This is a code to simulate the walk of the tip of an MT. binding to the motor is treated as an Possion distributed event.
 
"""

import numpy as np
from scipy.optimize import fsolve

#Choosing the time interval between binding events and the number of steps an ABP should run between these binding events
def Time_step(rate, num_events, dt):
    "Gives the number of time steps for each ABP"
    inter_arrival_time = np.empty(num_events)
    event_time = np.empty(num_events+1)
    event_time[0] = 0
    step = 0
    n_steps = np.empty(0)
    p = np.random.random(num_events)
    for i in range(num_events):
        inter_arrival_time[i] = -np.log(1.0 - p[i])/rate
        event_time[i+1] = event_time[i] + inter_arrival_time[i]
        step = int(inter_arrival_time[i]/dt)
        n_steps = np.append(n_steps,step)
    return inter_arrival_time, event_time, n_steps

def ABP(nsteps, theta, theta_0, dt, vel, Drot, Dtrans):
    "This returns all the positions, orientaitons and the total distance traveled by an ABP"
    number_steps = 0;
# initialize arrays that store x,y and theta values, as well as initial particle position and angle
    xvec=np.zeros(0); yvec=np.zeros(0); thetavec = np.zeros(0) ;
    x=0.0; y = 0.0; delta_theta = np.zeros(0); x_old = 0; y_old = 0; length =0;
    xvec = np.append(xvec,x); yvec = np.append(yvec,y); thetavec = np.append(thetavec, theta)
    length = 0
    while number_steps < nsteps :
        # diffusive/random steps
        s = np.random.normal(0,0.05, 3)
        dx = np.sqrt(2*Dtrans*dt)*2*s[0]; 
        dy= np.sqrt(2*Dtrans*dt)*2*s[1]; 
        dtheta = np.sqrt(2*Drot*dt)*(2*np.pi)*s[2];
        if abs(dtheta) < theta_0:
            # update coordinates (including ballistic step)
            delta_theta = np.append(delta_theta, dtheta)
            x_old = x
            y_old = y
            x += vel*dt*np.cos(theta) + dx 
            y += vel*dt*np.sin(theta) + dy
            xvec = np.append(xvec,x); yvec = np.append(yvec,y)  # store successive positions in arrays
            theta += dtheta
            thetavec = np.append(thetavec, theta)
            number_steps += 1
            length += np.sqrt((x-x_old)**2 + (y-y_old)**2)
    return length, xvec, yvec, thetavec, delta_theta

def func(x,dis_len):
    return dis_len*x - np.sin(x)

def To_get_theta_final(theta_initial, guess, d_by_l):
    "Returns the inital theta value for the next ABP"
    x0 = fsolve(func, guess, args=(d_by_l),full_output=0)
    root = x0.item()
    if (root == 0 or d_by_l >= 1.0):
        theta_final = theta_initial
    elif (abs(2*root + theta_initial) < abs(-2*root + theta_initial)):
        theta_final = 2*root + theta_initial
    else:
        theta_final = -2*root + theta_initial
    return theta_final

def Curvature(theta, theta_final, length):
    "Using the inital and final theta of an ABP, gives the curvature"
    curvature = (abs(theta_final-theta))/length
    return curvature

#Main
rate = 2
num_events = 1000
vel = 2.0; dt = 0.005; Drot = 0.1; Dtrans = 0.1;
theta = 0.0
theta_0 = 0.001
guess = 3.0

file1 = open("Position_and_orientation.txt", "w+")
file2 = open("Curvature.txt","w+")
k = 0
inter_arrival_time, event_time, n_steps = Time_step(rate, num_events, dt)

for i in range(len(n_steps)):
    steps = int(n_steps[i])
    length, xvec, yvec, thetavec, delta_theta = ABP(steps, theta, theta_0, dt, vel, Drot, Dtrans)
    
    for j in range(len(xvec)):
        file1.write('%f   %f   %f \n' %(xvec[j],yvec[j],thetavec[j]))
    
    distance = np.sqrt((xvec[0]-xvec[-1])**2 + (yvec[0]-yvec[-1])**2)
    if length != 0:
        d_by_l = distance/length 
        theta_final_ABP = To_get_theta_final(theta, guess, d_by_l)
        curvature_ABP = Curvature(theta, theta_final_ABP, length)
        delta_theta_ABP = abs(theta_final_ABP - theta)
        file2.write('%d %f  %f  %f\n' %(k ,curvature_ABP, theta_final_ABP, delta_theta_ABP))
        k += 1
        theta = theta_final_ABP

print("DONE")
file1.close()
file2.close()