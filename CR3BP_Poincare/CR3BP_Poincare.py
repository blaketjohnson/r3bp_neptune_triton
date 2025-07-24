'''
! Purpose:

! This program is used to generate Poincaré Surface of Section (or Poincaré maps)
! of systems in the Circular Restricted Three-Body Problem (CR3BP).

! Author: Diogo Merguizo Sanchez, 2023
! Based on the 2015 version of the program in FORTRAN.
!
! Based on orginal routines (Fortran) written by Roger Broucke
! Some references:
! * Broucke, R. A., Numerical Integration of Periodic Orbits in
!   the Main Problem of Artificial Satellite Theory, Celestial
!   Mechanics & Dynamical Astronomy, Volume 58, Issue 2, pp.99-123
!   DOI: 10.1007/BF00695787
! * Broucke, R. A., Periodic Orbits in the Restricted Three-Body Problem
!   With Earth-Moon Masses, NASA Technical Report 32-1168, 1968.

'''
import numpy as np
from scipy.integrate import solve_ivp
import timeit
from datetime import datetime
import os
# import csv
from functions import fix2rot, rot2fix
from control import pause
import matplotlib.pyplot as plt
 
print("\nInit", datetime.now())

#### Reading the input parameters from file ####
file = open('/Users/blakejohnson/Documents/rtbp_neptune_triton/source/CR3BP_Poincare/input_py.in', 'r')

folder = file.readline().split()[0]

C_range = file.readline().split()[0:3]
C0 = float(C_range[0])
CF = float(C_range[1])
dC = float(C_range[2])

mu = float(file.readline().split()[0])
mu_star = 1.0 - mu
r = float(file.readline().split()[0])

r_max = float(file.readline().split()[0])
r1_min = float(file.readline().split()[0])
r2_min = float(file.readline().split()[0])

t_range = file.readline().split()[0:2]
tlim = float(t_range[0])
dt = float(t_range[1])

x_range = file.readline().split()[0:3]
xi = float(x_range[0])
xf = float(x_range[1])
dx = float(x_range[2])

save_traj = int(file.readline().split()[0])
step_traj = int(file.readline().split()[0])

plot_traj = int(file.readline().split()[0])

print ('--------------- Initial parameters ---------------')
print(f'Folder to keep the results: {folder}')
print(f'C0 = {C0}, CF = {CF}, DC = {dC}')
print(f'mu = {mu}')
print(f'Distance between the primaries: {r} km')
print(f'Distance to consider as escape: {r_max} km')
print(f'Minimum distance from M1: {r1_min} km')
print(f'Minimum distance from M2: {r2_min} km')
print(f'tlim = {tlim}, dt = {dt}')
print(f'xi = {xi}, xf = {xf}, dx = {dx}')
if save_traj != 0:
  print(f'Saving trajectory every {step_traj} points')
print('\nDouble check all these values.')
print('-----------------------------------------------')

print("\nRunning...")
start = timeit.default_timer()

# Loop to vary the Jacobi constant
ni = round((CF - C0)/dC) # number of points

for ij in range(0, ni + 1):
  
  CJ = C0 + float(ij)*dC
  # print (f"Running C = {CJ}\n")

  # Check the existence of the folder where the results will be stored.
  # Create the folder in the case it doesn't exist.
  isExist = os.path.exists(folder)
  if not isExist:
    os.mkdir(folder)

  # number of point in each trajectory
  nx  = round((xf - xi)/dx)

  for ii in range(0, nx + 1):
    x0 = xi + float(ii)*dx

    x = np.zeros(6)
    x[0] = x0
    x[1] = 0.0 # this is the hyperplane
    x[2] = 0.0 # considering planar orbits

    # Radius from the first primary (M1)
    R1 = np.sqrt((x[0] + mu)**2 + x[1]**2 + x[2]**2)

    # Radius from the second primaty (M2)
    R2 = np.sqrt((x[0] - mu_star)**2 + x[1]**2 + x[2]**2)

    R = np.sqrt(x[0]**2 + x[1]**2 + x[2]**2)

    # Critical conditions (events)
    cond  = 0
    if R1 <= r1_min/r:
        # print(f"Collision with M1, R1_0 = {R1*r} km")
        cond = 1
    if R2 <= r2_min/r:
        # print(f"Collision with M2, R2_0 = {R2*r} km")
        cond = 2
    if R >= r_max/r:
        # print(f"Escape from the system, R_0 = {R*r} km")
        cond = 9

    # Calculation of the velocity. Since the hyperplane is in y = 0, 
    # the massless 3rd body has only vy, calculated using x and CJ
    arg = - CJ + x[0]**2 + x[1]**2 + 2.0*(mu_star/R1 + mu/R2)
    if arg < 0:
      continue
    vy0 = np.sqrt(arg)
    x[3] = 0.0
    x[4] = vy0
    x[5] = 0.0

    # print(f"CJ = {CJ}, x_0 = {x0}, vy_0 = {vy0}")

    # set array with the initial conditions (to be used by solve_ivp)

    f0 = np.zeros(6)
    f0[0] = x[0]
    f0[1] = x[1]
    f0[2] = x[2]
    f0[3] = x[3]
    f0[4] = x[4]
    f0[5] = x[5]

    # critical = 1

    def critical_events(t, f):

      global cond, critical

      # Distance from M1
      R1 = np.sqrt((mu + f[0])**2 + f[1]**2 + f[2]**2)

      # Distance from M2
      R2 = np.sqrt((f[0] - mu_star)**2 + f[1]**2 + f[2]**2)

      # Distance from the center of mass
      R = np.sqrt(f[0]**2 + f[1]**2 + f[2]**2)

      # cond = 0
      if R1 <= r1_min/r:
        # print(f"Collision with M1, R1 = {R1} n.d.")
        cond = 1
      if R2 <= r2_min/r:
        # print(f"Collision with M2, R2 = {R2} n.d.")
        cond = 2
      if R >= r_max/r:
        # print(f"Escape from the system, R = {R} n.d.")
        cond = 9
      
      if cond != 0:
        critical = 0
      else:
        critical = 1
      
      return critical

    critical_events.direction = 1
    critical_events.terminal = True

    def dFdt (t, f):

      # Distance from M1
      R1 = np.sqrt((mu + f[0])**2 + f[1]**2 + f[2]**2)

      # Distance from M2
      R2 = np.sqrt((f[0] - mu_star)**2 + f[1]**2 + f[2]**2)

      dx = f[3]
      dy = f[4]
      dz = f[5]

      ddx = 2.0*f[4] + f[0] - (mu_star*(f[0] + mu))/R1**3 -(mu*(f[0] - mu_star))/R2**3
      ddy =-2.0*f[3] + f[1] - (mu_star*f[1])/R1**3 - (mu*f[1])/R2**3
      ddz =-(mu_star*f[2])/R1**3 - (mu*f[2])/R2**3

      return [dx, dy, dz, ddx, ddy, ddz]
    
    nt = int(tlim/dt)
    tspan = np.arange (0, tlim + dt, dt)
#DOP853
    solution = solve_ivp (dFdt, [0.0, tlim + dt], f0, events=[critical_events], method='DOP853',\
                          t_eval=tspan, first_step =dt/100.0, rtol = 1e-12, atol = 1e-12)
    # solution = solve_ivp (dFdt, [0.0, tlim + dt], f0, method='LSODA',\
    #                       t_eval=tspan, rtol = 1e-12, atol = 1e-12)
    
    state = solution.y

    # print(f"x_0 = {round(x[0], 5):.5f}, cond = {cond}, critical = {critical}, sol_status = {solution.status}")
    print(f"x_0 = {round(x[0], 5):.5f}, cond = {cond}")

    # pause()

    # print(state.shape)
    # print(solution.status)

################################################################################
# Data analysis
    def Poincare (state):
        
        # global CJ, x0, vy0
        
        xa = np.zeros(6)
        xm = np.zeros(6)

        xa[0] = state[0, 0]
        xa[1] = state[1, 0]
        xa[2] = state[2, 0]
        xa[3] = state[3, 0]
        xa[4] = state[4, 0]
        xa[5] = state[5, 0]

        for it in range(1, nt + 1):

          if state[1, it]*xa[1] < 0 and state[1, it] > 0:
            xm[0] = (state[0, it] + xa[0])/2.0
            xm[1] = (state[1, it] + xa[1])/2.0
            xm[2] = (state[2, it] + xa[2])/2.0
            xm[3] = (state[3, it] + xa[3])/2.0
            xm[4] = (state[4, it] + xa[4])/2.0
            xm[5] = (state[5, it] + xa[5])/2.0

            with open(file1, "a") as file_PS:
                np.savetxt(file_PS, xm, newline=' ')
                file_PS.write(str(round(x0, 5)) + ' ' + str(round(vy0, 14)) + \
                              ' ' + str(round(CJ, 5)) + "\n")
            file_PS.close()

          xa[0] = state[0, it]
          xa[1] = state[1, it]
          xa[2] = state[2, it]
          xa[3] = state[3, it]
          xa[4] = state[4, it]
          xa[5] = state[5, it]

    if solution.status == 1:
      # print("There was a critical event -- excluding orbit from PS.\n")
      # print("Check PS_log.dat file for details.\n")
      continue
    else:
      # Calculate and store the Poincaré Section
      # Create a string with the Jacobi constant to use it in a file name
      aux1 = str(round(CJ, 5))
      aux2 = str(round(x0, 5))
      file1 = folder + '/' + 'PY-C' + aux1 + 'Xi' + aux2 + '.dat'

      # print (it, state[9, it], xa[1], state[9, it]*xa[1])
      if x0 == xi:
        # Test if the file alreay exists.
        isExist = os.path.exists(file1)
        if isExist:
          print (file1, "already exists. Overwriting.")
          os.remove(file1)
      
      # Generate the Poincare Section
      Poincare(state)

#-------------------------------------------------------------------------------
      # Diminishing the number of points. Save every step_traj points.
      aux_state = np.transpose(state[:, ::step_traj])
      # print(aux_state.shape)
      
      # Plot trajectory (rotating frame) if desired
      if plot_traj == 1:
        plt.scatter(aux_state[:,0],aux_state[:,1], s=0.1)
        plt.show()
        # pause()


      # Save trajectories if desired
      if save_traj == 1:
          file2 = folder + '/' + 'TR-C' + aux1 + 'Xi' + aux2 + '.dat'
          isExist = os.path.exists(file2)
          if isExist:
            print (file2, "already exists. Overwriting.")
            os.remove(file2)
          with open(file2, "a") as file_traj:
                np.savetxt(file_traj, aux_state)
          file_traj.close()
          # pause()

#-------------------------------------------------------------------------------
# End of execution - print time of execition and date       
print("End:", datetime.now())
stop = timeit.default_timer()
runtime = stop - start
if runtime < 60:
    print(f"Runtime = {runtime:.2f} seconds.\n")
elif runtime >= 60 and  runtime < 3600:
    print(f"Runtime = {runtime/60:.2f} minutes.\n")
else:
    print(f"Runtime = {runtime/3600:.2f} hours.\n")
################################################################################