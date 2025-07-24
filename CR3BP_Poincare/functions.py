import numpy as np
# Function to convert from the inertial system to the rotating sytem
def fix2rot(xf, t):
  ct = np.cos(t)
  st = np.sin(t)

  xr = np.zeros[6]

  xr[0] = xf[0]*ct + xf[1]*st
  xr[1] = xf[1]*ct - xf[0]*st
  xr[2] = xf[2]

  xr[3] = xf[3]*ct + xf[1]*ct - xf[0]*st + xf[4]*st
  xr[4] =-xf[0]*ct + xf[4]*ct - xf[3]*st - xf[1]*st
  xr[5] = xf[5]

  return xr

# Function to convert from the rotating sytem to the inertial system

def rot2fix(xr, t):
  ct = np.cos(t)
  st = np.sin(t)

  xf = np.zeros[6]

  xf[0] = xr[0]*ct - xr[1]*st
  xf[1] = xr[1]*ct + xr[0]*st
  xf[2] = xr[2]
  xf[3] =-xr[0]*st - xr[4]*st + xr[3]*ct - xr[1]*ct
  xf[4] = xr[3]*st - xr[1]*st + xr[0]*ct + xr[4]*ct
  xf[5] = xr[5]

  return xf
