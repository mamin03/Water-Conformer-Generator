import sys
import numpy as np
from pyquaternion import Quaternion
from random import random
from math import acos, sin, radians, cos

def get_Point(P1,P2,Radius):
    """
    This function returns a point on a line that
    connects the metal and the oxygen for placing
    the hydrogens
    """

    V = P1 - P2
    V = V / np.linalg.norm(V)
    e = np.array([0,0,0])
    e[np.argmin(np.abs(V))] = 1
    V1 = np.cross(e, V)
    V1 = V1 / np.linalg.norm(V)
    V2 = np.cross(V, V1)
    s=0
    P = P1 + Radius*( np.cos(s)*V1 + np.sin(s)*V2 )
    return P

def calc_t(a1,a2):
   """
   This function returns the value of the required
   extention to the line that connects the metal and
   the oxygen to generate the correct hydrogen position
   """

   a = 0
   b = 0
   c = 0
   for i,j in zip(a1,a2):
     a += (j-i)**2
     b += 2*(j-i)*(i-j)
     c += (j-i)**2
   l = 0.96*sin(radians(54))
   return np.roots([a,b,c-l**2])

def make_H(a1, a2, theta=0., n=100, rand = False):
   """
   This function generates the hydrogen positions
   """

   t = max(calc_t(a1,a2))
   v = np.array([a1[0] + (a2[0]-a1[0])*t,
                 a1[1] + (a2[1]-a1[1])*t,
                 a1[2] + (a2[2]-a1[2])*t])
   s2 = open("step2_out.pdb","w")
   v1 = v-a2
   r = 0.96*cos(radians(54))
   h1 = get_Point(v,a2,r)
   h1 = h1-a2
   print(2*n+3)
   print()
   print("Mg", a1[0], a1[1], a1[2])
   print("O", a2[0], a2[1], a2[2])
   for i in range(n):
      if rand:
          q = Quaternion.random()
      else:
          q = Quaternion(axis=v1/np.linalg.norm(v1), angle=2*np.pi*random())
      qinv = Quaternion(axis=v1/np.linalg.norm(v1), angle=np.pi)
      h2 = qinv.rotate(h1) 
      x1 = h1[0]+a2[0];y1 = h1[1]+a2[1]; z1 = h1[2]+a2[2]
      x2 = h2[0]+a2[0];y2 = h2[1]+a2[1]; z2 = h2[2]+a2[2]

      AtomO = '{:07.3f}'.format(a2[0])+" "+'{:07.3f}'.format(a2[1])+" "+'{:07.3f}'.format(a2[2])
      AtomH1 = '{:07.3f}'.format(x1)+" "+'{:07.3f}'.format(y1)+" "+'{:07.3f}'.format(z1)
      AtomH2 = '{:07.3f}'.format(x2)+" "+'{:07.3f}'.format(y2)+" "+'{:07.3f}'.format(z2)
     
      print("H1 "+AtomH1)
      print("H2 "+AtomH2)

      s2.write("ATOM      1  O   HOH X1111_"+'{:03}'.format(i)+" "+AtomO+ "   1.520      -0.800      01O000M000"+"\n")
      s2.write("ATOM      2 1H   HOH X1111_"+'{:03}'.format(i)+" "+AtomH1+"   1.100       0.400      01O000M000"+"\n")
      s2.write("ATOM      3 2H   HOH X1111_"+'{:03}'.format(i)+" "+AtomH2+"   1.100       0.400      01O000M000"+"\n")

      h1 = q.rotate(h1)

   return np.round(h1,3), np.round(h2,3)


if __name__ == "__main__":
   f = open(sys.argv[1])
   t = f.readline().split()
   a1 =  np.array([float(t[6]), float(t[7]), float(t[8])])
   t = f.readline().split()
   a2 =  np.array([float(t[6]), float(t[7]), float(t[8])])
   h1, h2 = make_H(a1, a2, n=100)
