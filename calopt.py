from numpy import *

def beta(gamma):
  return sqrt(1-1/gamma**2)
#--- transfer matrices
def mat_drift(l,gamma):
  """transfer matrix drift"""
  return matrix([[1,l,0,0,0,0],[0,1,0,0,0,0],[0,0,1,l,0,0],[0,0,0,1,0,0],[0,0,0,1,l/(beta(gamma)*gamma)**2,0],[0,0,0,0,0,1]]) 
def mat_quad(l,k,gamma):
  """transfer matrix quadrupole"""
  g=sqrt(abs(k))
  cs=cos(g*l)
  ss=sin(g*l)
  ch=cosh(g*l)
  sh=sinh(g*l)
  if k>0: #focusing
    return matrix([[cs,ss/g,0,0,0,0],[-g*ss,cs,0,0,0,0],[0,0,ch,sh/g,0,0],[0,0,g*sh,ch,0,0],[0,0,0,1,l/(beta(gamma)*gamma)**2,0],[0,0,0,0,0,1]]) 
  if k<0: #focusing
    return matrix([[ch,sh/g,0,0,0,0],[g*sh,ch,0,0,0,0],[0,0,cs,ss/g,0,0],[0,0,-g*ss,cs,0,0],[0,0,0,1,l/(beta(gamma)*gamma)**2,0],[0,0,0,0,0,1]]) 

class tmat:
  """class for tranfer matrices"""
  def __init__(self,elem,l,k,m,gamma):
    self.gamma=gamma
    self.beta = beta(gamma)
    self.elem = elem
    self.l    = l
    self.k    = k
    self.mat=matrix(m)
  @classmethod
  def getmat(cls,l,k,elem,gammarel):
    if(elem=='drift'):
      return cls(elem,l,0,mat_drift(l,gammarel),gammarel)
    if(elem=='quad'):
      return cls(elem,l,k,mat_quad(l,k,gammarel),gammarel)
  def twiss(self,ax0,bx0,ay0,by0):
    m=self.mat
    t=matrix([]
    bx1=0;ax1=0;gx1=0;by1=0;ay1=0;gy1=0;
    for a0,b0,a1,b1,g1,i0,i1 in [[ax0,bx0,ax1,bx1,gx1,0,1],[ay0,by0,ay1,by1,gy1,2,3]]:
      g0=(1+a0**2)/b0
      b1=m[i0,i0]**2*b0-2*m[i0,i0]*m[i0,i1]*a0+m[i0,i1]**2*g0
      a1=-m[i0,i0]*m[i1,i0]*b0+(m[i0,i0]*m[i1,i1]+m[i0,i1]*m[i1,i0])*a0-m[i0,i1]*m[i1,i1]*g0
      g1=m[i1,i0]**2*b0-2*m[i1,i0]*m[i1,i1]*a0+m[i1,i1]**2*g0
      print a1,b1,g1
      if(g1!=(1+a1**2)/b1 ):#make a sanity check
        print 'WARNING: g1!=(1+a1**2)/b1'
    return ax1,bx1,gx1,ay1,by1,gy1
