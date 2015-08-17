from numpy import *

def betarel(gammarel):
  return sqrt(1-1/gammarel**2)

#--- transfer matrices
def mat_drift(l,gammarel=7460.5228):
  """transfer matrix drift"""
  return matrix([[1,l,0,0,0,0],[0,1,0,0,0,0],[0,0,1,l,0,0],[0,0,0,1,0,0],[0,0,0,0,1,l/(betarel(gammarel)*gammarel)**2],[0,0,0,0,0,1]]) 
def mat_quad(l,k,gammarel=7460.5228):
  """transfer matrix quadrupole"""
  g=sqrt(abs(k))
  cs=cos(g*l)
  ss=sin(g*l)
  ch=cosh(g*l)
  sh=sinh(g*l)
  if k>0: #focusing
    return matrix([[cs,ss/g,0,0,0,0],[-g*ss,cs,0,0,0,0],[0,0,ch,sh/g,0,0],[0,0,g*sh,ch,0,0],[0,0,0,1,l/(betarel(gammarel)*gammarel)**2,0],[0,0,0,0,0,1]]) 
  if k<0: #focusing
    return matrix([[ch,sh/g,0,0,0,0],[g*sh,ch,0,0,0,0],[0,0,cs,ss/g,0,0],[0,0,-g*ss,cs,0,0],[0,0,0,1,l/(betarel(gammarel)*gammarel)**2,0],[0,0,0,0,0,1]]) 

def gamma(a,b):
  """calculate the twiss parameter gamma from
  alpha (a) and beta (b)"""
  return (1+a**2)/b
  

def test_tw(ax,bx,cx,ay,by,cy):
  if cx!=gamma(ax,bx):
    print 'WARNING: cx!=(1+ax**2)/bx'
  if cy!=gamma(ay,by):
    print 'WARNING: cy!=(1+ay**2)/by'
def m_to_tw(m):
  """get twiss parameters from matrix"""
  cx=m[0,0];ax=m[1,0];bx=m[1,1]
  cy=m[2,2];ay=m[3,2];by=m[3,3]
  test_tw(ax,bx,cx,ay,by,cy)
  return ax,bx,ay,by
class mtwiss(matrix):
  """class for transfering twiss parameters
   class object is the matrix:
   t=[[gx,ax],[ax,bx]]"""
  def __new__(self,ax,bx,ay,by):
    """create the twiss matrix 
    t=[[cx,ax],[ax,bx]]
    """
    cx=gamma(ax,bx);cy=gamma(ay,by);
    self.ax=ax;self.bx=bx;self.cx=cx;
    self.ay=ay;self.by=by;self.cy=cy;
    return super(mtwiss,self).__new__(self,[[cx,ax,0,0,0,0],[ax,bx,0,0,0,0],[0,0,cy,ay,0,0],[0,0,ay,by,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]]) 
 # @classmethod
 # def from_matrix(cls,m):
 #   cx=m[0,0];ax=m[1,0];bx=m[1,1]
 #   cy=m[2,2];ay=m[3,2];by=m[3,3]
 #   test_tw(ax,bx,cx,ay,by,cy)
 #   print 'a'  
 #   return cls(ax,bx,ay,by)
  def trans(self,r):
    """transport twiss parameters by transfer matrix r"""
    a1=((r.getT()).getI())*self*(r.getI())
    ax,bx,ay,by=m_to_tw(a1)
    return mtwiss(ax,bx,ay,by)
