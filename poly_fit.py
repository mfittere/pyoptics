from numpy import *

def solveconst(A,b,C,d):
  """Solve constrained least square problem
    minimize   || A x - b ||_2
    subject to Cx-d=0
    x: n
    A: m x n
    b: m
    C: l x n
    d: l
  using Lagrange multiplier
    L(x,l) = x.T A.T A x - 2 b.T A x + b.T b + l.T C x-l.T d
  where l is called the Lagrange multiplier. 
  The Kharush-Kuhn-Tucker equations
    nabla_x L(x,l) = 0, nabla_l L(x,l)=0
  yield:
  (A.T A   C.T )  (x)  = (A.T b)
  (C       0   )  (l)  = (d)
  see http://floatium.stanford.edu/ee263/notes/gradient-lagrange.pdf
  """
  nl,nx=C.shape
  m=hstack([dot(A.T,A),C.T])
  m=vstack([m,hstack([C,zeros((nl,nl))])])
  n=hstack([dot(A.T,b),d])
  sol=linalg.solve(m,n)
  return sol[:nx]

def makeA(x,N):
  """return [1,x0,x0**2,..,x0**N]
            [1,x1,x1**2,..,x1**N]
            ...
  """
  return column_stack([x**i for i in range(N+1)])

def makeAp(x,N):
  return column_stack([zeros(len(x))]+[i*x**(i-1) for i in range(1,N+1)])

def makeApp(x,N):
  z=[zeros(len(x))]*2
  z+=[i*(i-1)*x**(i-2) for i in range(2,N+1)]
  return column_stack(z)

def poly_val(p,x):
  return sum([p[i]*x**i for i in range(len(p))],axis=0)

def poly_print(p,x='x',power='**',mul='*'):
  res=['%+.10e%s%s%s%d'%(p[i],mul,x,power,i) for  i in range(len(p))]
  return ''.join(res)

def poly_fit(N,xdata,ydata,x0=[],y0=[],xp0=[],yp0=[],xpp0=[],ypp0=[]):
  """least square polynomial fit with constraints
  using Lagrangian multipliers.
  Parameters
  ----------
  xdata: n
  ydata: n
  N    : Degree of the fitting polynomial.
  x0,y0    : Fixed points.
  xp0,yp0  : Derivative at fixed points.
  xpp0,ypp0: Second derivative at fixed points.

  Returns
  -------
  p: Polynomial coefficients, highest power first.
  
  Notes
  -----
  Solves the least square problem:
    minimize   || A p - b ||_2
    subject to Cp-d=0
  with:
  p: polynomial
  A: matrix of [xdata**N ... xdata 1]
  b: ydata
  C: matrix of fixed points z=x0,xp0,xpp0
     with C=[C0,C1,C2
            [C2]
     with C0=[x0**N     ... x0 1]
          C1=[xp0**N-1  ... 1  0]
          C2=[xpp0**N-2 ... 0  0]
  d: fixed points with d=[y0,yp0,ypp0]
  """
  A=makeA(xdata,N)
  b=ydata
  C0=makeA(array(x0),N)
  C1=makeAp(array(xp0),N)
  C2=makeApp(array(xpp0),N)
  C=vstack([C0,C1,C2])
  d=hstack([y0,yp0,ypp0])
  p=solveconst(A,b,C,d)
  return p

if __name__=='__main__':
 x=linspace(0,1.0,101)
 p=array([2,3,-8,5])
 y=pol_val(p,x)
 n=.02*sin(2*pi*8*x)
 N=len(p)-1
 clf()
 plot(x,y+n)
 plot(x,poly_val(poly_fit(4,x,y+n,[0,1],[2,2]),x))
 A=makeA(x,N)
 lstsq(A,y+n)
 xx=array([0,1])
 yy=2+3*xx-8*xx**2+5*xx**3
 C=makeA(xx,N)
 d=yy
 p=solveconst(A,y+n,C,d)
 A=array([[2.,3],[1,2],[4,5]])
 b=array([2,3,1])
 C=array([[1,0],[0,1]])
 d=array([2,3])
 C=array([[1,0]])
 d=array([2])
 p=solveconst(A,b,C,d)






