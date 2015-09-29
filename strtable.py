"""module to plot the squeeze tables"""

import re

import matplotlib.pyplot as pl
from numpy import *
import os as os

from pydataobj import dataobj
import tfsdata

from poly_fit import poly_fit, poly_print, poly_val

def iter(a):
  """make an iterable out of a"""
  return array(a).flatten()
def bdump(bet,alf,beam='1'):
  al_dump=761
  if beam=='1':
    bdump=bet-2*al_dump*alf+al_dump**2*(1+alf**2)/bet
  if beam=='2':
    bdump=bet+2*al_dump*alf+al_dump**2*(1+alf**2)/bet
  return bdump
def def_subplot(nsub):
  """define a grid of subplots for
  *nsub* subplots"""
  ps=3#dimensions of each subplot
  if nsub==1:
    pl.figure()
    lsub=[[1,1,1]]
  if nsub==2:
    pl.figure(figsize=(2*ps,ps))
    lsub=[[2,1,1],[2,1,2]]
  if 2<nsub<=4:
    pl.figure(figsize=(2*ps,2*ps))
    lsub=[[2,2,n] for n in range(1,nsub+1)]
  if 4<nsub<=6:
    pl.figure(figsize=(3*ps,2*ps))
    lsub=[[3,2,n] for n in range(1,nsub+1)]
  if 6<nsub<=9:
    pl.figure(figsize=(3*ps,3*ps))
    lsub=[[3,3,n] for n in range(1,nsub+1)]
  if 9<nsub<=12:
    pl.figure(figsize=(4*ps,3*ps))
    lsub=[[4,3,n] for n in range(1,nsub+1)]
  return lsub
def mk_subplot(l):
  a,b,c=l
  return pl.subplot(a,b,c)

class StrTable(dataobj):
  scale=23348.89927
  @classmethod
  def open(cls,fn):
    obj=cls(tfsdata.open(fn))
    return obj
  def get_vars(self,reg):
    rxp=re.compile(reg)
    return sorted(l for l in self.keys() if rxp.search(l))
  def get_kq(self,n):
    return self.get_vars(r'kq[xt]?l?%da?\.'%n)
  def get_phases(self):
    out=self.get_vars(r'mu[xy]ip[1-8]b[12]$')
    out+=self.get_vars(r'mu[xy]ip[1-8]b[12]_l')
    return out
  def get_acb(self,n,knob='on_sep'):
    out=[]
    for n in self.get_vars('acb.*%d\.[lr][1-8].*%s'%(n,knob)):
      if sum(abs(self[n]))>0:
        out.append(n)
    return sorted(out)
  def plot_acb(self,n,knob,n1,n2,x=None,scale=1,brho=None):
    if brho is None:
        scale*=1e6
    else:
        scale*=brho
    if x is None:
      xv=arange(len(self[self.keys()[0]][n1:n2]))
    else:
      xv=self[x][n1:n2]
    ks=self.get_acb(n,knob)
    for k in ks:
      pl.plot(xv,self[k][n1:n2]*scale,label=k)
    pl.legend(loc=0,frameon=False)
    if x is not None:
      pl.xlabel(x)
    if  brho is None:
      pl.ylabel(r'angle [$\mu$rad]')
    else:
      pl.ylabel('k0l [Tm]')
  def get_triplet(self):
    return self.get_vars('kqx\.r[2815]|kt?qx[123].*r[15]')
  def plot_triplet(self,n1,n2,x=None):
    scale=StrTable.scale
    if x is None:
      xv=arange(len(self[self.keys()[0]][n1:n2]))
    else:
      xv=self[x][n1:n2]
    ks=self.get_triplet()
    for k in ks:
      pl.plot(xv,abs(self[k][n1:n2]*scale),label=k)
    pl.legend(loc=0,frameon=False)
    if x is not None:
      pl.xlabel(x)
    pl.ylabel('k [T/m]')
  def plot_2in1(self,kq,n1,n2,x=None,sign=False,ylab='k [T/m]'):
    scale=StrTable.scale
    if x is None:
      xv=arange(len(self[self.keys()[0]][n1:n2]))
    else:
      xv=self[x][n1:n2]
    for k in self.get_kq(kq):
      kv=self[k][n1:n2]*scale
      if sign or kv[0]>0:
        pl.plot(xv,kv,label=k)
      else:
        pl.plot(xv,-kv,label='-'+k)
    pl.legend(loc=0,frameon=False)
    if x is not None:
      pl.xlabel(x)
    pl.ylabel(ylab)
  def plot_ipbeta(self,n1=None,n2=None,x=None):
    if x is None:
      xv=arange(len(self[self.keys()[0]][n1:n2]))
    else:
      xv=self[x][n1:n2]
    for k in self.get_vars('bet'):
      kv=self[k][n1:n2]
      pl.plot(xv,kv,label=k)
    pl.legend(loc=0,frameon=False)
    if x is not None:
      pl.xlabel(x)
    pl.ylabel(r"$\beta$ [m]")
  def plot_phase(self,n1=None,n2=None,x=None):
    if x is None:
      xv=arange(len(self[self.keys()[0]][n1:n2]))
    else:
      xv=self[x][n1:n2]
    colors='bgrc'
    for k in self.get_phases():
      kv=self[k][n1:n2]
      pl.plot(xv,kv,color=colors[0],label=k)
      colors=colors[1:]+colors[0]
    pl.legend(loc=0,frameon=False)
    if x is not None:
      pl.xlabel(x)
    pl.ylabel(r'mu [2$\pi$]')
  def plot_squeeze(self,n1=0,n2=None,x=None,title='squeeze'):
    fig=pl.figure(title,figsize=(16,12))
    fig.canvas.mpl_connect('button_release_event',self.button_press)
    pl.clf()
    if len(self.get_vars('kqx'))>0:
      pl.subplot(3,4,1)
      self.plot_triplet(n1,n2,x=x)
    for n in range(4,11):
      if len(self.get_kq(n))>0:
        pl.subplot(3,4,n-2)
        self.plot_2in1(n,n1,n2,x=x,sign=False)
    for n in range(11,14):
      if len(self.get_kq(n))>0:
        pl.subplot(3,4,n-2)
        self.plot_2in1(n,n1,n2,x=x,sign=True)
    #pl.subplot(3,4,12)
    #self.plot_ipbeta(n1,n2,x=x)
    pl.tight_layout()
    self.xvar=x
    return self
  def _plot_beta(self,reg,n1=0,n2=None,lim=None,s=(3,4,1),lbl=''):
    """plot beta function for squeeze points
    *n1* to *n2* for regular expression *reg*
    and in subplot *s*, where *s* is a tuple"""
    if len(self.get_vars(reg))>0:
      s1,s2,s3=s
      pl.subplot(s1,s2,s3)
      pl.cla()
      if lim!=None:
        for ll in lim:
          if n2==None: n2=len(self[self.get_vars(reg)[0]])#n2=tablelength
          pl.plot([n1,n2],[ll,ll],'k-')
      for k in self.get_vars(reg):
        pl.plot(self[k][n1:n2],label=k)
      pl.ylabel(lbl)
      pl.legend(loc='best')
  def plot_squeeze_ir6(self,n1=0,n2=None,x=None,title='squeeze'):
    self.plot_squeeze(n1=n1,n2=n2,x=x,title=title)
    self._plot_beta('b[xy]dumpb[12]',n1=n1,n2=n2,s=(3,4,1),lim=[5012,3955,3698],lbl=r'$\beta_{\rm dump}$ [m]')
#    self._plot_beta('betxmkdb[12]',n1=n1,n2=n2,lim=[380],s=(3,4,4),lbl=r'$\beta_{\rm MKD}$ [m]')
    self._plot_beta('dmuxkickb[12]',n1=n1,n2=n2,lim=[0.25],s=(3,4,4),lbl=r'$\Delta\mu_{\rm (MKD.H -> TCSG)}$ [m]')
    self._plot_beta('bet[xy]tcdqb[12]',n1=n1,n2=n2,lim=[100],s=(3,4,5),lbl=r'$\beta_{\rm TCDQ}$ [m]')#lim=[160,495]
    pl.ylim(ymin=90)
    self._plot_beta('bet[xy]_ip[15]',n1=n1,n2=n2,s=(3,4,12),lbl=r'$\beta_{\rm IP[15]}$ [m]')
  def poly_load(self,fn='',pol=None):
    """return dictionary with dic[var]=p , where
    var is the variable and p the polynomial defined as
    p=[p[0],p[1],...,p[n]]"""
    p={}
    if os.path.isfile(fn):
      for s in open(fn).read().split('\n'):
        if s.find(':=')>-1:
          k,o=self.poly_def(s)
          p[k]=o
    if pol!=None:
      if pol.find(':=')>-1:
        k,o=self.poly_def(pol)
        p[k]=o
    return p
  def plot_dump(self,n1=0,n2=None,title='Twiss @ Dump',fn=None):
    """plot the beta function at the dump. *fn* can be given to
    plot in addition the beta function at the dump from the fit,
    where fn contains:
    betxip6b1:=+2.4683372155e+02*k^0-8.8694921897e-01*k^1;
    alfxip6b1:=-1.0363147736e+00*k^0+7.4283762955e-03*k^1;
    ...
    """
    cls={'x1':'b','x2':'g','y1':'r','y2':'c'}
    p={}
    if os.path.isfile(fn):
      p=self.poly_load(fn)
    for k in ['x','y']:
      for b in ['1','2']:
        bet,alf=[ self.get_vars('%s%sip.b%s'%(a,k,b))[0] for a in ['bet','alf'] ]
        pl.plot(bdump(self[bet],self[alf],b),'%s-'%cls[k+b],label='b%sdumpb%s'%(k,b))
        if os.path.isfile(fn):
          x=arange(len(self[bet]))
          pl.plot(bdump(poly_val(p[bet],x),poly_val(p[alf],x),b),'%s--'%cls[k+b],label='b%sdumpb%s fit'%(k,b))
    for ll in [5012,3955,3698]:
      if n2==None: n2=len(self[self.get_vars('betxip.b1')[0]])#n2=tablelength
      pl.plot([n1,n2],[ll,ll],'k-')
    pl.legend() 
  def plot_betip(self,n1=0,n2=None):
    fig=pl.figure(figsize=(10,8))
    pl.subplot(2,2,1)
    for k in self.get_vars('bet[xy]ip.b[12]'):
      pl.plot(self[k][n1:n2],label=k)
      pl.ylabel(r'$\beta_{\rm IP}$ [m]')
      pl.legend(loc='best')
    pl.subplot(2,2,2)
    for k in self.get_vars('alf[xy]ip.b[12]'):
      pl.plot(self[k][n1:n2],label=k)
      pl.ylabel(r'$\alpha_{\rm IP}$ [m]')
      pl.legend(loc='best')
    pl.subplot(2,2,3)
    for k in self.get_vars('dxip.b[12]'):
      pl.plot(self[k][n1:n2],label=k)
      pl.ylabel(r'$D_{x,\rm IP}$ [m]')
      pl.legend(loc='best')
    pl.subplot(2,2,4)
    for k in self.get_vars('dpxip.b[12]'):
      pl.plot(self[k][n1:n2],label=k)
      pl.ylabel(r"$D'_{x,\rm IP}$ [m]")
      pl.legend(loc='best')
    pl.tight_layout()
  def plot_betsqueeze(self,n1=0,n2=None,figname=None):
    x=self.get_vars('betxip')[0]
    if figname is None:
      fig=pl.figure(x,figsize=(16,12))
    else:
      fig=pl.figure(figname,figsize=(16,12))
    fig.canvas.mpl_connect('button_release_event',self.button_press)
    pl.clf()
    pl.subplot(3,4,1)
    if len(self.get_vars('kqx'))>0:
      self.plot_triplet(n1,n2,x=x)
    for n in range(4,11):
      pl.subplot(3,4,n-2)
      self.plot_2in1(n,n1,n2,x=x,sign=False)
    for n in range(11,14):
      pl.subplot(3,4,n-2)
      self.plot_2in1(n,n1,n2,x=x,sign=True)
    pl.subplot(3,4,12)
    self.plot_phase(n1,n2,x=x)
    pl.tight_layout()
    self.xvar=x
    return self
  def plot_knobs(self,n1=0,n2=None,figname=None,scales=[1,1]):
    x=self.get_vars('betxip')[0]
    if figname is None:
      fig=pl.figure('knobs',figsize=(16,12))
    else:
      fig=pl.figure(figname,figsize=(16,12))
    fig.canvas.mpl_connect('button_release_event',self.button_press)
    pl.clf()
    for ii,(knob,scale) in enumerate(zip(['on_x','on_sep'],scales)):
       pl.subplot(2,3,1+ii*3)
       pl.title('%s=%g'%(knob,scale))
       self.plot_acb(1,knob,n1,n2,x=x,scale=scale)
       self.plot_acb(2,knob,n1,n2,x=x,scale=scale)
       self.plot_acb(3,knob,n1,n2,x=x,scale=scale)
       pl.ylim(-100,100)
       pl.subplot(2,3,2+ii*3)
       pl.title('%s=%g'%(knob,scale))
       self.plot_acb(4,knob,n1,n2,x=x,scale=scale)
       pl.ylim(-100,100)
       pl.subplot(2,3,3+ii*3)
       pl.title('%s=%g'%(knob,scale))
       self.plot_acb(5,knob,n1,n2,x=x,scale=scale)
       if self.get_acb(6):
         self.plot_acb(6,knob,n1,n2,x=x,scale=scale)
       pl.ylim(-100,100)
    pl.tight_layout()
    self.xvar=x
    return self
  def linint(t,name,val1,val2,x,par='ttt'):
    n1=where(t[x]==val1)[0][-1]
    n2=where(t[x]==val2)[0][-1]
    tmp= "%s:=(%g)*(1-%s)+(%g)*(%s);"
    print tmp%(name,t[name][n1],par,t[name][n2],par)
  def button_press(self,event):
    self.event=event
  def poly_def(self,s):
    """converts var:=p[0]*k^0+...+p[n]*k^n to
    an array p=[p[0],...,p[n]] and return (var,p)"""
    if s.find(':=')>-1:
      var,p=s.split(':=')
      return var,[ float(ss) for ss in re.sub('.k.[0-9]',' ',p).split(' ')[:-1] ]
  def poly_fit(self,var='betxip.b1',order=2,n1=0,n2=None,param=None,fn=None,force=False):
    """plot *param* vs *var* and the polynomial fit
    If *fn* is None, a new polynomial of order *order*
    is fitted with fixed start and endpoints at *n1*
    *n2*. 
    If *fn* is given the polynomial defined in *fn* is
    plotted."""
    if param==None: param_plot='k'
    else: param_plot=param
    if fn!=None:
      p=self.poly_load(fn)
      var=p.keys()
    var=iter(var)
    lsub=def_subplot(len(var))#define the grid of subplots
    if fn==None or force==True:
      out=[]
    else:
      out=p
    for v,s in zip(var,lsub):
      mk_subplot(s)
      if param==None:
        if n2==None: n2=len(self[v])
        x=arange(n1,n2)
      y=self[v]
      pl.plot(x,y,'r')
      pl.xlabel(param_plot)
      pl.ylabel(v)
      if fn==None or force==True:
        pol=self.poly_fit_var(var=v,order=order,n1=n1,n2=n2,param=param)
        out.append("%s:=%s;"%(v, poly_print(self.poly_load(pol=pol)[v],x=param_plot,power='^')))
      else:
        pl.plot(x,poly_val(p[v],x),'k--')
    pl.tight_layout()
    return out
  def poly_fit_var(self,var,order,n1=0,n2=None,param=None,x0_idx=[0,-1],y0_idx=[],xp0_idx=[],yp0_idx=[]):
    """makes a polynomial fit of order *order*
    to parameter *param* vs *var* between [n1,n2].
    Parameters
    ---------
    x0_idx: indices of fixed points
    y0_idx: value at x0_idx. if y0_idx=[], 
            If y0_idx=[] is set to the
            y-value of *var* for index x0_idx:
              y0_idx=var(x0_idx)
    xp0_idx,yp0_idx: as x0_idx and y0_idx,
            only that the slope is fixed.
            yp0_idx=[]: yp0_idx is set to the
              the slope of *var* for index x0_idx:
                yp0_idx=var(x0_idx+1)-var(x0_idx)
            yp0_idx=k=const.: yp0_idx is set to
              the constant value k. E.g. for dervative
              0, set yp0_idx=0
    """
    if n2==None or n2==-1: n2=len(self[var])
    if param==None:
      x=arange(n1,n2)
      param='k'#set some dummy parameter name to print later the polynomial
    else: x=self[param][n1:n2];
    y=self[var][n1:n2]
    x0=[];y0=[];xp0=[];yp0=[]
    for idx in iter(x0_idx):
      if idx==-1: idx=-1
      else: idx=idx-n1#shift by n1
      x0.append(x[idx])
      if len(y0_idx)==0:
        y0.append(y[idx])
    for idx in iter(xp0_idx):
      idx=idx-n1#shift by n1
      xp0.append(x[idx])
      if len(yp0_idx)==0:
        if idx==n2-n1-1:#endpoint
          yp0.append(y[idx]-y[idx-1])
        else:#all other points
          yp0.append(y[idx+1]-y[idx])
      if len(yp0_idx)==1:
        yp0.append(yp0_idx)
    pol=poly_fit(order,x,y,x0,y0,xp0,yp0)
    out="%s:=%s;"%(var, poly_print(pol,x=param,power='^'))
    yv=poly_val(pol,x)
    pl.plot(x,yv,'k--')
    return out
  def list_betip(self):
    t=['bet','alf']; d=['dx','dpx']; p=['x','y']; q=['1','2']
    tip=['%s%sip.b%s'%(a,b,c) for a in t for b in p for c in q]
    dip  =['%sip.b%s'%(a,c) for a in d for c in q]
    return tip+dip
  def poly_fit_betip(self,order=1,n1=0,n2=None,param=None,x0_idx=[0,-1],xp0_idx=[],fn='fit.out',force=False):
    """make a polinomial fit to the twiss parameters at the ip
    if force=false the given file fit.out is used, if force=true
    a new fit with a polynomial of order *order* is performed"""
    self.plot_betip(n1=n1,n2=n2)
    out=[]
    if os.path.isfile(fn) and force==False:
      pol=self.poly_load(fn)
      out=open(fn,'r')
    for s,p in zip([1,2,3,4],['bet[xy]ip.b[12]','alf[xy]ip.b[12]','dxip.b[12]','dpxip.b[12]']):
      pl.subplot(2,2,s)
      for k in self.get_vars(p):
        if os.path.isfile(fn) and force==False:
          if k in pol.keys():
            x=arange(len(self[k]))
            pl.plot(x,poly_val(pol[k],x),'k--')
        else:
          out.append(self.poly_fit_var(k,order,param=param,n1=n1,n2=n2,x0_idx=x0_idx,xp0_idx=xp0_idx))
    if os.path.isfile(fn) and force==False:
      out=open(fn).read().split('\n')
    else:
      open(fn,'w').write('\n'.join(out))
    return out
  def poly_fit_k(self,var,order,n1=None,n2=None,param="betxip8b1",slope0=[]):
    scale=StrTable.scale
    x=self[param][n1:n2]; y=self[var][n1:n2]
    x0=[x[0],x[-1]]
    y0=[y[0],y[-1]]
    xp0=[]; yp0=[]
    for idx in slope0:
        xp0.append(x[idx])
        yp0.append(0)
    pol=poly_fit(order,x,y,x0,y0,xp0,yp0)
    out="%s:=%s;"%(var, poly_print(pol,x=param,power='^'))
    n=re.match('[kqtxl]+([0-9]+)\.',var)
    if n is None:
        if re.match('kqt?x\.',var):
            n=3
        else:
            n=14;scale=1
    else:
        n=int(n.groups()[0])
    pl.subplot(3,4,n-2)
    yv=poly_val(pol,x)
    print ' '.join(['%2d'%i for i in  sign(diff(yv))])
    if n<11:
        yv=abs(yv)
    pl.plot(self[self.xvar],yv*scale)
    return out
  def poly_fit_all(self,order,param,fn):
    out=[]
    for kq in self.get_triplet():
        out.append(self.poly_fit_k(kq,order,param=param))
    for n in range(4,14):
        for kq in self.get_kq(n):
            out.append(self.poly_fit_k(kq,order,param=param))
    open(fn,'w').write('\n'.join(out))
  def check_slopes(self):
    for kq in range(4,14):
        for name in self.get_kq(kq):
          v=self[name]
          slopes=sign(diff(v))
          if sum(abs(diff(slopes)))>0:
             print name,' '.join(['%2d'%i for i in  slopes])
  def set_log(self):
    fig=pl.gcf()
    for ax in fig.axes:
        ax.set_xscale('log')
    pl.draw()
    return self
  def set_xlim(self,a,b):
    fig=pl.gcf()
    for ax in fig.axes:
        ax.set_xlim(a,b)
    pl.draw()
    return self
  def savefig(self):
    self.plot_betsqueeze()
    pl.savefig(self.filename.replace('.tfs','.png'))
    if len(self.get_acb(1))>0:
      self.plot_knobs()
      pl.savefig(self.filename.replace('.tfs','_knobs.png'))
    return self



