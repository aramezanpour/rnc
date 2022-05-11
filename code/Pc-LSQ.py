from numpy import *           
from random import uniform,randint           
from scipy.optimize import least_squares


#################

def make_I(N):
  
  r=0
  for i in range(N):
    I[r][0]=i
    r+=1
    
#################

def make_II(N):
  
  L=int(N*(N-1)/2)
  Li=zeros(L)
  Lj=zeros(L)
  
  ###
  
  l=0
  for i in range(N):
    for j in range(i):
        Li[l]=i
        Lj[l]=j
        l+=1  
  ###
  
  r=0
  while(r<M):
    ll=randint(0,l-1)
    i=Li[ll]
    j=Lj[ll]
    l=l-1
    Li[ll]=Li[l]
    Lj[ll]=Lj[l]
    I[r][0]=i
    I[r][1]=j
    r+=1
  
  
#################

def make_R(P):
  
  for a in range(P):
    R[a]=uniform(0,1)/3
    if(uniform(0,1)<0.5):
      R[a]=2.0/3+uniform(0,1)/3
      
  
#################

def make_S(P,N):
  
  for a in range(P):
    for i in range(N):
      S[a][i]=uniform(0,1)/3
      if(uniform(0,1)<0.5):
        S[a][i]=2.0/3+uniform(0,1)/3
  
#################

def make_A(P,M,L):
  
  for a in range(P):
    for r in range(M):
      
      prod=1
      for l in range(L):
	i=I[r][l]
	prod=prod*S[a][i]
      
      A[a][r]=prod
  
#################

def fun(x):
  
  y=zeros(P)
  for a in range(P):
    
    summ=0
    for r in range(M):      
      summ+=A[a][r]*x[r]
    summ-=R[a]
    
    y[a]=summ
  
  return y

#################

def check(P,M):
  
  sat=True
  for r in range(M):
      if(k[r]<0):sat=False
  
  for a in range(P):
    
    summ=0
    for r in range(M):
      summ+=A[a][r]*k[r]
      
    if(summ>1.0):sat=False  
    if((R[a]<1.0/3)and(summ>1.0/3)):sat=False  
    if((R[a]>2.0/3)and(summ<2.0/3)):sat=False  
      
      
  return sat

#################

Nsample=10000

N=40
M=N
L=2

output=open("LS.dat",'w')
output.close()

P=12
while(P<13):
  
  A=zeros((P,M),float)
  S=zeros((P,N),float)
  I=zeros((M,L),int)
  R=zeros(P,float)
  
  make_II(N)
  
  K1=0.0
  K2=0.0
  K3=0.0
  K4=0.0
  ps=0.0
  for sample in range(Nsample):
    
      make_R(P)
      make_S(P,N)
      make_A(P,M,L)
    
      try:
        x0=[1]*M
        LS=least_squares(fun,x0,bounds=(0.0001,10))
        k=LS.x
      except:
        k=[-1]*M
        pass
    
    
      if(check(P,M)):
          ps+=1
          sum1=0.0
          sum2=0.0
          kmax=-1
          for r in range(M):
            sum1+=k[r]
            sum2+=k[r]*k[r]
            if(k[r]>kmax):kmax=k[r]
          sum1=sum1/M
          sum2=sum2/M
          K1+=sum1
          K2+=sqrt(sum2-sum1*sum1)
          K3+=kmax
          K4+=kmax*kmax
      
      q1=0
      q2=0
      q3=0
      q4=0
      qs=0
      if(ps>0):     
        q1=K1/ps
        q2=K2/ps
        q3=K3/ps
        q4=K4/ps
        q4=sqrt((q4-q3*q3)/ps)
      qs=ps/(sample+1)
  
  
      print P,qs,q1,q2,q3,q4  
  
      output=open("LS.dat",'a')
      output.write(str(P)+" "+str(qs)+" "+str(q1)+" "+str(q2)+" "+str(q3)+" "+str(q4)+" "+str(sample)+"\n")
      output.close()
  
  P+=1