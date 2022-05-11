#include <parallel/algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <math.h>
#include <vector>
#include <time.h>
#include <string>
using namespace std;


#define pow2(n) (1 << (n))

#define bit_map(x, n) (((x) & (1 << (n))) != 0)

#define bx(x)        ((x) - (((x)>>1) & 0x77777777)  \
			   - (((x)>>2) & 0x33333333)  \
	 		   - (((x)>>3) & 0x11111111))
#define count_bits(x)  (((bx(x)+(bx(x)>>4)) & 0x0F0F0F0F) % 255)

#define hdist(x, y) (count_bits((x)^(y)))


////////////////////// functions


void alloc_mem();

void make_i();
void make_O();

void make_R11();
void make_R21();
void make_R12();
void make_R22();

void make_rev();

void initial();

void initial_X(int p);

void update_X();

void compute_ave();

void check_response();

void report(int sample);

double rang(int& seed);


////////////////////// end functions



///////////////////// variables

const int  NS=6;
const int  NR=1;
const int  N=NS+NR;
const int  ND=0;
const int  K=2;
const int  Nmax=1000;
const int  Pat=8;

const int  TL=2000;
const int  To=2000;
const int  Tav=1000;
const int  Nav=10;

const double  Cmin=0.01;
const double  Cmax=10.0;
const double  eta=0.2*(Cmax/Pat);

//int  Mir=NS*NR;
int  Mir=NR*(NS*(NS-1))/2;
//int  Mir=NS*(NR*(NR-1))/2;
//int  Mir=((NS*(NS-1))/2)*((NR*(NR-1))/2);
int  M=2*Mir;
int  kmax=M;

const double  pr=1.0;

int seed;

double X0;
double Et,Emin;
double Pt,Pmax;

vector<double> avE1,avE2,avP1,avP2;


vector<int>  O;
vector<int>  Dr;
vector<int>  Xi;
vector<double> Xr;

vector<short int>  L,R;
vector<short int>  Di,Si,Ri;
vector<vector<short int> >  Sig,Res;

vector<int>  ki;
vector<vector<int> >  ri;

vector<double>  Cr,Copt;
vector<vector<int> >  Il,Ir;
vector<vector<short int> >  vl,vr;

vector<double>  avX,Var;
vector<vector<double> > CI;


////////////////////// end variables





////////////////////// main

int main()
{

int nav,p;
int i1,i2,i3,i4,id,ir,is;
double dC;


seed=time(NULL);

alloc_mem();

make_i();

for(int t=0;t<TL;t++){
avE1[t]=0;
avE2[t]=0;
avP1[t]=0;
avP2[t]=0;
}

int Nsample=10000;
for(int sample=0;sample<Nsample;sample++){


//make_R11();
make_R21();
//make_R12();
//make_R22();

make_rev();

initial();

/////////////////////////////////////// optimization

Emin=1e+6;
Pmax=0;
for(int r=0;r<M;r++)Copt[r]=Cr[r];


for(int t=0;t<TL;t++){

if(rang(seed)<0.25){  
for(int r=0;r<M;r++){
    if(rang(seed)<0.5){
      Cr[r]=0.5*Cr[r]+0.5*Copt[r];
    }else{
      Cr[r]=0.5*Cr[r]+0.5*(Cmin+rang(seed)*(Cmax-Cmin));
    }  
}
}


make_O();

for(int lp=0;lp<Pat;lp++){
  p=O[lp];
  
  initial_X(p);

  for(int i=0;i<N;i++){
    avX[i]=0;  
    Var[i]=0;
    for(int j=0;j<i+1;j++)CI[i][j]=0; 
  }
    
///
  
  nav=0;
  for(int to=0;to<To;to++){ 
    
    for(int ls=0;ls<NS;ls++){
      is=Si[ls];
      if(Sig[p][is]<0)Xi[is]=min(Xi[is],Nmax/3);	
      if(Sig[p][is]>0)Xi[is]=max(Xi[is],2*Nmax/3);	
    }    
    for(int ld=0;ld<ND;ld++){
      id=Di[ld];
      Xi[id]=Dr[id];	
    }
        
    
    for(int r=0;r<M;r++)update_X();
        
    
    if(to>Tav){
      nav+=1;
      compute_ave();
    }
    
  }

///

  for(int i=0;i<N;i++){
    avX[i]=avX[i]/nav;
    Var[i]=Var[i]/nav;
    Var[i]=sqrt(Var[i]-avX[i]*avX[i])+1;
    for(int j=0;j<i+1;j++)CI[i][j]=CI[i][j]/nav;
  }
  
  
///
  
  for(int r=0;r<M;r++){

    dC=0;
    for(int l=0;l<R[r];l++){
      i1=Ir[r][l];
      for(int ll=0;ll<NR;ll++){
        i2=Ri[ll];
        i3=max(i1,i2);
        i4=min(i1,i2);
        dC+=Res[p][i2]*(CI[i3][i4]-avX[i3]*avX[i4])/(Var[i3]*Var[i4]);
      }
    }
    Cr[r]+=eta*dC;
    if(Cr[r]<Cmin)Cr[r]=Cmin;
    if(Cr[r]>Cmax)Cr[r]=Cmax;
    
  }
  
    
}

//////////////


Et=0;
Pt=0;
for(int na=0;na<Nav;na++)check_response();
Et=Et/Nav;
Pt=Pt/Nav;


Et=Et/Pat;
Pt=Pt/Pat;
if(Et<Emin){
  Emin=Et;
  Pmax=Pt;
  for(int r=0;r<M;r++)Copt[r]=Cr[r];
}
avE1[t]+=Emin;
avE2[t]+=Emin*Emin;
avP1[t]+=Pmax;
avP2[t]+=Pmax*Pmax;


}
//cout<<sample<<endl;

report(sample);
}    



 return 0;
}
////////////////////// end  main




/////////////////// alloc_mem

void alloc_mem()
{

 
O.resize(Pat);

Dr.resize(ND);
Di.resize(ND);
Si.resize(NS);
Ri.resize(NR);
  

Sig.resize(Pat);
Res.resize(Pat);
for(int p=0;p<Pat;p++){
Sig[p].resize(N);
Res[p].resize(N);
}


Xi.resize(N);
Xr.resize(M);


Cr.resize(M);
Copt.resize(M);


L.resize(M);
R.resize(M);

Il.resize(M);
Ir.resize(M);
vl.resize(M);
vr.resize(M);
for(int r=0;r<M;r++){
Il[r].resize(K);
Ir[r].resize(K);
vl[r].resize(K);
vr[r].resize(K);
}


ki.resize(N);
ri.resize(N);
for(int i=0;i<N;i++){
ri[i].resize(kmax);
}

avX.resize(N);
Var.resize(N);

CI.resize(N);
for(int i=0;i<N;i++){
CI[i].resize(i+1);
}


avE1.resize(TL);
avE2.resize(TL);
avP1.resize(TL);
avP2.resize(TL);




}
/////////////////// end  alloc_mem




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




////////////////////// make_i

void make_i()
{
  
for(int ls=0;ls<NS;ls++)Si[ls]=ls;
for(int lr=0;lr<NR;lr++)Ri[lr]=NS+lr;
for(int ld=0;ld<ND;ld++)Di[ld]=rang(seed)*NS;
 
}
////////////////////// end make_i



////////////////////// make_O

void make_O()
{
  
int n,m,l;
int Lp[Pat];

n=Pat;
for(int p=0;p<Pat;p++)Lp[p]=p;


m=0;
while(n>0){
l=rang(seed)*n;
O[m]=Lp[l];
n-=1;
Lp[l]=Lp[n];
m+=1;
}



}
////////////////////// end make_O




////////////////////// make_R11

void make_R11()
{

int i1,i2,li;  
int r;


r=0;
for(int j=0;j<N;j++)ki[j]=0;


///////// layer 1

for(int l1=0;l1<NS;l1++){
for(int l2=0;l2<NR;l2++){
if(rang(seed)<pr){  
i1=Si[l1];
i2=Ri[l2];

L[r]=1;
R[r]=1;
Cr[r]=Cmin+rang(seed)*(Cmax-Cmin);

Il[r][0]=i1;
vl[r][0]=-1;
li=ki[i1];
ri[i1][li]=r;
ki[i1] += 1;

Ir[r][0]=i2;
vr[r][0]=+1;

r+=1;
}
}
}
Mir=r;
M=r;


}
////////////////////// end make_R11




////////////////////// make_R21

void make_R21()
{

int i1,i2,i3,i4,li;  
int r;


r=0;
for(int j=0;j<N;j++)ki[j]=0;


///////// layer 1

for(int l1=0;l1<NS;l1++){
for(int l2=0;l2<l1;l2++){
for(int l3=0;l3<NR;l3++){
if(rang(seed)<pr){  
i1=Si[l1];
i2=Si[l2];
i3=Ri[l3];

L[r]=2;
R[r]=1;
Cr[r]=Cmin+rang(seed)*(Cmax-Cmin);

Il[r][0]=i1;
vl[r][0]=-1;
li=ki[i1];
ri[i1][li]=r;
ki[i1] += 1;

Il[r][1]=i2;
vl[r][1]=-1;
li=ki[i2];
ri[i2][li]=r;
ki[i2] += 1;

Ir[r][0]=i3;
vr[r][0]=+1;

r+=1;
}
}
}
}
Mir=r;
M=r;


}
////////////////////// end make_R21




////////////////////// make_R12

void make_R12()
{

int i1,i2,i3,i4,li;  
int r;


r=0;
for(int j=0;j<N;j++)ki[j]=0;


///////// layer 1

for(int l1=0;l1<NS;l1++){
for(int l2=0;l2<NR;l2++){
for(int l3=0;l3<l2;l3++){
if(rang(seed)<pr){  
i1=Si[l1];
i2=Ri[l2];
i3=Ri[l3];

L[r]=1;
R[r]=2;
Cr[r]=Cmin+rang(seed)*(Cmax-Cmin);

Il[r][0]=i1;
vl[r][0]=-1;
li=ki[i1];
ri[i1][li]=r;
ki[i1] += 1;

Ir[r][0]=i2;
vr[r][0]=+1;

Ir[r][1]=i3;
vr[r][1]=+1;

r+=1;
}
}
}
}
Mir=r;
M=r;


}
////////////////////// end make_R12




////////////////////// make_R22

void make_R22()
{

int i1,i2,i3,i4,li;  
int r;


r=0;
for(int j=0;j<N;j++)ki[j]=0;


///////// layer 1

for(int l1=0;l1<NS;l1++){
for(int l2=0;l2<l1;l2++){
for(int l3=0;l3<NR;l3++){
for(int l4=0;l4<l3;l4++){
if(rang(seed)<pr){  
i1=Si[l1];
i2=Si[l2];
i3=Ri[l3];
i4=Ri[l4];

L[r]=2;
R[r]=2;
Cr[r]=Cmin+rang(seed)*(Cmax-Cmin);

Il[r][0]=i1;
vl[r][0]=-1;
li=ki[i1];
ri[i1][li]=r;
ki[i1] += 1;

Il[r][1]=i2;
vl[r][1]=-1;
li=ki[i2];
ri[i2][li]=r;
ki[i2] += 1;

Ir[r][0]=i3;
vr[r][0]=+1;

Ir[r][1]=i4;
vr[r][1]=+1;

r+=1;
}
}
}
}
}
Mir=r;
M=r;


}
////////////////////// end make_R22




////////////////////// make_rev

void make_rev()
{
  
int i,li,rr;
 

for(int r=0;r<Mir;r++){
  
  rr=Mir+r;
  
  L[rr]=R[r];
  R[rr]=L[r];
  Cr[rr]=Cmin+rang(seed)*(Cmax-Cmin);
  
  ///
  
  for(int l=0;l<R[r];l++){
        i=Ir[r][l];
        Il[rr][l]=i;
        vl[rr][l]=-1;
        li=ki[i];
        ri[i][li]=rr;
        ki[i] += 1;
  }
  
  for(int l=0;l<L[r];l++){
       i=Il[r][l];
       Ir[rr][l]=i;
       vr[rr][l]=+1;
  }
  
  ///
  
}
M=2*Mir;


  
}
////////////////////// end make_rev




////////////////////// initial

void initial()
{

int id,is,ir;  
int l,np,nq,pp,qq;
int Pm=(1 << NS);
int Qm=(1 << NR);
int Lp[Pm];
short int check;


for(int ld=0;ld<ND;ld++){
  id=Di[ld];
  Dr[id]=rang(seed)*Nmax;
}

//////////

np=Pm;
for(int p=0;p<Pm;p++)Lp[p]=p;

///

for(int p=0;p<Pat;p++){

  
for(int i=0;i<N;i++){
  Sig[p][i]=0;
  Res[p][i]=0;
}

///

l=rang(seed)*np;
pp=Lp[l];
np=np-1;
Lp[l]=Lp[np];
qq=rang(seed)*Qm;

///

for(int ls=0;ls<NS;ls++){
  is=Si[ls];
  Sig[p][is]=2*bit_map(pp,ls)-1;  
}
for(int lr=0;lr<NR;lr++){
  ir=Ri[lr];
  Res[p][ir]=2*bit_map(qq,lr)-1;  
}



}



}
////////////////////// end initial




////////////////////// initial_X

void initial_X(int p)
{

int i,is;    
double prod;



for(int i=0;i<N;i++){
     Xi[i]=rang(seed)*Nmax; 
}

for(int ls=0;ls<NS;ls++){
  is=Si[ls];
  if(Sig[p][is]<0)Xi[is]=2+rang(seed)*Nmax/3.0;	
  if(Sig[p][is]>0)Xi[is]=(2+rang(seed))*Nmax/3.0;	
}    

///

X0=0;
for(int r=0;r<M;r++){
    prod=1.0;
    for(int l=0;l<L[r];l++){
      i=Il[r][l];
      prod=(prod*Xi[i])/Nmax;
    }
    Xr[r]=Cr[r]*prod;
    X0+=Xr[r];  
}
 
///


 
}
////////////////////// end initial_X





///////////////////  update_X
void update_X()
{

int i,io,r,ro;    
double u;
double sum,prod;
short int checkr,checkx;

u=rang(seed);

r=0;
sum=0;
checkr=0;
while((checkr==0)&&(r<M)){
   sum+= Xr[r];
   if(sum>u*X0){
     checkr=1;
   }else{  
     r+=1;
   }
}

    
if(checkr==1){
checkx=1;
for(int l=0;l<L[r];l++){
  i=Il[r][l];  
  if(Xi[i]+vl[r][l]<0)checkx=0;
}
}

if((checkr==1)&&(checkx==1)){


///

for(int l=0;l<L[r];l++){
  i=Il[r][l];  
  Xi[i]=Xi[i]+vl[r][l];
}
for(int l=0;l<R[r];l++){
  i=Ir[r][l];  
  Xi[i]=min(Nmax-1,Xi[i]+vr[r][l]);
}

///

for(int l=0;l<L[r];l++){
  i=Il[r][l];
  for(int li=0;li<ki[i];li++){
    ro=ri[i][li];
    prod=1.0;
    for(int lo=0;lo<L[ro];lo++){
      io=Il[ro][lo];  
      prod=(prod*Xi[io])/Nmax;  
    }
    X0=X0-Xr[ro];
    Xr[ro]=Cr[ro]*prod;
    X0=X0+Xr[ro];
  }
}
for(int l=0;l<R[r];l++){
  i=Ir[r][l];
  for(int li=0;li<ki[i];li++){
    ro=ri[i][li];
    prod=1.0;
    for(int lo=0;lo<L[ro];lo++){
      io=Il[ro][lo];  
      prod=(prod*Xi[io])/Nmax;  
    }
    X0=X0-Xr[ro];
    Xr[ro]=Cr[ro]*prod;
    X0=X0+Xr[ro];
  }
}

///

}



}
/////////////////// update_X




/////////////////// compute_ave
void compute_ave()
{


for(int i=0;i<N;i++){
  avX[i]+=Xi[i];
  Var[i]+=Xi[i]*Xi[i];
  for(int j=0;j<i+1;j++)CI[i][j]+=Xi[i]*Xi[j];    
}


}
/////////////////// compute_ave



  
////////////////////// check_response

void check_response()
{

int id,is;


for(int p=0;p<Pat;p++){

  initial_X(p);

  for(int i=0;i<N;i++){
    avX[i]=0;  
    Var[i]=0;
    for(int j=0;j<i+1;j++)CI[i][j]=0; 
  }
  
///
  
  int nav=0;
  for(int to=0;to<To;to++){ 
    
    for(int ls=0;ls<NS;ls++){
      is=Si[ls];
      if(Sig[p][is]<0)Xi[is]=min(Xi[is],Nmax/3);	
      if(Sig[p][is]>0)Xi[is]=max(Xi[is],2*Nmax/3);	
    }    
    for(int ld=0;ld<ND;ld++){
      id=Di[ld];
      Xi[id]=Dr[id];	
    }    
    
    for(int r=0;r<M;r++)update_X();
    
    if(to>Tav){
      nav+=1;
      compute_ave();
    }
    
  }

///

  int E=0; 
  for(int i=0;i<N;i++){
    avX[i]=avX[i]/nav;
    Var[i]=Var[i]/nav;
    Var[i]=sqrt(Var[i]-avX[i]*avX[i]);
    for(int j=0;j<i+1;j++)CI[i][j]=CI[i][j]/nav;

    if(Res[p][i]>0){
      Et+=1-avX[i]/Nmax;
      if(avX[i]<2*Nmax/3.0)E+=1;
    }
    if(Res[p][i]<0){
      Et+=avX[i]/Nmax;
      if(avX[i]>Nmax/3.0)E+=1;  
    }
  }
  if(E==0)Pt+=1;
  
  

}
  

}
////////////////////// end check_response



/////////////////// report
void report(int sample)
{

double e1,e2,p1,p2;

ofstream output("avE-R21-S6R1-P8.dat");
output<<"#sample="<<sample+1<<endl;


for(int t=0;t<TL;t++){
  
e1=avE1[t]/(sample+1);
e2=avE2[t]/(sample+1);
e2=sqrt((e2-e1*e1)/sample);

p1=avP1[t]/(sample+1);
p2=avP2[t]/(sample+1);
p2=sqrt((p2-p1*p1)/sample);


output<< t << "  "<< e1 << "  "<< e2 << "  "<< p1<< "  "<< p2<<endl;
}


}
/////////////////// report




////////////////////// rang

double rang(int& seed)
{

      int a, m, q, r, l;
      double conv, rand;


      a = 16807;
      m = 2147483647;
      q = 127773;
      r = 2836;

      conv = 1.0 / (m - 1);

      l = seed / q;
      seed = a * (seed - q * l) - r * l;
      if (seed < 0) {
	      seed += m;
      }
      rand = conv * (seed - 1);



      return rand;
} 
////////////////////// end rang


