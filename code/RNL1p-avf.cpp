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
void initial();

void make_R11();
void make_R12();
void make_R21();
void make_R22();

void make_RI();
void make_RII();

void make_rev();

void initial_X(int p);

void update_X();

void compute_P(int sample,int p);
void compute_Q();
void compute_ave();

void report(int sample);

double rang(int& seed);


////////////////////// end functions



////////////////////// variables

const int  NS=10;
const int  NI=5;
const int  NR=5;
const int  N=NS+NI+NR;
const int  K=2;
const int  Nmax=1000;
const int  Pat=(1<<NS);

const int  To=5000;
const int  Tav=4000;

const double  Cmin=0.01;
const double  Cmax=10.0;

//int  Mir=NI*NS+NR*NI+(NI*(NI-1))/2;
int  Mir=NS*(NI*(NI-1))/2+NI*(NR*(NR-1))/2+(NI*(NI-1))/2;
//int  Mir=NI*(NS*(NS-1))/2+NR*(NI*(NI-1))/2+(NI*(NI-1))/2;
//int  Mir=(NS*(NS-1)*NI*(NI-1))/4+(NI*(NI-1)*NR*(NR-1))/4+(NI*(NI-1))/2;
int  M=2*Mir;
int  kmax=M;

const double  pS=1.0;
const double  pR=1.0;
const double  pI=1.0;

int seed;

double X0;

vector<int>  Xi;
vector<double> Xr;
vector<vector<double> >  Yr;

vector<short int>  L,R;
vector<short int>  Si,Ri,Hi;
vector<vector<short int> >  Sig;
vector<vector<double> >  Res;

vector<int>  ki;
vector<vector<int> >  ri;

vector<double>  Cr;
vector<vector<int> >  Il,Ir;
vector<vector<short int> >  vl,vr;

vector<double>  avX,Var;
vector<vector<double> > CI;

vector<double>  P1,P2,Q1,Q2,V1,V2;
vector<vector<double> > QQ1,QQ2;

////////////////////// end variables





////////////////////// main

int main()
{

int nav;
int i1,i2,i3,i4,id,ir,is;


seed=time(NULL);

alloc_mem();

make_i();

initial();
  
  
int Nsample=10000;
for(int sample=0;sample<Nsample;sample++){


//make_R11();
make_R12();
//make_R21();
//make_R22();

//make_RI();
//make_RII();

make_rev();


/////////////////////////////////////// evolution


for(int p=0;p<Pat;p++){

  
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
  
  compute_P(sample,p);
  
    
}
//////////////


compute_Q();

cout<<sample<<endl;

report(sample);
}    



 return 0;
}
////////////////////// end  main




/////////////////// alloc_mem

void alloc_mem()
{

 
Si.resize(NS);
Ri.resize(NR);
Hi.resize(NI);
  
P1.resize(Pat);
P2.resize(Pat);
Q1.resize(Pat);
Q2.resize(Pat);
V1.resize(Pat);
V2.resize(Pat);

Yr.resize(Pat);
for(int p=0;p<Pat;p++){
Yr[p].resize(NR);
}

QQ1.resize(Pat);
QQ2.resize(Pat);
for(int p=0;p<Pat;p++){
QQ1[p].resize(p);
QQ2[p].resize(p);
}

Sig.resize(Pat);
for(int p=0;p<Pat;p++){
Sig[p].resize(N);
}


Res.resize(NR);
for(int l=0;l<NR;l++){
Res[l].resize(l+1);
}


Xi.resize(N);
Xr.resize(M);


Cr.resize(M);


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



}
/////////////////// end  alloc_mem




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




////////////////////// make_i

void make_i()
{
  
for(int ls=0;ls<NS;ls++)Si[ls]=ls;
for(int li=0;li<NI;li++)Hi[li]=NS+li;
for(int lr=0;lr<NR;lr++)Ri[lr]=NS+NI+lr;
 
}
////////////////////// end make_i




////////////////////// initial

void initial()
{

int is;  


for(int p=0;p<Pat;p++){
  
for(int ls=0;ls<NS;ls++){
  is=Si[ls];
  Sig[p][is]=2*bit_map(p,ls)-1;  
}

P1[p]=0;
P2[p]=0;
Q1[p]=0;
Q2[p]=0;
V1[p]=0;
V2[p]=0;

for(int lr=0;lr<NR;lr++)Yr[p][lr]=rang(seed)*Nmax;

for(int pp=0;pp<p;pp++){
QQ1[p][pp]=0;  
QQ2[p][pp]=0;  
}

}

for(int l=0;l<NR;l++){
  for(int ll=0;ll<l+1;ll++){
    Res[l][ll]=0;
  }
}
    
 
 
}
////////////////////// end initial




////////////////////// make_R11

void make_R11()
{

int i1,i2,i3,i4,li;  
int r;


r=0;
for(int j=0;j<N;j++)ki[j]=0;


///////// layer 1

for(int l1=0;l1<NS;l1++){
for(int l2=0;l2<NI;l2++){
if(rang(seed)<pS){  
i1=Si[l1];
i2=Hi[l2];

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

///////// layer 2

for(int l1=0;l1<NI;l1++){
for(int l2=0;l2<NR;l2++){
if(rang(seed)<pR){  
i1=Hi[l1];
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




////////////////////// make_R12

void make_R12()
{

int i1,i2,i3,i4,li;  
int r;


r=0;
for(int j=0;j<N;j++)ki[j]=0;


///////// layer 1

for(int l1=0;l1<NS;l1++){
for(int l2=0;l2<NI;l2++){
for(int l3=0;l3<l2;l3++){
if(rang(seed)<pS){  
i1=Si[l1];
i2=Hi[l2];
i3=Hi[l3];

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

///////// layer 2

for(int l1=0;l1<NI;l1++){
for(int l2=0;l2<NR;l2++){
for(int l3=0;l3<l2;l3++){
if(rang(seed)<pR){  
i1=Hi[l1];
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
for(int l3=0;l3<NI;l3++){
if(rang(seed)<pS){  
i1=Si[l1];
i2=Si[l2];
i3=Hi[l3];

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

///////// layer 2

for(int l1=0;l1<NI;l1++){
for(int l2=0;l2<l1;l2++){
for(int l3=0;l3<NR;l3++){
if(rang(seed)<pR){  
i1=Hi[l1];
i2=Hi[l2];
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
for(int l3=0;l3<NI;l3++){
for(int l4=0;l4<l3;l4++){
if(rang(seed)<pS){  
i1=Si[l1];
i2=Si[l2];
i3=Hi[l3];
i4=Hi[l4];

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

///////// layer 2

for(int l1=0;l1<NI;l1++){
for(int l2=0;l2<l1;l2++){
for(int l3=0;l3<NR;l3++){
for(int l4=0;l4<l3;l4++){
if(rang(seed)<pR){  
i1=Hi[l1];
i2=Hi[l2];
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



////////////////////// make_RI

void make_RI()
{

int i1,i2,li;  
int r;


r=Mir;

///////// interactions

for(int l1=0;l1<NI;l1++){
for(int l2=0;l2<l1;l2++){
if(rang(seed)<pI){  
i1=Hi[l1];
i2=Hi[l2];

  
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
////////////////////// end make_RI




////////////////////// make_RII

void make_RII()
{

int i1,i2,li;  
int r;


r=Mir;

///////// interactions

for(int l1=0;l1<NI;l1++){
for(int l2=0;l2<l1;l2++){
if(rang(seed)<pI){  
i1=Hi[l1];
i2=Hi[l2];

if(rang(seed)<0.5){  
  
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

Ir[r][0]=i1;
vr[r][0]=+1;

r+=1;

}else{

   
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

Ir[r][0]=i2;
vr[r][0]=+1;

r+=1;

}
}
}
}
Mir=r;
M=r;


}
////////////////////// end make_RII



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




////////////////////// initial_X

void initial_X(int p)
{

int i,is;    
double prod;
double bas,ex;


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
      bas=Xi[i]/double(Nmax);
      ex=fabs(vl[r][l]);
      prod=prod*pow(bas,ex);
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
double bas,ex;
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
      bas=Xi[io]/double(Nmax);
      ex=fabs(vl[ro][lo]);
      prod=prod*pow(bas,ex);  
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
      bas=Xi[io]/double(Nmax);
      ex=fabs(vl[ro][lo]);
      prod=prod*pow(bas,ex);  
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




/////////////////// compute_P
void compute_P(int sample,int p)
{

int i1,i2,i3,i4,ir;  
double sum1,sum2,sum3;

  
for(int l=0;l<NR;l++){
      i1=Ri[l];
      for(int ll=0;ll<l+1;ll++){
        i2=Ri[ll];
        i3=max(i1,i2);
        i4=min(i1,i2);
	if(i3==i4){
	  Res[l][ll]+=avX[i3]/Nmax;	  
	}else{
	  Res[l][ll]+=(CI[i3][i4]-avX[i3]*avX[i4])/(Var[i3]*Var[i4]);
	}
      }
}
  
///

sum1=0;
sum2=0;
sum3=0;
for(int lr=0;lr<NR;lr++){
    ir=Ri[lr];
    sum1+=fabs(avX[ir]/Nmax-0.5);
    sum2+=fabs(avX[ir]/Nmax-0.5)*fabs(Yr[p][lr]/Nmax-0.5);
    sum3+=Var[ir]/(avX[ir]+1);
    Yr[p][lr]=avX[ir];
}
sum1=sum1/NR;
sum2=sum2/NR;
sum3=sum3/NR;
P1[p]+=sum1;
P2[p]+=sum1*sum1;
if(sample>0){
Q1[p]+=sum2;
Q2[p]+=sum2*sum2;
}
V1[p]+=sum3;
V2[p]+=sum3*sum3;


  
}
/////////////////// compute_P




/////////////////// compute_Q
void compute_Q()
{

double sum;


for(int p=0;p<Pat;p++){
for(int pp=0;pp<p;pp++){
  sum=0;
  for(int lr=0;lr<NR;lr++){
    sum+=(Yr[p][lr]/Nmax-0.5)*(Yr[pp][lr]/Nmax-0.5);
  }
  sum=sum/NR;
  QQ1[p][pp]+=sum;
  QQ2[p][pp]+=sum*sum;
}
}


}
/////////////////// compute_Q



/////////////////// report
void report(int sample)
{


int l1,l2;  
double m,p1,p2,q1,q2,v1,v2;
double sum1,sum2,sum3,sum4;

ofstream output1("avF.dat");
output1<<"#sample="<<sample+1<<endl;

sum1=0;
sum2=0;
sum3=0;
sum4=0;
for(int l=0;l<NR;l++){
  m=Res[l][l]/(sample+1);
  m=m/Pat;
  sum1+=m;
  sum2+=m*m;
  for(int ll=0;ll<NR;ll++){
    l1=max(l,ll);
    l2=min(l,ll);
    m=Res[l1][l2]/(sample+1);
    m=m/Pat;
    if(l<ll){
      sum3+=m;
      sum4+=m*m; 
    }
    output1<< m <<" ";
  }
  output1<<endl;
}
sum1=sum1/NR;
sum2=sum2/NR;
sum2=sqrt(sum2-sum1*sum1);
sum3=2*sum3/(NR*(NR-1));
sum4=2*sum4/(NR*(NR-1));
sum4=sqrt(sum4-sum3*sum3);
output1<<"#"<<" "<<sum1<<" "<<sum2<<" "<<sum3<<" "<<sum4<<endl;

///

ofstream output2("avP.dat");
output2<<"#sample="<<sample+1<<endl;


for(int p=0;p<Pat;p++){
  p1=P1[p]/(sample+1);
  p2=P2[p]/(sample+1);
  p2=sqrt((p2-p1*p1)/sample);
  q1=Q1[p]/sample;
  q2=Q2[p]/sample;
  q2=sqrt((q2-q1*q1)/(sample-1));
  v1=V1[p]/(sample+1);
  v2=V2[p]/(sample+1);
  v2=sqrt((v2-v1*v1)/sample);
  m=0;
  for(int ls=0;ls<NS;ls++)m+=(2*bit_map(p,ls)-1);
  m=m/NS;
  output2<< p <<" "<< p1 <<" "<< p2 <<" "<< q1 <<" "<<q2 <<" "<< v1 <<" "<<v2 <<" "<<m<<endl;
}

///

ofstream output3("avQ.dat");
output3<<"#sample="<<sample+1<<endl;

for(int p=0;p<Pat;p++){
  for(int pp=0;pp<p;pp++){
    p1=QQ1[p][pp]/(sample+1);
    p2=QQ2[p][pp]/(sample+1);
    p2=sqrt((p2-p1*p1)/sample);
    m=0;
    for(int ls=0;ls<NS;ls++){
      m+=(2*bit_map(p,ls)-1)*(2*bit_map(pp,ls)-1);
    }
    m=m/NS;
    output3<< p1 <<" "<< p2 <<" "<< m <<endl;
  }
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


