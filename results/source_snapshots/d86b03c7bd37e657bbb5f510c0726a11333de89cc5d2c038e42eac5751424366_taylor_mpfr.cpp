// Independent arbitrary-precision Taylor integration of state and sensitivities.
// Coefficients are normalized derivatives, computed from s*g'=alpha*g*s'.
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>
#include "mpfr_abi.hpp"

static long bits=256;
static int digits=70;
struct R {
    mpfr_t x;
    R() { mpfr_init2(x,bits);mpfr_set_si(x,0,MPFR_RNDN); }
    R(long a):R(){mpfr_set_si(x,a,MPFR_RNDN);}
    R(const std::string& a):R(){if(mpfr_set_str(x,a.c_str(),10,MPFR_RNDN))throw std::runtime_error("Invalid decimal");}
    R(const R& a):R(){mpfr_set(x,a.x,MPFR_RNDN);}
    R& operator=(const R& a){mpfr_set(x,a.x,MPFR_RNDN);return *this;}
    ~R(){mpfr_clear(x);}
    R& operator+=(const R& a){mpfr_add(x,x,a.x,MPFR_RNDN);return *this;}
    R& operator-=(const R& a){mpfr_sub(x,x,a.x,MPFR_RNDN);return *this;}
    std::string str()const {std::vector<char>b(digits+40);mpfr_snprintf(b.data(),b.size(),"%.*Re",digits,x);return b.data();}
    double dbl()const{return mpfr_get_d(x,MPFR_RNDN);}
};
R operator+(const R&a,const R&b){R z;mpfr_add(z.x,a.x,b.x,MPFR_RNDN);return z;}
R operator-(const R&a,const R&b){R z;mpfr_sub(z.x,a.x,b.x,MPFR_RNDN);return z;}
R operator*(const R&a,const R&b){R z;mpfr_mul(z.x,a.x,b.x,MPFR_RNDN);return z;}
R operator/(const R&a,const R&b){R z;mpfr_div(z.x,a.x,b.x,MPFR_RNDN);return z;}
R abs(const R&a){R z;mpfr_abs(z.x,a.x,MPFR_RNDN);return z;}
R sqrt(const R&a){R z;mpfr_sqrt(z.x,a.x,MPFR_RNDN);return z;}
bool operator<(const R&a,const R&b){return mpfr_cmp(a.x,b.x)<0;}
bool operator>(const R&a,const R&b){return b<a;}
using V=std::vector<R>;
using Series=std::vector<V>;
R conv(const V&a,const V&b,int k){R z;for(int l=0;l<=k;l++)z+=a[l]*b[k-l];return z;}

Series coefficients(const V&x,const V&mass,int n,int columns) {
    int size=12*(columns+1);
    Series y(size,V(n+1));for(int i=0;i<size;i++)y[i][0]=x[i];
    Series d(6,V(n)),s(3,V(n)),g(3,V(n)),w(3,V(n)),dd(9,V(n)),B(9,V(n));
    int ii[3]={0,0,1},jj[3]={1,2,2};
    R alpha=R(-3)/R(2),beta=R(-5)/R(2);
    for(int k=0;k<n;k++){
        for(int a=0;a<6;a++)y[a][k+1]=y[a+6][k]/R(k+1);
        for(int col=0;col<columns;col++)for(int a=0;a<6;a++)y[12+12*col+a][k+1]=y[18+12*col+a][k]/R(k+1);
        for(int pair=0;pair<3;pair++){
            int i=ii[pair],j=jj[pair];
            d[2*pair][k]=y[2*j][k]-y[2*i][k];d[2*pair+1][k]=y[2*j+1][k]-y[2*i+1][k];
            s[pair][k]=conv(d[2*pair],d[2*pair],k)+conv(d[2*pair+1],d[2*pair+1],k);
            if(k==0){g[pair][0]=R(1)/(s[pair][0]*sqrt(s[pair][0]));if(columns)w[pair][0]=g[pair][0]/s[pair][0];}
            else {
                R a,b;
                for(int l=1;l<=k;l++){
                    a+=((alpha+R(1))*R(l)-R(k))*s[pair][l]*g[pair][k-l];
                    if(columns)b+=((beta+R(1))*R(l)-R(k))*s[pair][l]*w[pair][k-l];
                }
                g[pair][k]=a/(R(k)*s[pair][0]);if(columns)w[pair][k]=b/(R(k)*s[pair][0]);
            }
            for(int a=0;a<2;a++){
                R force=conv(d[2*pair+a],g[pair],k)/R(k+1);
                y[6+2*i+a][k+1]+=mass[j]*force;y[6+2*j+a][k+1]-=mass[i]*force;
            }
            if(columns){
                dd[3*pair][k]=conv(d[2*pair],d[2*pair],k);
                dd[3*pair+1][k]=conv(d[2*pair],d[2*pair+1],k);
                dd[3*pair+2][k]=conv(d[2*pair+1],d[2*pair+1],k);
                B[3*pair][k]=g[pair][k]-R(3)*conv(dd[3*pair],w[pair],k);
                B[3*pair+1][k]=R(0)-R(3)*conv(dd[3*pair+1],w[pair],k);
                B[3*pair+2][k]=g[pair][k]-R(3)*conv(dd[3*pair+2],w[pair],k);
                for(int col=0;col<columns;col++)for(int a=0;a<2;a++){
                    R value;
                    for(int l=0;l<=k;l++)for(int b=0;b<2;b++){
                        int idx=(a==0&&b==0)?0:((a==1&&b==1)?2:1);
                        value+=B[3*pair+idx][k-l]*(y[12+12*col+2*j+b][l]-y[12+12*col+2*i+b][l]);
                    }
                    value=value/R(k+1);
                    y[18+12*col+2*i+a][k+1]+=mass[j]*value;y[18+12*col+2*j+a][k+1]-=mass[i]*value;
                }
            }
        }
    }
    return y;
}

int main(int argc,char**argv){
 try {
    if(argc!=11)throw std::runtime_error("Usage: taylor_mpfr digits order tolerance mu u v T columns output sample_stride");
    digits=std::stoi(argv[1]);bits=long(std::ceil((digits+8)*std::log2(10.)));
    int order=std::stoi(argv[2]),columns=std::stoi(argv[8]),stride=std::stoi(argv[10]);
    if(order<8 || (columns!=0&&columns!=2&&columns!=12) || stride<1)throw std::runtime_error("Invalid order/columns/stride");
    R tol(argv[3]),mu(argv[4]),u(argv[5]),v(argv[6]),T(argv[7]);
    V mass={R(1),mu,R(1)},x(12*(columns+1));
    x[0]=R(-1);x[4]=R(1);x[6]=u;x[7]=v;x[8]=R(-2)*u/mu;x[9]=R(-2)*v/mu;x[10]=u;x[11]=v;
    if(columns==12)for(int j=0;j<12;j++)x[12+12*j+j]=R(1);
    if(columns==2){for(int j=0;j<2;j++){x[18+12*j+j]=R(1);x[20+12*j+j]=R(-2)/mu;x[22+12*j+j]=R(1);}}
    V initial=x;R t(0);int steps=0,rejected=0;
    std::ofstream out(argv[9]);if(!out)throw std::runtime_error("Cannot open output");
    out<<"# MPFR "<<mpfr_get_version()<<" digits="<<digits<<" bits="<<bits<<" order="<<order<<" tol="<<argv[3]<<" columns="<<columns<<"\n";
    auto emit=[&](){out<<t.str();for(const R&a:x)out<<","<<a.str();out<<"\n";};emit();
    while(t<T){
        if(steps>200000)throw std::runtime_error("Step limit reached");
        Series a=coefficients(x,mass,order,columns);
        double proposal=.2;
        // Control the last TWO terms of every state/variational component.
        // Floating estimates choose h only; acceptance and evolution are MPFR.
        for(int i=0;i<int(x.size());i++)for(int k=order-1;k<=order;k++){
            R ratio=tol*(R(1)+abs(x[i]))/(abs(a[i][k])+R("1e-1000"));
            double logratio=std::log(ratio.dbl());
            double cand=.65*std::exp(logratio/k);
            if(std::isfinite(cand)&&cand>0)proposal=std::min(proposal,cand);
        }
        R h(std::to_string(proposal));
        // to_string gives six decimals; preserve small steps with scientific formatting.
        char buf[64];std::snprintf(buf,sizeof(buf),"%.17e",proposal);h=R(buf);
        if(T-t<h)h=T-t;
        for(;;){
            R hn(1);for(int k=0;k<order-1;k++)hn=hn*h;
            bool ok=true;
            for(int i=0;i<int(x.size());i++){
                R tail=abs(a[i][order-1]*hn)+abs(a[i][order]*hn*h);
                if(tail>tol*(R(1)+abs(x[i]))){ok=false;break;}
            }
            if(ok)break;
            h=h/R(2);rejected++;
        }
        if(h<R("1e-30"))throw std::runtime_error("Step underflow / collision candidate");
        for(int i=0;i<int(x.size());i++){
            R z=a[i][order];for(int k=order-1;k>=0;k--)z=z*h+a[i][k];x[i]=z;
        }
        t+=h;steps++;if(steps%stride==0||!(t<T))emit();
    }
    out<<"# success=1 steps="<<steps<<" rejected="<<rejected<<"\n";
    R closure;for(int i=0;i<12;i++)if(abs(x[i]-initial[i])>closure)closure=abs(x[i]-initial[i]);
    std::cout<<"success steps="<<steps<<" rejected="<<rejected<<" closure="<<closure.str()<<"\n";
    return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}
}
