/**
 * @author Victor Mingru Bai
 * @version 3.0 Oct21
 * Modeling Sugar Allocations in Plants using Radioisotope Tracer Data
 */

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "TMatrixD.h"
#include "TMath.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TCanvas.h"


//function readFile() that reads data from .txt files to array:
void readFile(Int_t len, Double_t myArray[len], string fileName){
    ifstream file(fileName.c_str());
    if(file.is_open()){
        for(Int_t i = 0; i < len; ++i){
            string tempStr;
            file >> tempStr;
            stringstream temp(tempStr);
            temp >> myArray[i];
        }
    }
    file.close();
}

//array contains:
bool contains(Double_t *array,Int_t num){
    bool check = false;
    int len = sizeof(array);
    for (int i=0; i<len; i++) {
        if (num==array[i]) {
            check=true;
        }
    }
    return check;
}

//chi-square:
Double_t chiSq(Int_t start, Int_t end, Int_t totaltime, Double_t xtotal[totaltime], Double_t mtotal[totaltime],Double_t noise[]){
    Double_t chichi=0;
    for (Int_t j=start;j<end;j++){
        if (contains(noise,j)==false){
            chichi=chichi+(mtotal[j]-xtotal[j])*(mtotal[j]-xtotal[j])/(mtotal[j]);
        }
    }
    return chichi;
}



//main simulation function:
void sim(Int_t totaltime,Double_t xloop[totaltime],Int_t dt,Double_t lambda,Double_t voldet,Double_t bw,Double_t mtotalder[totaltime],Int_t delay,Double_t jvnorm,Double_t xsectb[5],Double_t xsectf[4],Double_t xsectv[4],Double_t xsectvd[4],Double_t xtotal[totaltime]){
    
    //leafbins:
    Int_t leafbin=int(xsectb[1]/bw);
    TMatrixD xleafbin(leafbin,totaltime);
    
    //sugar in each section:
    Double_t xleaf[totaltime];
    Double_t xupperbin[totaltime];
    Double_t xupper[totaltime];
    Double_t xlowerbin[totaltime];
    Double_t xlower[totaltime];
    Double_t xrootbin[totaltime];
    Double_t xroot[totaltime];
    
    //sugar packets:
    Int_t totalpacket=leafbin*totaltime;
    TMatrixD xpacket(totalpacket,totaltime);
    TMatrixD xfront(totalpacket,totaltime);
    TMatrixD xback(totalpacket,totaltime);
    
    //sugar packet position:
    for (Int_t k=0; k<totalpacket; k++) {
        Int_t jmin=int(k/leafbin);
        Int_t i=k-int(k/leafbin)*leafbin;
        
        for (Int_t j=jmin;j<=jmin+delay && j<totaltime;j++){
            xfront(k,j)=(i+1)*bw;
            xback(k,j)=i*bw;
        }
        
        Int_t n1=0;
        for (Int_t j=jmin+delay; j<totaltime-1; j++) {
            if (xfront(k,j)<xsectb[n1+1]) {
                xfront(k,j+1)=xfront(k,j)+(xsectv[n1]+0.5*xsectvd[n1])*dt;
            } else {
                n1=n1+1;
                xfront(k,j+1)=xfront(k,j)+(xsectv[n1]+0.5*xsectvd[n1])*dt;
            }
        }
        
        Int_t n2=0;
        for (Int_t j=jmin+delay; j<totaltime-1; j++) {
            if (xback(k,j)<xsectb[n2+1]) {
                xback(k,j+1)=xback(k,j)+(xsectv[n2]-0.5*xsectvd[n2])*dt;
            } else {
                n2=n2+1;
                xback(k,j+1)=xback(k,j)+(xsectv[n2]-0.5*xsectvd[n2])*dt;
            }
        }
    }
    
    //inside leaf:
    for (Int_t i=0; i<leafbin; i++) {
        xleafbin(i,0)=mtotalder[0]*jvnorm*xloop[0]*1.0/(voldet*lambda*leafbin);
        for (Int_t j=0; j<totaltime-1; j++) {
            xleafbin(i,j+1)=xleafbin(i,j)*(1-lambda*dt)+mtotalder[j+1]*jvnorm*xloop[j+1]*1.0*(1-xsectf[0]*dt)/(voldet*lambda*leafbin);
        }
    }
    
    for (Int_t k=0; k<totalpacket; k++) {
        Int_t jmin=int(k/leafbin);
        xpacket(k,jmin)=mtotalder[jmin]*jvnorm*xloop[jmin]*1.0*xsectf[0]*dt/(voldet*lambda*leafbin);
        for (Int_t j=jmin; j<totaltime-1; j++) {
            if (xback(k,j)<xsectb[1]) {
                xpacket(k,j+1)=xpacket(k,j)*(1-lambda*dt);
            }
        }
    }
    
    for (Int_t j=0; j<totaltime; j++) {
        Double_t sum1=0;
        for (Int_t i=0; i<leafbin; i++) {
            sum1=sum1+xleafbin(i,j);
        }
        Double_t sum2=0;
        for (Int_t k=0; k<totalpacket; k++) {
            if (xback(k,j)<xsectb[1]) {
                sum2=sum2+xpacket(k,j)*(min(xfront(k,j),xsectb[1])-xback(k,j))*1.0/(max(xfront(k,j)-xback(k,j),double(2.0)));
            }
        }
        xleaf[j]=sum1+sum2;
    }
    
    //inside stem and root:
    for (Int_t k=0; k<totalpacket; k++) {
        for (Int_t j=0; j<totaltime-1; j++) {
            if (xfront(k,j)>xsectb[1]) {
                if (xfront(k,j)<xsectb[2]) {
                    xpacket(k,j+1)=xpacket(k,j)*(1-lambda*dt)-xpacket(k,j)*(1-lambda*dt)*xsectf[1]*dt*(xfront(k,j)-max(xback(k,j),xsectb[1]))*1.0/(xfront(k,j)-xback(k,j));
                } else if (xback(k,j)<xsectb[2]) {
                    xpacket(k,j+1)=xpacket(k,j)*(1-lambda*dt)-xpacket(k,j)*(1-lambda*dt)*xsectf[1]*dt*(xsectb[2]-xback(k,j))*1.0/(xfront(k,j)-xback(k,j))-xpacket(k,j)*xsectf[2]*dt*(xfront(k,j)-xsectb[2])*1.0/(xfront(k,j)-xback(k,j));
                } else if (xfront(k,j)<xsectb[3]) {
                    xpacket(k,j+1)=xpacket(k,j)*(1-lambda*dt)*(1-xsectf[2]*dt);
                } else if (xback(k,j)<xsectb[3]) {
                    xpacket(k,j+1)=xpacket(k,j)*(1-lambda*dt)-xpacket(k,j)*(1-lambda*dt)*xsectf[2]*dt*(xsectb[3]-xback(k,j))*1.0/(xfront(k,j)-xback(k,j))-xpacket(k,j)*xsectf[3]*dt*(xfront(k,j)-xsectb[3])*1.0/(xfront(k,j)-xback(k,j));
                } else {
                    xpacket(k,j+1)=xpacket(k,j)*(1-lambda*dt)*(1-xsectf[3]*dt);
                }
            }
        }
    }
    
    for (Int_t j=0; j<totaltime-1; j++) {
        Double_t xupperdif=0;
        Double_t xlowerdif=0;
        Double_t xrootdif=0;
        Double_t xupperdecayed=xupperbin[j]*(1-lambda*dt);
        Double_t xlowerdecayed=xlowerbin[j]*(1-lambda*dt);
        Double_t xrootdecayed=xrootbin[j]*(1-lambda*dt);
        
        for (Int_t k=0; k<totalpacket; k++) {
            if (xfront(k,j)>xsectb[1]) {
                if (xfront(k,j)<xsectb[2]) {
                    xupperdif=xupperdif+xpacket(k,j)*(1-lambda*dt)*xsectf[1]*dt*(xfront(k,j)-max(xback(k,j),xsectb[1]))*1.0/(xfront(k,j)-xback(k,j));
                } else if (xback(k,j)<xsectb[2]) {
                    xupperdif=xupperdif+xpacket(k,j)*(1-lambda*dt)*xsectf[1]*dt*(xsectb[2]-xback(k,j))*1.0/(xfront(k,j)-xback(k,j));
                    xlowerdif=xlowerdif+xpacket(k,j)*(1-lambda*dt)*xsectf[2]*dt*(xfront(k,j)-xsectb[2])*1.0/(xfront(k,j)-xback(k,j));
                } else if (xfront(k,j)<xsectb[3]) {
                    xlowerdif=xlowerdif+xpacket(k,j)*(1-lambda*dt)*xsectf[2]*dt;
                } else if (xback(k,j)<xsectb[3]) {
                    xlowerdif=xlowerdif+xpacket(k,j)*(1-lambda*dt)*xsectf[2]*dt*(xsectb[3]-xback(k,j))*1.0/(xfront(k,j)-xback(k,j));
                    xrootdif=xrootdif+xpacket(k,j)*(1-lambda*dt)*xsectf[3]*dt*(xfront(k,j)-xsectb[3])*1.0/(xfront(k,j)-xback(k,j));
                } else {
                    xrootdif=xrootdif+xpacket(k,j)*(1-lambda*dt)*xsectf[3]*dt;
                }
            }
        }
        xupperbin[j+1]=xupperdecayed+xupperdif;
        xlowerbin[j+1]=xlowerdecayed+xlowerdif;
        xrootbin[j+1]=xrootdecayed+xrootdif;
    }
    
    for (Int_t j=0; j<totaltime; j++) {
        Double_t xupperph=0;
        Double_t xlowerph=0;
        Double_t xrootph=0;
        for (Int_t k=0; k<totalpacket; k++) {
            if (xfront(k,j)>xsectb[1]) {
                if (xfront(k,j)<xsectb[2]) {
                    xupperph=xupperph+xpacket(k,j)*(xfront(k,j)-max(xback(k,j),xsectb[1]))*1.0/(xfront(k,j)-xback(k,j));
                } else if (xback(k,j)<xsectb[2]) {
                    xupperph=xupperph+xpacket(k,j)*(xsectb[2]-xback(k,j))*1.0/(xfront(k,j)-xback(k,j));
                    xlowerph=xlowerph+xpacket(k,j)*(xfront(k,j)-xsectb[2])*1.0/(xfront(k,j)-xback(k,j));
                } else if (xfront(k,j)<xsectb[3]) {
                    xlowerph=xlowerph+xpacket(k,j);
                } else if (xback(k,j)<xsectb[3]) {
                    xlowerph=xlowerph+xpacket(k,j)*(xsectb[3]-xback(k,j))*1.0/(xfront(k,j)-xback(k,j));
                    xrootph=xrootph+xpacket(k,j)*(xfront(k,j)-xsectb[3])*1.0/(xfront(k,j)-xback(k,j));
                } else {
                    xrootph=xrootph+xpacket(k,j);
                }
            }
        }
        xupper[j]=xupperbin[j]+xupperph;
        xlower[j]=xlowerbin[j]+xlowerph;
        xroot[j]=xrootbin[j]+xrootph;
    }
    
    //producing output data array:
    xtotal[0]=0.0;
    for (Int_t j=0; j<totaltime-1; j++) {
        xtotal[j+1]=xupper[j]+xlower[j]+xroot[j];
    }    
}


void exe(){
    //importing data:
    Int_t totaltime=250;
    Double_t xloop[totaltime];
    Double_t mleaf[totaltime];
    Double_t mupper[totaltime];
    Double_t mlower[totaltime];
    Double_t mroot[totaltime];
    readFile(totaltime,xloop,"/Users/mb298/Root/1499/1499_line.txt");
    readFile(totaltime,mleaf,"/Users/mb298/Root/1499/1499_leaf.txt");
    readFile(totaltime,mupper,"/Users/mb298/Root/1499/1499_uppershoot.txt");
    readFile(totaltime,mlower,"/Users/mb298/Root/1499/1499_lowershoot.txt");
    readFile(totaltime,mroot,"/Users/mb298/Root/1499/1499_root.txt");
    
    //experiment setup parameters:
    Int_t dt=1;
    Double_t lambda=0.034;
    Double_t voldet=20;
    Double_t bw=2;
    
    //jv ratio:
    Double_t jmtotal[totaltime];
    for (Int_t j=0; j<totaltime; j++) {
        jmtotal[j]=mleaf[j]+mupper[j]+mlower[j]+mroot[j];
    }
    Double_t mtotalder[totaltime];
    mtotalder[0]=jmtotal[1]-jmtotal[0];
    for (Int_t j=0; j<totaltime-1; j++) {
        mtotalder[j]=(jmtotal[j+1]-jmtotal[j]+lambda*jmtotal[j])/xloop[j];
    }
    
    //normalization:
    Double_t jvnorm=0.6801;
    
    //imager noise data points:
    Double_t noise[]={48,59,71,83,94,106,118,129,141,152};
    
    //section parameters:
    Int_t delay=0;
    Double_t xsectb[5]={0,40,115,190,10000};
    Double_t xsectf[4]={0.00725,0.007,0.055,0.05};
    Double_t xsectv[4]={3.96,0,3.5,2.01};
    Double_t xsectvd[4]={0,0,0,0};
    
    //data arrays:
    Double_t timeindex[totaltime];
    Double_t xtotal[totaltime];
    Double_t mtotal[totaltime];
    for (Int_t j=0;j<totaltime;j++){
        timeindex[j]=j;
        mtotal[j]=mupper[j]+mlower[j]+mroot[j];
    }
    
    //search range:
    Int_t sstart=50;
    Int_t send=110;
    
    //performing imulation:
    Double_t minChi=10000000;
    Double_t fmin;
    Double_t vmin;
    Double_t dmin;
    for (int q1=0;q1<10;q1++){
        for (int q2=0; q2<10; q2++) {
            for (int q3=0; q3<1; q3++) {
                xsectf[0]=0.5+0.01*q1;
                xsectv[0]=0.4+0.01*q2;
                delay=5+1*q3;
                sim(totaltime,xloop,dt,lambda,voldet,bw,mtotalder,delay,jvnorm,xsectb,xsectf,xsectv,xsectvd,xtotal);
                Double_t chichi=chiSq(sstart, send, totaltime, xtotal, mtotal, noise);
                if (chichi<minChi) {
                    minChi=chichi;
                    fmin=xsectf[0];
                    vmin=xsectv[0];
                    dmin=delay;
                }
            }
        }
    }
        
    xsectf[0]=fmin;
    xsectv[0]=vmin;
    delay=dmin;
    
    cout<<xsectf[0]<<endl;
    cout<<xsectv[0]<<endl;
    cout<<delay<<endl;
    cout<<minChi<<endl;

    sim(totaltime,xloop,dt,lambda,voldet,bw,mtotalder,delay,jvnorm,xsectb,xsectf,xsectv,xsectvd,xtotal);
    
    //draw graphs:
    TCanvas *sim=new TCanvas("sim","simulation_1499");
    TGraph xtotalgraph(totaltime,timeindex,xtotal);
    TGraph mtotalgraph(totaltime,timeindex,mtotal);
    mtotalgraph.SetLineColor(kRed);
    mtotalgraph.DrawClone("CA");
    xtotalgraph.DrawClone("SAME");
    sim->Print("simulation_1499");
    
    

}