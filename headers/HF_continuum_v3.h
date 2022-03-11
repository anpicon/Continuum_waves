#include <cstdlib>
#include <cmath>
#include <complex>
#include <armadillo>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_coulomb.h>
#include <gsl/gsl_sf_legendre.h>

using namespace arma;

//function to calculate if integer is odd or even
bool isOdd( int integer )
{
    if ( integer % 2== 0 )
        return false;
    else
        return true;
}

//function to calculate integer position of mq
int mqi( int integer , int mmax )
{
    if (integer >=0) {
        return integer;
    }
    else
        return (2*mmax-integer);
}

//create wave in the continuum
void Create_Pelm(vec3x& Pelm,int& Lmax,int& mmax,int& iRmax)
{
    for (int L=0; L<=Lmax; L++) {
        int M1;
        if (L<=mmax) {M1=L;} else {M1=mmax;}
        for (int M=-M1; M<=M1; M++) {
            Pelm.push_back(vecE2x);
            for (int l=0; l<=Lmax; l++) {
                int M1;
                if (l<=mmax) {M1=l;} else {M1=mmax;}
                for (int m=-M1; m<=M1; m++) {
                    Pelm[SpH(L,M,mmax)].push_back(vecE1x);
                    for (int iR=0; iR<iRmax; iR++) Pelm[SpH(L,M,mmax)][SpH(l,m,mmax)].push_back(0.);
                }
            }
        }
    }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
//%%%%%%%%%%%%%%%% RADIAL INTEGRATION OF PARTIAL WAVES %%%%%%%%%%%%%%%%%//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
complex<double> ZZk (int k,double R, double rho1,complex<double> Zk,complex<double> Pab){
    complex<double> alpha;
    alpha= -(k/R)*Zk/rho1 + Pab/rho1;
    return alpha;
}
complex<double> YYk (int k,double R, double rho1,complex<double> Yk,complex<double> Zk){
    complex<double> alpha;
    alpha= ((k+1)/R)*Yk/rho1 - ((2*k+1)/R)*(Zk/rho1);
    return alpha;
}
//perform radial Yk integration given two Plm (P1,P2) radial functions
void do_Yk(cx_vec &Yk,vec1d &R,double dh,int iRmax,vec1d &rho_1,vec1d &rho_2,vec1d &rho_3,int k,int l1,int m1,int l2, int m2,int mmax,vec2x &P1,vec2x &P2,int mode){
    
    double value;
    bool itest=false;
    
    if (itest==true) {
        cout << "\nBeginning test:\n";
        for (int iR=0; iR<iRmax; iR++) {
            P1[SpH(l1,m1,mmax)][iR] = exp(-(R[iR]-3.)*(R[iR]-3.));
            P2[SpH(l2,m2,mmax)][iR] = exp(-(R[iR]-3.)*(R[iR]-3.));
        }
    }
    
    double Rth=4.; //we look for the distance that makes the integration steps negligible
    int iRlim=0;
    bool iP1P2_toolow=true;
    
    if (mode==0) {
        for (int iR=0; iR<iRmax; iR++) {
            if (R[iR]<Rth) {
                value = abs(conj(P1[SpH(l1,m1,mmax)][iR])*P2[SpH(l2,m2,mmax)][iR]);//*dh/rho_1(iR);
                if (value>1.e-13) {
                    iP1P2_toolow=false;
                }
            }
            if (R[iR]>Rth && iRlim==0) {
                if (iP1P2_toolow==true) return;
                value = abs(P1[SpH(l1,m1,mmax)][iR]);//*dh/rho_1(iR);
                if (value<1.e-14) {
                    iRlim=iR;
                }
            }
        }
    }
    else if (mode==1) {
        for (int iR=0; iR<iRmax; iR++) {
            if (R[iR]<Rth) {
                value = abs(conj(P1[SpH(l1,m1,mmax)][iR])*P2[SpH(l2,m2,mmax)][iR]);//*dh/rho_1(iR);
                if (value>1.e-12) {
                    iP1P2_toolow=false;
                }
            }
            if (R[iR]>Rth && iRlim==0) {
                if (iP1P2_toolow==true) return;
                value = abs(conj(P1[SpH(l1,m1,mmax)][iR])*P2[SpH(l2,m2,mmax)][iR]);//*dh/rho_1(iR);
                if (value<1.e-14) {
                    iRlim=iR;
                }
            }
        }
    }
    
    if (iRlim==0) iRlim=iRmax-2;
    if (isOdd(iRlim)==true) iRlim--;
    
    if (k>=0) { //Two first-order ordinary diff equations, with Adams-Moulton (3rd order) method
        vec1x Zk(iRlim+3,0.);
        for (int iR=1; iR<=1; iR++) {
            Zk[iR+1]=3.*ZZk(k,R[iR],rho_1[iR],Zk[iR],conj(P1[SpH(l1,m1,mmax)][iR])*P2[SpH(l2,m2,mmax)][iR]);
            Zk[iR+1]-=ZZk(k,R[iR-1],rho_1[iR-1],Zk[iR-1],conj(P1[SpH(l1,m1,mmax)][iR-1])*P2[SpH(l2,m2,mmax)][iR-1]);
            Zk[iR+1]*=dh/2.;
            Zk[iR+1]+=Zk[iR];
        }
        for (int iR=2; iR<=iRlim-1; iR++) {
            Zk[iR+1]=23.*ZZk(k,R[iR],rho_1[iR],Zk[iR],conj(P1[SpH(l1,m1,mmax)][iR])*P2[SpH(l2,m2,mmax)][iR]);
            Zk[iR+1]-=16.*ZZk(k,R[iR-1],rho_1[iR-1],Zk[iR-1],conj(P1[SpH(l1,m1,mmax)][iR-1])*P2[SpH(l2,m2,mmax)][iR-1]);
            Zk[iR+1]+=5.*ZZk(k,R[iR-2],rho_1[iR-2],Zk[iR-2],conj(P1[SpH(l1,m1,mmax)][iR-2])*P2[SpH(l2,m2,mmax)][iR-2]);
            Zk[iR+1]*=dh/12.;
            Zk[iR+1]+=Zk[iR];
        }
        Yk[iRlim]=Zk[iRlim];
        for (int iR=iRlim+1; iR<=iRlim+2; iR++) Zk[iR]=pow(R[iRlim]/R[iR],k)*Zk[iRlim];
        for (int iR=iRlim+1; iR<iRmax; iR++) Yk[iR]=pow(R[iRlim]/R[iR],k)*Yk[iRlim];
        for (int iR=iRlim; iR>=1; iR--) {
            Yk[iR-1]=23.*YYk(k,R[iR],rho_1[iR],Yk[iR],Zk[iR]);
            Yk[iR-1]-=16.*YYk(k,R[iR+1],rho_1[iR+1],Yk[iR+1],Zk[iR+1]);
            Yk[iR-1]+=5.*YYk(k,R[iR+2],rho_1[iR+2],Yk[iR+2],Zk[iR+2]);
            Yk[iR-1]*=-dh/12.;
            Yk[iR-1]+=Yk[iR];
        }
        
        //We have calculated Yk(rho)=r*y_k(r); we have to return y_k(r)
        for (int iR=0; iR<iRmax; iR++){
            Yk[iR]/=R[iR];
        }
    }
} //---End Do_Yk function



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
//%%%%%%%%%%%%%%%% CONTINUUM WAVE ASSUMING A HOLE IN THE GROUND STATE %%%%%%%%%%%%%%%%%//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
// epsilon -> energy continuum wave
// se -> spin continuum wave
// (kc,sc) -> orbital and spin of the hole state
void continuum_wave_1h(vec3x &Pelm,double epsilon,int se,cx_mat &Kmat,int iRj,vec1d &R,vec1d &RadialGrid_P,int Lmax,int mmax,vec2d &Atom_xyz,vec2i& Occj,vec3x &Pj) {
    bool iPrint=true;
    //for (int j=0; j<Occj.size(); j++) for (int s=0; s<2; s++) cout << "j " << j << " s " << s << " " << Occj[j][s] << endl;
    
    int socca=Occj.size();
    int natoms=Atom_xyz.size();
    int iRmax=R.size();
    double dh=RadialGrid_P[2];
    double Ar=RadialGrid_P[3];
    double Br=RadialGrid_P[4];
    vec1d rho_1(iRmax,0.); for (int iR=0; iR<iRmax; iR++) rho_1[iR]=dRx(Ar,Br,R[iR]);
    vec1d rho_2(iRmax,0.); for (int iR=0; iR<iRmax; iR++) rho_2[iR]=d2Rx(Ar,Br,R[iR]);
    vec1d rho_3(iRmax,0.); for (int iR=0; iR<iRmax; iR++) rho_3[iR]=d3Rx(Ar,Br,R[iR]);
    
    //Print for testing
    clock_t clock1,clock2,clock3,clock4,clock5; clock1=clock();
    if (iPrint==true) cout << "Radial grid: (" << R[0] << "," << R[iRmax-1] << ")" << endl;
    if (iPrint==true) cout << "Lmax: " << Lmax << " -> #angular basis: " << SpH(Lmax,Lmax,Lmax)+1 << endl;
    if (iPrint==true) cout << "Mmax: " << mmax << " -> #angular basis (constraint): " << SpH(Lmax,mmax,mmax)+1 << endl;
    
    //vec2i Occj(socca,vec1i(2,1)); for (int j=0; j<socca; j++) for (int s=0; s<2; s++) if(j==kc && s==sc) Occj[j][s]=0;
    
    
    //%%%%%%%% Calculation of Jlm functions
    cx_cube Jlm; Jlm.zeros(SpH(Lmax,mmax,mmax)+1,SpH(Lmax,mmax,mmax)+1,iRmax);
    
    bool iJlm_initial=true;
    if (iJlm_initial==true) {
        ifstream fp_initial;
        fp_initial.open("initial_Jlm.bin", ios::binary | ios::in);
        bool iJlm_read=false;
        //if (fp_initial) iJlm_read=true;
        fp_initial.close();
        if (iJlm_read==true) {
            fp_initial.open("initial_Jlm.bin", ios::binary | ios::in);
            for (int l=0; l<=Lmax; l++) {
                int M1;
                if (l<=mmax) {M1=l;} else {M1=mmax;}
                for (int m=-M1; m<=M1; m++) {
                    for (int ll=0; ll<=Lmax; ll++) {
                        int M2;
                        if (ll<=mmax) {M2=ll;} else {M2=mmax;}
                        for (int mm=-M2; mm<=M2; mm++) {
                            for (int iR=0; iR<iRmax; iR++) {
                                double value;
                                fp_initial.read((char*) &value, sizeof(double));
                                Jlm(SpH(l,m,mmax),SpH(ll,mm,mmax),iR)=value;
                                //real(Jlm(SpH(l,m,mmax),SpH(ll,mm,mmax),iR))=value;
                                fp_initial.read((char*) &value, sizeof(double));
                                Jlm(SpH(l,m,mmax),SpH(ll,mm,mmax),iR)+=c1*value;
                                //imag(Jlm(SpH(l,m,mmax),SpH(ll,mm,mmax),iR))=value;
                            }
                        }
                    }
                }
            }
            fp_initial.close();
            if (iPrint==true) printf("The potential Jlm has been read\n");
        }
        else{
            vec3x Jlm_temp(SpH(Lmax,mmax,mmax)+1,vec2x(SpH(Lmax,mmax,mmax)+1,vec1x(iRmax,0.)));
            vec2x Radial_part(SpH(2*Lmax,2*mmax,2*mmax)+1,vec1x(iRmax,0.)); //Radial_part[k/mq][iR]
            //Calculation of the radial part
            double th_W=1.e-14;
#pragma omp parallel for schedule(dynamic) default(shared)
            for (int k=0; k<=(2*Lmax); k++){
                cx_vec Yk; Yk.zeros(iRmax);//if this vector is out of the loop, it gives a numerical error -> omp.h
                for (int lj=0; lj<=Lmax; lj++) { //Loop lj and llj
                    for (int llj=0; llj<=Lmax; llj++) {
                        double Wj_lj_llj_k=gsl_sf_coupling_3j(2*lj,2*llj,2*k,0,0,0);
                        if (abs(Wj_lj_llj_k)>th_W)
                        {
                            double Sq_lj_llj=sqrt((2*lj+1)*(2*llj+1));
                            int M3=(lj<=mmax ? lj : mmax );
                            for (int mj=-M3; mj<=M3; mj++) { //Loop mj
                                int Mk=(k<=2*mmax ? k : 2*mmax );
                                for (int mq=-Mk; mq<=Mk; mq++) {
                                //int M4=(llj<=mmax ? llj : mmax );
                                //for (int mmj=-M4; mmj<=M4; mmj++) { //Loop mmj
                                    //int mq=mj-mmj; //mq is fixed with mmj and mj
                                    int mmj=mj-mq; //mmj is fixed with mq and mj
                                    double sign=(isOdd(mj+mq)==true ? -1. : 1. );
                                    double Wj_lj_llj_k_mj=gsl_sf_coupling_3j(2*lj,2*llj,2*k,-2*mj,2*mmj,2*mq);
                                    if (abs(Wj_lj_llj_k_mj)>th_W && abs(mmj)<=mmax)
                                    {
                                        for (int j=0; j<socca; j++) { //Loop j
                                            Yk.zeros();
                                            do_Yk(Yk,R,dh,iRmax,rho_1,rho_2,rho_3,k,lj,mj,llj,mmj,mmax,Pj[j],Pj[j],1);
#pragma omp simd
                                            for (int iR=0; iR<iRmax; iR++){
                                                Radial_part[SpH(k,mq,2*mmax)][iR]+=sign*Sq_lj_llj*Wj_lj_llj_k*Wj_lj_llj_k_mj*((Occj[j][0]+Occj[j][1])/(rho_1[iR]*rho_1[iR]))*Yk[iR];
                                            }
                                        }//-- End j loop
                                    }
                                } //end mq/mmj loop
                            } //end mj loop
                        }//fi
                    } //end llj loop
                }//end lj loop
            }//end k loop
            printf("Radial part of Jlm calculated...\n");
            
            //Sum over the angular part
#pragma omp parallel for schedule(dynamic) default(shared)
            for (int l=0; l<=Lmax; l++) {  //Loop l and ll
                for (int ll=0; ll<=Lmax; ll++) {
                    double Sq_l_ll=sqrt((2*l+1)*(2*ll+1));
                    for (int k=abs(l-ll); k<=l+ll; k++) {
                        double Wj_l_ll_k=gsl_sf_coupling_3j(2*ll,2*k,2*l,0,0,0);
                        if (abs(Wj_l_ll_k)>th_W)
                        {
                            int M1=(l<=mmax ? l : mmax );
                            for (int m=-M1; m<=M1; m++) {
                                double sign=(isOdd(m)==true ? -1. : 1. );
                                int M2=(ll<=mmax ? ll : mmax );
                                for (int mm=-M2; mm<=M2; mm++) { //Loop m and mm
                                    int mq=mm-m;//mq is fixed with m and mm
                                    double Wj_l_ll_k_m=gsl_sf_coupling_3j(2*ll,2*k,2*l,2*mm,-2*mq,-2*m);
                                    if (abs(Wj_l_ll_k_m)>th_W && abs(mq)<=2*mmax)
                                    {
#pragma omp simd
                                        for (int iR=0; iR<iRmax; iR++)
                                        {
                                            Jlm_temp[SpH(l,m,mmax)][SpH(ll,mm,mmax)][iR]+=sign*Sq_l_ll*Wj_l_ll_k*Wj_l_ll_k_m*Radial_part[SpH(k,mq,2*mmax)][iR];
                                        }
                                    
                                    }//fi
                                }//end mm loop
                            }//end m loop
                        }//fi
                    }//end k loop
                }//end ll loop
            }//end l loop
 
            printf("Jlm calculated...\n");
            
            ofstream fp_Jlm;
            fp_Jlm.open("initial_Jlm.bin", ios::binary | ios::out);
            for (int l=0; l<=Lmax; l++) {
                int M1;
                if (l<=mmax) {M1=l;} else {M1=mmax;}
                for (int m=-M1; m<=M1; m++) {
                    for (int ll=0; ll<=Lmax; ll++) {
                        int M2;
                        if (ll<=mmax) {M2=ll;} else {M2=mmax;}
                        for (int mm=-M2; mm<=M2; mm++) {
                            for (int iR=0; iR<iRmax; iR++) {
                                Jlm(SpH(l,m,mmax),SpH(ll,mm,mmax),iR)=Jlm_temp[SpH(l,m,mmax)][SpH(ll,mm,mmax)][iR];
                                double value;
                                value=real(Jlm(SpH(l,m,mmax),SpH(ll,mm,mmax),iR)); fp_Jlm.write((char*) &value, sizeof(double));
                                value=imag(Jlm(SpH(l,m,mmax),SpH(ll,mm,mmax),iR)); fp_Jlm.write((char*) &value, sizeof(double));
                            }
                        }
                    }
                }
            }
            fp_Jlm.close();
        }
    }
    if(Jlm.has_inf() || Jlm.has_nan()) printf("WARNING: Jlm has NaN/inf\n");
    
    
    
    //%%%%%%%% Calculation of Vlm functions
    cx_cube Vlm; Vlm.zeros(SpH(Lmax,mmax,mmax)+1,SpH(Lmax,mmax,mmax)+1,iRmax);
    bool iVlm_initial=true;
    if (iVlm_initial==true) {
        ifstream fp_initial;
        fp_initial.open("initial_Vlm.bin", ios::binary | ios::in);
        bool iVlm_read=false;
        //if (fp_initial) iVlm_read=true;
        fp_initial.close();
        if (iVlm_read==true) {
            fp_initial.open("initial_Vlm.bin", ios::binary | ios::in);
            for (int l=0; l<=Lmax; l++) {
                int M1;
                if (l<=mmax) {M1=l;} else {M1=mmax;}
                for (int m=-M1; m<=M1; m++) {
                    for (int ll=0; ll<=Lmax; ll++) {
                        int M2;
                        if (ll<=mmax) {M2=ll;} else {M2=mmax;}
                        for (int mm=-M2; mm<=M2; mm++) {
                            for (int iR=0; iR<iRmax; iR++) {
                                double value;
                                fp_initial.read((char*) &value, sizeof(double));
                                Vlm(SpH(l,m,mmax),SpH(ll,mm,mmax),iR)=value;
                                fp_initial.read((char*) &value, sizeof(double));
                                Vlm(SpH(l,m,mmax),SpH(ll,mm,mmax),iR)+=c1*value;
                            }
                        }
                    }
                }
            }
            fp_initial.close();
            if (iPrint==true) printf("The potential Vlm has been read\n");
        }
        else{
            vec3x Vlm_temp(SpH(Lmax,mmax,mmax)+1,vec2x(SpH(Lmax,mmax,mmax)+1,vec1x(iRmax,0.)));
            for (int i=0; i<natoms; i++) { //Note: we still need to define the # electron per MO!
                double Ri=sqrt(Atom_xyz[i][2]*Atom_xyz[i][2]+Atom_xyz[i][3]*Atom_xyz[i][3]+Atom_xyz[i][4]*Atom_xyz[i][4]);
                double Zi=Atom_xyz[i][1];
                double th_W=1.e-14;
                if (Ri<dh/100.) {
                    for (int iR=0; iR<iRmax; iR++) {
                        double r1=rho_1[iR];
                        for (int l=0; l<=Lmax; l++) {
                            int M1;
                            if (l<=mmax) {M1=l;} else {M1=mmax;}
                            for (int m=-M1; m<=M1; m++) {
                                Vlm_temp[SpH(l,m,mmax)][SpH(l,m,mmax)][iR]+=-1.*Zi/R[iR]/(r1*r1);
                            }
                        }
                    }
                }
                else {
                    double thetai=acos(Atom_xyz[i][4]/Ri);
                    double phii=atan2(Atom_xyz[i][3],Atom_xyz[i][2]);
                    for (int l=0; l<=Lmax; l++) {  //Loop l and ll
                        for (int ll=0; ll<=Lmax; ll++) {
                            double Sq_l_ll=sqrt((2*l+1)*(2*ll+1));
                            for (int k=abs(l-ll); k<=l+ll; k++) {
                                double Ak=sqrt(4.*pi/(2*k+1));
                                if (isOdd(l+ll+k)==false) { //condition that l+ll+k has to be even
                                    double Wj_l_ll_k=gsl_sf_coupling_3j(2*ll,2*k,2*l,0,0,0);
                                    if (abs(Wj_l_ll_k)>th_W) {
                                        int M1;
                                        if (l<=mmax) {M1=l;} else {M1=mmax;}
                                        for (int m=-M1; m<=M1; m++) {
                                            int M2;
                                            if (ll<=mmax) {M2=ll;} else {M2=mmax;}
                                            for (int mm=-M2; mm<=M2; mm++) { //Loop m and mm
                                                int mq=m-mm;//mq is fixed with m and mm
                                                double Wj_l_ll_k_m=gsl_sf_coupling_3j(2*ll,2*k,2*l,2*mm,2*mq,-2*m);
                                                if (abs(Wj_l_ll_k_m)>th_W) {
                                                    double sign=1.;
                                                    if (mq<0) {
                                                        if(isOdd(mq)==true) sign=-1.;
                                                        else sign=1.;
                                                    }
                                                    complex<double> Sph_k_mq=(sign*gsl_sf_legendre_sphPlm(k,abs(mq),cos(thetai)))*exp(-c1*(phii*mq));
                                                    if(isOdd(m)==true) sign=-1.; else sign=1.;
#pragma omp simd
                                                    for (int iR=0; iR<iRmax; iR++) {
                                                        double r1=rho_1[iR];
                                                        double temp=(R[iR]>Ri ? pow(Ri/R[iR],k)/R[iR]/(r1*r1) : pow(R[iR]/Ri,k)/Ri/(r1*r1) );
                                                        Vlm_temp[SpH(l,m,mmax)][SpH(ll,mm,mmax)][iR]+=-1.*Zi*sign*Ak*Sq_l_ll*Wj_l_ll_k*Wj_l_ll_k_m*Sph_k_mq*temp;
                                                        //if(R[iR]>Ri) Vlm_temp[SpH(l,m,mmax)][SpH(ll,mm,mmax)][iR]+=-1.*Zi*sign*Ak*Sq_l_ll*Wj_l_ll_k*Wj_l_ll_k_m*Sph_k_mq*pow(Ri/R[iR],k)/R[iR]/(r1*r1);
                                                        //else Vlm_temp[SpH(l,m,mmax)][SpH(ll,mm,mmax)][iR]+=-1.*Zi*sign*Ak*Sq_l_ll*Wj_l_ll_k*Wj_l_ll_k_m*Sph_k_mq*pow(R[iR]/Ri,k)/Ri/(r1*r1);
                                                    }
                                                } //--End if
                                            }
                                        }
                                    }
                                }//--End if
                            }//-- End k loop
                        }//-- End ll loop
                    }//-- End l loop
                }
            }//-- End i loop
            printf("Vlm calculated...\n");
            
            ofstream fp_Vlm;
            fp_Vlm.open("initial_Vlm.bin", ios::binary | ios::out);
            for (int l=0; l<=Lmax; l++) {
                int M1;
                if (l<=mmax) {M1=l;} else {M1=mmax;}
                for (int m=-M1; m<=M1; m++) {
                    for (int ll=0; ll<=Lmax; ll++) {
                        int M2;
                        if (ll<=mmax) {M2=ll;} else {M2=mmax;}
                        for (int mm=-M2; mm<=M2; mm++) {
                            for (int iR=0; iR<iRmax; iR++) {
                                Vlm(SpH(l,m,mmax),SpH(ll,mm,mmax),iR)=Vlm_temp[SpH(l,m,mmax)][SpH(ll,mm,mmax)][iR];
                                double value;
                                value=real(Vlm(SpH(l,m,mmax),SpH(ll,mm,mmax),iR)); fp_Vlm.write((char*) &value, sizeof(double));
                                value=imag(Vlm(SpH(l,m,mmax),SpH(ll,mm,mmax),iR)); fp_Vlm.write((char*) &value, sizeof(double));
                            }
                        }
                    }
                }
            }
            fp_Vlm.close();
            
        }
    } //--- end Vlm initial
    
    
    
    //%%%%%%%% Construction Flm matrix
    cx_cube Flm; Flm.zeros(SpH(Lmax,mmax,mmax)+1,SpH(Lmax,mmax,mmax)+1,iRmax);
    for (int l=0; l<=Lmax; l++) {
        int M1;
        if (l<=mmax) {M1=l;} else {M1=mmax;}
        for (int m=-M1; m<=M1; m++) {
            for (int iR=0; iR<iRmax; iR++) {
                double r1 = rho_1[iR];
                Flm(SpH(l,m,mmax),SpH(l,m,mmax),iR)=((l*(l+1)/R[iR]/R[iR])-2.*epsilon +0.5*rho_3[iR]/r1 - (3./4.)*rho_2[iR]*rho_2[iR]/(r1*r1))/(r1*r1);
            }
        }
    }
    Flm+=2.*Vlm+2.*Jlm;
    if(Flm.has_inf() || Flm.has_nan()) printf("WARNING: Flm has NaN/inf\n");
    
    clock2=clock();
    float diff=((float)clock2-(float)clock1);
    if(iPrint==true) cout<< "Potentials calculated in: " << diff/CLOCKS_PER_SEC/60. << " min\n";
    
    //%%%%%%%% Construction K-matrix
    //Rmat = - pi*Kmat
    //cx_mat Kmat;
    //Kmat.zeros(SpH(Lmax,mmax,mmax)+1,SpH(Lmax,mmax,mmax)+1); //K matrix
    
    //%%%%%%%%%%%%%%%% TERMS FOR TRIDIAGONAL MATRIX %%%%%%%%%%%%%%%%//
    //Defining Numerov coefficients  -> A(n+1) P(n+1) - B(n)P(n) + A(n-1)P(n-1) = fn + O(h^6)
    cx_mat Id; Id.eye(SpH(Lmax,mmax,mmax)+1,SpH(Lmax,mmax,mmax)+1);
    Vlm.zeros(); Jlm.zeros();
    for (int iR=0; iR<iRmax; iR++) Vlm.slice(iR)=(Id-dh*dh*(1./12.)*Flm.slice(iR));  //A coefficients - modify: symmatu
    for (int iR=0; iR<iRmax; iR++) Jlm.slice(iR)=(2.*Id+dh*dh*(10./12.)*Flm.slice(iR)); //B coefficients - modify: symmatu
    
    cx_cube Vtri; Vtri.zeros(SpH(Lmax,mmax,mmax)+1,SpH(Lmax,mmax,mmax)+1,iRmax); //NOTE: we are calculating Vtri again
    
    for (int iR=1; iR<iRmax-1; iR++) {
        Id=inv(Jlm.slice(iR)-(Vlm.slice(iR-1))*(Vtri.slice(iR-1)));
        Vtri.slice(iR)=Id*Vlm.slice(iR+1);
    }
    
    if(Vtri.has_inf() || Vtri.has_nan()) printf("WARNING: Vtri has NaN/inf\n");
    
    //%%%%%%%%%%%%%%%% LOOP WITH EXCHANGE POTENTIAL %%%%%%%%%%%%%%%%//
    cx_cube Klm; Klm.zeros(SpH(Lmax,mmax,mmax)+1,SpH(Lmax,mmax,mmax)+1,iRmax);
    bool iPelm_read=false; //we can control the initial Pelm
    if (iPrint==true) cout << "\nStarting loop to calculate exchange potential:" << endl;
    double thres_Kmat=1.e-7;
    int ixcMax=1;
    for (int ixc=0; ixc<ixcMax; ixc++) {
        double th_W=1.e-14;
        //Construction Klm matrix
        Klm.zeros();
        if (ixc!=0 || iPelm_read==true) {
#pragma omp parallel for schedule(dynamic) default(shared)
            for (int L=0; L<=Lmax; L++) {
                int M4;
                if (L<=mmax) {M4=L;} else {M4=mmax;}
                for (int M=-M4; M<=M4; M++) {
                    cx_cube Yg; Yg.zeros(SpH(2*Lmax,2*mmax,2*mmax)+1,socca,iRj);
                    cx_vec Yk; Yk.zeros(iRmax);
                    
                    //Generalized potential calculation Yg[j][kq][R]
                    for (int k=0; k<=2*Lmax; k++) {
                        int Mq;
                        if (k<=2*abs(mmax)) {Mq=k;} else {Mq=2*abs(mmax);}
                        for (int mq=-Mq; mq<=Mq; mq++) {
                            
                            for (int lj=0; lj<=Lmax; lj++) {
                                for (int ll=abs(k-lj); ll<=k+lj; ll++) {
                                    double Wj_lj_ll_k=gsl_sf_coupling_3j(2*lj,2*ll,2*k,0,0,0);
                                    if (isOdd(lj+ll+k)==false && ll<=Lmax && abs(Wj_lj_ll_k)>th_W) {
                                        double Sq_lj_ll=sqrt((2*lj+1)*(2*ll+1));
                                        int M3;
                                        if (lj<=mmax) {M3=lj;} else {M3=mmax;}
                                        for (int mj=-M3; mj<=M3; mj++) { //Loop mj
                                            int mm=mj-mq;
                                            double Wj_lj_ll_k_mm=gsl_sf_coupling_3j(2*lj,2*ll,2*k,-2*mj,2*mm,2*mq);
                                            if (abs(mm)<=ll && abs(mm)<=mmax && abs(Wj_lj_ll_k_mm)>th_W) {
                                                double sign=1.;
                                                if(isOdd(mj)==true) sign=-1.;
                                                
                                                for (int iR=0; iR<iRj; iR++) {
                                                    Pelm[SpH(L,M,mmax)][SpH(ll,mm,mmax)][iR]/=sqrt(rho_1[iR]);
                                                }
                                                
                                                for (int j=0; j<socca; j++) {
                                                    Yk.zeros();
                                                    do_Yk(Yk,R,dh,iRmax,rho_1,rho_2,rho_3,k,lj,mj,ll,mm,mmax,Pj[j],Pelm[SpH(L,M,mmax)],0);
                                                    
                                                    for (int iR=0; iR<iRj; iR++) {
                                                        Yg(SpH(k,mq,2*mmax),j,iR)+=sign*Sq_lj_ll*Wj_lj_ll_k*Wj_lj_ll_k_mm*Yk[iR];
                                                    }
                                                }//End j
                                                
                                                for (int iR=0; iR<iRj; iR++) {
                                                    Pelm[SpH(L,M,mmax)][SpH(ll,mm,mmax)][iR]*=sqrt(rho_1[iR]);
                                                }
                                                
                                            }
                                        }//End mj
                                    }
                                }//--End ll
                            }//--End lj
                            
                        }// End mq
                    }//--End k
                    //End calculation Yg
                    
                    
                    for (int l=0; l<=Lmax; l++) {
                        int M1;
                        if (l<=mmax) {M1=l;} else {M1=mmax;}
                        for (int m=-M1; m<=M1; m++) {
                            
                            for (int llj=0; llj<=Lmax; llj++) {
                                double Sq_l_llj=sqrt((2*l+1)*(2*llj+1));
                                int M2;
                                if (llj<=mmax) {M2=llj;} else {M2=mmax;}
                                for (int mmj=-M2; mmj<=M2; mmj++) {
                                    for (int k=abs(l-llj); k<=l+llj; k++) {
                                        if (isOdd(llj+l+k)==false) {
                                            double Wj_l_llj_k=gsl_sf_coupling_3j(2*llj,2*k,2*l,0,0,0);
                                            int mq=mmj-m;
                                            double Wj_l_llj_k_m=gsl_sf_coupling_3j(2*llj,2*k,2*l,2*mmj,-2*mq,-2*m);
                                            if (abs(Wj_l_llj_k_m)>th_W && abs(Wj_l_llj_k)>th_W) { //abs(mq)<=2*mmax &&
                                                double sign=1.;
                                                if(isOdd(mq+m)==true) sign=-1.;
                                                for (int j=0; j<socca; j++) {
                                                    for (int iR=0; iR<iRj; iR++) {
                                                        double r1=rho_1[iR];
                                                        Klm(SpH(L,M,mmax),SpH(l,m,mmax),iR)+=-2.*(Occj[j][se])*sign*Sq_l_llj*Wj_l_llj_k*Wj_l_llj_k_m*Yg(SpH(k,mq,2*mmax),j,iR)*Pj[j][SpH(llj,mmj,mmax)][iR]*sqrt(r1)/(r1*r1); //  Check negative sign
                                                    }
                                                }//End j
                                            }
                                        }
                                    }//--End k
                                }//End mmj
                            }//--End llj
                            
                        }//End m
                    }//--End l
                    
                }//--End M
            }//-- End L
            if (ixc==1) {
                clock3=clock();
                float diff=((float)clock3-(float)clock2);
                if(iPrint==true) cout<< "\n1st loop with exchange: " << diff/CLOCKS_PER_SEC/60. << " min\n";
            }
            if (ixc==2) {
                clock4=clock();
                float diff=((float)clock4-(float)clock3);
                if(iPrint==true) cout<< "\n2nd loop with exchange: " << diff/CLOCKS_PER_SEC/60. << " min\n";
            }
        }
        
        //%%%%%%%%%%%%%%%% TRIDIAGONAL MATRIX: FORWARD PROPAGATION %%%%%%%%%%%%%%%%//
        //Defining Numerov coefficients  -> A(n+1) P(n+1) - B(n)P(n) + A(n-1)P(n-1) = fn + O(h^6)
        cx_cube fn; fn.zeros(SpH(Lmax,mmax,mmax)+1,SpH(Lmax,mmax,mmax)+1,iRmax);
        for (int iR=1; iR<iRj-1; iR++) fn.slice(iR)=dh*dh*(1./12.)*((Klm.slice(iR+1))+10.*(Klm.slice(iR))+(Klm.slice(iR-1)));
        
        //Finding initial conditions
        cx_cube Utri; Utri.zeros(SpH(Lmax,mmax,mmax)+1,SpH(Lmax,mmax,mmax)+1,iRmax);
        
        //Forward sweep
        for (int iR=1; iR<iRmax-1; iR++) {
            Id=arma::inv(Jlm.slice(iR)-(Vlm.slice(iR-1))*(Vtri.slice(iR-1)));
            Utri.slice(iR)=((Utri.slice(iR-1))*strans(Vlm.slice(iR-1))-fn.slice(iR))*strans(Id);
        }
        
        if(Utri.has_inf() || Utri.has_nan()) printf("ixc %i WARNING: Utri has NaN/inf\n",ixc);
        
        //Boundaries at large distance & K-matrix  -> Vtri(iRmax-2) is required
        int iKmat=0;
        if (iKmat==0) { //1st method
            //Calculations of Coulomb waves in the boundaries
            double *fc_array1, *fc_array2, *gc_array1, *gc_array2;
            double *fc_exp1, *fc_exp2, *gc_exp1, *gc_exp2;
            fc_array1=new double [Lmax+1];
            gc_array1=new double [Lmax+1];
            fc_array2=new double [Lmax+1];
            gc_array2=new double [Lmax+1];
            fc_exp1=new double [Lmax+1];
            gc_exp1=new double [Lmax+1];
            fc_exp2=new double [Lmax+1];
            gc_exp2=new double [Lmax+1];
            double L_min=0;
            int kmax=Lmax;
            double ke=sqrt(2.*epsilon);
            double eta=-1./ke; //Z1*Z2/k Z stands for charge, k for linear momentum of electron in continuum
            gsl_sf_coulomb_wave_FG_array (L_min,kmax,eta,ke*R[iRmax-1],fc_array1,gc_array1,fc_exp1,gc_exp1);
            gsl_sf_coulomb_wave_FG_array (L_min,kmax,eta,ke*R[iRmax-2],fc_array2,gc_array2,fc_exp2,gc_exp2);
            
            cx_mat Gl1; Gl1.zeros(SpH(Lmax,mmax,mmax)+1,SpH(Lmax,mmax,mmax)+1);
            cx_mat Gl2; Gl2.zeros(SpH(Lmax,mmax,mmax)+1,SpH(Lmax,mmax,mmax)+1);
            cx_mat Fl1; Fl1.zeros(SpH(Lmax,mmax,mmax)+1,SpH(Lmax,mmax,mmax)+1);
            cx_mat Fl2; Fl2.zeros(SpH(Lmax,mmax,mmax)+1,SpH(Lmax,mmax,mmax)+1);
            for (int l=0; l<=Lmax; l++) {
                int M1;
                if (l<=mmax) {M1=l;} else {M1=mmax;}
                for (int m=-M1; m<=M1; m++) {
                    Gl1(SpH(l,m,mmax),SpH(l,m,mmax))=sqrt(rho_1[iRmax-1])*sqrt(2./pi/ke)*gc_array1[l];//*exp(gc_exp1[l]);
                    Gl2(SpH(l,m,mmax),SpH(l,m,mmax))=sqrt(rho_1[iRmax-2])*sqrt(2./pi/ke)*gc_array2[l];//*exp(gc_exp2[l]);
                    Fl1(SpH(l,m,mmax),SpH(l,m,mmax))=sqrt(rho_1[iRmax-1])*sqrt(2./pi/ke)*fc_array1[l];//*exp(fc_exp1[l]);
                    Fl2(SpH(l,m,mmax),SpH(l,m,mmax))=sqrt(rho_1[iRmax-2])*sqrt(2./pi/ke)*fc_array2[l];//*exp(fc_exp2[l]); //NOTE: check Gl1... for nonhomogeneous grid
                }
            }
            /*
             for (int l=0; l<=Lmax; l++) {
             cout << "gc..." << endl;
             cout << fc_array1[l] << " " << exp(fc_exp1[l]) << " " << fc_array2[l] << " " << exp(fc_exp2[l]) << " " << endl;
             cout << gc_array1[l] << " " << exp(gc_exp1[l]) << " " << gc_array2[l] << " " << exp(gc_exp2[l]) << " " << endl;
             }
             */
            delete[] fc_array1; delete[] fc_array2; delete[] gc_array1; delete[] gc_array2;
            delete[] fc_exp1; delete[] fc_exp2; delete[] gc_exp1; delete[] gc_exp2;
            
            //Calculating the K-matrix with the last matrices in the radial grid
            bool method_Kmat=true;
            if (method_Kmat==true) {
                cx_mat Kmat2; Kmat2.zeros(SpH(Lmax,mmax,mmax)+1,SpH(Lmax,mmax,mmax)+1); Kmat2=Kmat;
                
                if(Gl1.is_finite() && Gl2.is_finite() && Fl1.is_finite() && Fl2.is_finite()); else printf("ixc %i WARNING: Gl,Fl has NaN/inf\n",ixc);
                if((Vtri.slice(iRmax-2)).is_finite() && (Utri.slice(iRmax-2)).is_finite()); else printf("ixc %i WARNING: Vtri,Utri has NaN/inf\n",ixc);
                //Kmat=arma::inv(Gl2-(Vtri.slice(iRmax-2))*Gl1)*(strans(Utri.slice(iRmax-2))+(Vtri.slice(iRmax-2))*Fl1-Fl2);
                arma::solve(Kmat,(Gl2-(Vtri.slice(iRmax-2))*Gl1),(strans(Utri.slice(iRmax-2))+(Vtri.slice(iRmax-2))*Fl1-Fl2));
                if(Kmat.has_inf() || Kmat.has_nan()) printf("ixc %i WARNING: Kmat has NaN/inf\n",ixc);
                
                if (iPrint==true) cout << ixc << " "; //" Kmat_(i) - Kmat_(i-1)=";
                cout.precision(8); cout.setf( std::ios::fixed, std:: ios::floatfield );
                if (iPrint==true) cout << arma::norm(Kmat) << " " << arma::norm(Kmat2-Kmat) << endl;
                if(thres_Kmat>arma::norm(Kmat2-Kmat)) ixc=ixcMax;
                
                if (iPrint==true){
                    cx_mat tempKmat; tempKmat.zeros(SpH(Lmax,mmax,mmax)+1,SpH(Lmax,mmax,mmax)+1);
                    tempKmat=Kmat-Kmat.t();
                    if(abs(tempKmat.max())>1.e-13) {
                        printf("WARNING: Kmat no hermitian with max deviation %10.2e, max Kmat %10.2e\n",abs(tempKmat.max()),abs(Kmat.max()));
                        //cout << " WARNING: Kmat is not hermitian with max deviation " << tempKmat.max() << ", and max Kmat " << Kmat.max() << endl;
                    }
                }
            }
            
            //%%%%%%%%%%%%%%%% TRIDIAGONAL MATRIX: BACKWARDS PROPAGATION %%%%%%%%%%%%%%%%//
            //Calculation of Pelm in the last point of the radial grid
            bool iBackwards=true;
            if (iBackwards==true) {
                cx_mat tempKmat; tempKmat.zeros(SpH(Lmax,mmax,mmax)+1,SpH(Lmax,mmax,mmax)+1);
                tempKmat=symmatu(Kmat);
                
                cx_cube PelmMatrix; PelmMatrix.zeros(SpH(Lmax,mmax,mmax)+1,SpH(Lmax,mmax,mmax)+1,iRmax); //Plm for electron in continuum, Pelm.zeros() is this necessary?
                
                for (int L=0; L<=Lmax; L++) {
                    int M2;
                    if (L<=mmax) {M2=L;} else {M2=mmax;}
                    for (int M=-M2; M<=M2; M++) {
                        for (int l=0; l<=Lmax; l++) {
                            int M1;
                            if (l<=mmax) {M1=l;} else {M1=mmax;}
                            for (int m=-M1; m<=M1; m++) {
                                PelmMatrix(SpH(L,M,mmax),SpH(l,m,mmax),iRmax-1) = Fl1(SpH(L,M,mmax),SpH(l,m,mmax)) + tempKmat(SpH(l,m,mmax),SpH(L,M,mmax))*Gl1(SpH(l,m,mmax),SpH(l,m,mmax));
                            }
                        }
                    }
                }
                
                //Backward sweep
                for (int iR=iRmax-2; iR>=0; iR--) PelmMatrix.slice(iR)=Utri.slice(iR) + (PelmMatrix.slice(iR+1))*strans(Vtri.slice(iR));
                
                // NOTE: check pivot method
                
                if (ixc>=ixcMax-1) {
                    double testMax=0.;
                    cx_cube Vtest; Vtest.zeros(SpH(Lmax,mmax,mmax)+1,SpH(Lmax,mmax,mmax)+1,iRmax);
                    for (int iR=1; iR<iRmax-1; iR++){
                        Vtest.slice(iR)=(PelmMatrix.slice(iR-1))*strans(Vlm.slice(iR-1))-(PelmMatrix.slice(iR))*strans(Jlm.slice(iR))
                        +(PelmMatrix.slice(iR+1))*strans(Vlm.slice(iR+1))-fn.slice(iR);
                        if (testMax>arma::norm(Vtest.slice(iR),"inf")) {
                            testMax=arma::norm(Vtest.slice(iR),"inf");
                        }
                    }
                    
                    if(testMax<1.e-8) printf("Check of convergence: good\n");
                    else printf("Check of convergence: bad, with %e\n",testMax);
                }
                
                //Conversion from armadillo matrix to multi-vector
                for (int L=0; L<=Lmax; L++) {
                    int M2;
                    if (L<=mmax) {M2=L;} else {M2=mmax;}
                    for (int M=-M2; M<=M2; M++) {
                        for (int l=0; l<=Lmax; l++) {
                            int M1;
                            if (l<=mmax) {M1=l;} else {M1=mmax;}
                            for (int m=-M1; m<=M1; m++) {
                                for (int iR=0; iR<iRmax; iR++){
                                    Pelm[SpH(L,M,mmax)][SpH(l,m,mmax)][iR]=PelmMatrix(SpH(L,M,mmax),SpH(l,m,mmax),iR);
                                }
                            }
                        }
                    }
                }
                if(PelmMatrix.has_inf() || PelmMatrix.has_nan()) printf("ixc %i WARNING: PelmMatrix has NaN/inf\n",ixc);
            } //------ End backwards
        }
    } //--- END loop exchange potential calculation
    
    
    //Normalization of the solution with the incoming-wave condition
    bool iIncomeWave=true;
    if (iIncomeWave==true) {
        vec lambdaK; lambdaK.zeros(SpH(Lmax,mmax,mmax)+1);
        cx_mat Uk; Uk.zeros(SpH(Lmax,mmax,mmax)+1,SpH(Lmax,mmax,mmax)+1);
        Kmat=symmatu(Kmat); //we impose Kmat is exactly hermitian, make sure it is really hermitian
        arma::eig_sym(lambdaK,Uk,Kmat);
        
        cx_mat Id; Id.zeros(SpH(Lmax,mmax,mmax)+1,SpH(Lmax,mmax,mmax)+1);
        for (int l=0; l<=Lmax; l++) {
            int M1=(l<=mmax ? l : mmax);
            for (int m=-M1; m<=M1; m++) {
                double angle=atan(lambdaK(SpH(l,m,mmax)));
                Id(SpH(l,m,mmax),SpH(l,m,mmax))=exp(-c1*angle)*cos(angle);
            }//---end m
        }//---end l
        
        cx_cube PelmMatrix; PelmMatrix.zeros(SpH(Lmax,mmax,mmax)+1,SpH(Lmax,mmax,mmax)+1,iRmax);
        for (int L=0; L<=Lmax; L++) {
            int M2=(L<=mmax ? L : mmax);
            for (int M=-M2; M<=M2; M++) {
                for (int l=0; l<=Lmax; l++) {
                    int M1=(l<=mmax ? l : mmax);
                    for (int m=-M1; m<=M1; m++) {
                        for (int iR=0; iR<iRmax; iR++){
                            PelmMatrix(SpH(L,M,mmax),SpH(l,m,mmax),iR)=Pelm[SpH(L,M,mmax)][SpH(l,m,mmax)][iR];
                        }
                    }
                }
            }//---end M
        }//---end L
        for (int iR=0; iR<iRmax; iR++){
            PelmMatrix.slice(iR)=(strans(Uk*Id*trans(Uk))*PelmMatrix.slice(iR));
        }
        
        for (int L=0; L<=Lmax; L++) {
            int M2=(L<=mmax ? L : mmax);
            for (int M=-M2; M<=M2; M++) {
                for (int l=0; l<=Lmax; l++) {
                    int M1=(l<=mmax ? l : mmax);
                    for (int m=-M1; m<=M1; m++) {
                        for (int iR=0; iR<iRmax; iR++){
                            Pelm[SpH(L,M,mmax)][SpH(l,m,mmax)][iR]=PelmMatrix(SpH(L,M,mmax),SpH(l,m,mmax),iR);
                        }
                    }
                }
            }//---end M
        }//---end L
    }
    //cx_vec lambdaK; lambdaK.zeros(SpH(Lmax,Lmax)+1);
    //eig_gen(lambdaK,Uk,Kmat);
    
    //%%%%%%%% SAVING THE CONTINUUM WAVES
    ofstream fp_output; 
    for (int L=0; L<=Lmax; L++) {
        stringstream Lname;
        Lname.seekp(0,ios::beg); Lname << L;
        int M2=( L<=mmax ? L : mmax);
        //if (L<=mmax) {M2=L;} else {M2=mmax;}
        for (int M=-M2; M<=M2; M++) {
            stringstream Mname;
            Mname.seekp(0,ios::beg); Mname << M;
            string sname="Continuum_waves/wf_" + Lname.str() + "_" + Mname.str() + ".txt";
            fp_output.open(sname.c_str());
            fp_output << "R, ";
            for (int l=0; l<=Lmax; l++)
            {
                int M1=( l<=mmax ? l : mmax);
                for (int m=-M1; m<=M1; m++)
                {
                    fp_output << "l=" << l << " m=" << m << " , ";
                }//--- End m loop
            }//--- End l loop
            fp_output << endl;
            
            for (int iR=0; iR<iRmax; iR++) {
                //if (iR<iRj) {
                    fp_output << R[iR] << " ";
                    for (int l=0; l<=Lmax; l++)
                    {
                        int M1=( l<=mmax ? l : mmax);
                        for (int m=-M1; m<=M1; m++)
                        {
                            fp_output << Pelm[SpH(L,M,mmax)][SpH(l,m,mmax)][iR] << " ";
                        }//--- End m loop
                    }//--- End l loop
                    fp_output << endl;
                //}
            }//--- End R loop
            fp_output.close();
        } //--- End M loop
    } //--- End L loop
    
    
    clock5=clock();
    diff=((float)clock5-(float)clock2);
    if (iPrint==true) cout<< "Time for continuum diff. eqs.: " << diff/CLOCKS_PER_SEC/60. << " min" << endl;
    diff=((float)clock5-(float)clock1);
    if (iPrint==true) cout<< "Continuum part finished in   : " << diff/CLOCKS_PER_SEC/60. << " min" << endl;
    
} //---- END continuum_wave_1h
