#include <cstdlib>
#include <cmath>
#include <complex>
#include <armadillo>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_legendre.h>


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
//%%%%%%% ELECTRIC DIPOLE MATRIX ELEMENTS FOR TRANSITIONS BETWEEN TWO ORBITALS %%%%%%%%%%%%//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
// < 1| r | 2 >   where 2=(kc,sc) 1=(epsilon,se)
void Electric_dipoles(vec3x& Pelm,double epsilon,int se,int kc,int sc,int iRj,vec1d& R,vec1d& RadialGrid_P,int Lmax,int mmax,vec3x& Pj_ref,vec2x& dipoles) {
    
    if (sc!=se) {
        cout << "Photoelectron spin and core-hole spin are different" << endl;
    }
    
    int iRmax=R.size();
    double dh=RadialGrid_P[2];
    double Ar=RadialGrid_P[3];
    double Br=RadialGrid_P[4];
    vec1d rho_1(iRmax,0.); for (int iR=0; iR<iRmax; iR++) rho_1[iR]=dRx(Ar,Br,R[iR]);
    
    //%%%%%%%% CALCULATING RADIAL INTEGRALS
    vec3x Dip_R(SpH(Lmax,mmax,mmax)+1,vec2x(SpH(Lmax,mmax,mmax)+1,vec1x(SpH(Lmax,mmax,mmax)+1,complexd(0,0))));//array to store the single center expansion
    
    int iRlim; //it must be odd
    //if (isOdd(iRmax)==false) iRlim=iRmax-1;
    //else iRlim=iRmax-2;
    if (isOdd(iRj)==false) iRlim=iRj-1;
    else iRlim=iRj-2;
    
    for (int L=0; L<=Lmax; L++) {
        int M1=( L<=mmax ? L : mmax);
        //if (L<=mmax) {M1=L;} else {M1=mmax;}
        for (int M=-M1; M<=M1; M++) {
            
            for (int l=0; l<=Lmax; l++) {
                int M2=( l<=mmax ? l : mmax);
                //if (l<=mmax) {M2=l;} else {M2=mmax;}
                for (int m=-M2; m<=M2; m++) {
                    for (int lj=0; lj<=Lmax; lj++) {
                        int M3=( lj<=mmax ? lj : mmax);
                        //if (lj<=mmax) {M3=lj;} else {M3=mmax;}
                        for (int mj=-M3; mj<=M3; mj++) {
                            for (int iR=0; iR<=iRlim; iR++) { //Integration by Simpsons' rule, iRlim is odd, to have even points
                                double r1=rho_1[iR];
                                if (iR==0 || iR==iRlim) Dip_R[SpH(L,M,mmax)][SpH(l,m,mmax)][SpH(lj,mj,mmax)]+=R[iR]*conj(Pelm[SpH(L,M,mmax)][SpH(l,m,mmax)][iR]/sqrt(r1))*Pj_ref[kc][SpH(lj,mj,mmax)][iR]/r1;
                                else if (isOdd(iR)==true) Dip_R[SpH(L,M,mmax)][SpH(l,m,mmax)][SpH(lj,mj,mmax)]+=4.*R[iR]*conj(Pelm[SpH(L,M,mmax)][SpH(l,m,mmax)][iR]/sqrt(r1))*Pj_ref[kc][SpH(lj,mj,mmax)][iR]/r1;
                                else Dip_R[SpH(L,M,mmax)][SpH(l,m,mmax)][SpH(lj,mj,mmax)]+=2.*R[iR]*conj(Pelm[SpH(L,M,mmax)][SpH(l,m,mmax)][iR]/sqrt(r1))*Pj_ref[kc][SpH(lj,mj,mmax)][iR]/r1;
                            }
                            Dip_R[SpH(L,M,mmax)][SpH(l,m,mmax)][SpH(lj,mj,mmax)]*=dh/3.;
                        }
                    }
                }
            }
            //Once calculated the radial integrals...
    //%%%%%%%% CALCULATING DIPOLE TRANSITIONS
            for (int l=0; l<=Lmax; l++) {
                int M1=( l<=mmax ? l : mmax);;
                //if (l<=mmax) {M1=l;} else {M1=mmax;}
                for (int m=-M1; m<=M1; m++) {
                    for (int lj=abs(l-1); lj<=l+1; lj++) {
                        if (isOdd(l+lj+1)==false && lj<=Lmax) {
                            double Wj_l_lj=gsl_sf_coupling_3j(2*l,2*1,2*lj,0,0,0);
                            double Sq_l_lj=sqrt((2*l+1)*(2*lj+1));
                            int M2=( lj<=mmax ? lj : mmax);;
                            //if (lj<=mmax) {M2=lj;} else {M2=mmax;}
                            for (int mj=-M2; mj<=M2; mj++) {
                                int mq=mj-m;
                                if (std::abs(mq)<=1) {
                                    double Wj_l_lj_m=gsl_sf_coupling_3j(2*l,2*1,2*lj,-2*m,-2*mq,2*mj);
                                    if (abs(Wj_l_lj_m)>1.e-8) {
                                        int sign;
                                        if(isOdd(mq+m)==true) sign=-1;
                                        else sign=1;
                                        complex <double> temp1=sign*sqrt(4.*pi/3.)*Sq_l_lj*Wj_l_lj*Wj_l_lj_m*Dip_R[SpH(L,M,mmax)][SpH(l,m,mmax)][SpH(lj,mj,mmax)];
                                        sign=1;
                                        if (mq<0) {
                                            if(isOdd(mq)==true) sign=-1;
                                            else sign=1;
                                        }
                                        dipoles[0][SpH(L,M,mmax)]+=temp1*(sign*gsl_sf_legendre_sphPlm(1,abs(mq),0.)); //cos(pi/2.)=0.
                                        dipoles[1][SpH(L,M,mmax)]+=temp1*(sign*gsl_sf_legendre_sphPlm(1,abs(mq),0.))*exp(c1*(mq*pi)/2.); //cos(pi/2.)=0.
                                        dipoles[2][SpH(L,M,mmax)]+=temp1*(sign*gsl_sf_legendre_sphPlm(1,abs(mq),1.)); //cos(0)=1.
                                    }
                                }
                            }
                        }
                    }
                }
            }
            
        }//End M
    }//-- End L
    
    //%%%%%%%% SAVING DIPOLE TRANSITIONS
    ofstream fp_output;
    for (int i=0; i<3; i++) {
        string coord;
        if(i==0) coord="x";
        else if(i==1) coord="y";
        else coord="z";
        string sname="Photoionization/dipoles_" + coord + ".txt";
        fp_output.open(sname.c_str());
        fp_output << "Energy, ";
        for (int L=0; L<=Lmax; L++)
        {
            int M1=( L<=mmax ? L : mmax);
            for (int M=-M1; M<=M1; M++)
            {
                fp_output << "L=" << L << " M=" << M << " , ";
            }//--- End m loop
        }//--- End l loop
        fp_output << endl;
        
        fp_output << epsilon << " ";
        for (int L=0; L<=Lmax; L++)
        {
            int M1=( L<=mmax ? L : mmax);
            for (int M=-M1; M<=M1; M++)
            {
                fp_output << dipoles[i][SpH(L,M,mmax)] << " ";
                //if(i==0) fp_output << Dip_x(SpH(L,M,mmax)) << " ";
                //else if(i==1) fp_output << Dip_y(SpH(L,M,mmax)) << " ";
                //else fp_output << Dip_z(SpH(L,M,mmax)) << " ";
            }//--- End m loop
        }//--- End l loop
        fp_output << endl;
        fp_output.close();
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
}
