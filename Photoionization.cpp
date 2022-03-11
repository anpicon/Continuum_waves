/* Created by APA 05/07/2021
 Program to calculate a continuum wave based on k-matrix approach and static-exchange approximation
 */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <complex>
#include <ctime>
#include <vector>
#include <map>
#include <utility>
#include <armadillo>

#include "omp.h"
#include "../headers/Constants.h"
#include "../headers/typedef.h"
#include "../headers/Loading_wavefunction.h"
#include "../headers/HF_continuum.h"
#include "../headers/Observables.h"

using namespace std;

int main (int argc, char* argv[])
{
    //%%%%%%%%%%%%%%%%%%% READING INPUT
    #include "Source_Main/input.cpp"

    //%%%%%%%%%%%%%%%%%%% PRINTING INPUT
    #include "Source_Main/print_input.cpp"
    
    //%%%%%%%% DEFINING RADIAL GRID
    vec1d RadialGrid_P(5,0.); RadialGrid_P[0]=Rmax; RadialGrid_P[1]=Rj; RadialGrid_P[2]=dh; RadialGrid_P[3]=Ar; RadialGrid_P[4]=Br;
    double Rmin=0.001;
    double rhoMax=Rx(Ar,Br,0,Rmax,0);
    double rhojmax=Rx(Ar,Br,0,Rj,0);
    double rhoMin=Rx(Ar,Br,0,Rmin,0);
    int iRmax=int((rhoMax-rhoMin)/dh);
    int iRj=int((rhojmax-rhoMin)/dh);
    printf("#radial points iRmax: %d \n",iRmax);
    printf("#radial points iRj: %d \n\n",iRj);
    
    //%%%%%%%% ARRAYS AND VARIABLES
    int TMO=0; Sy.TMO(TMO); //count total molecular orbitals
    vec2d Geometry; //store geometry molecule
    vec3x Plm(TMO,vec2x(SpH(Lmax,mmax,mmax)+1,vec1x(iRmax,complexd(0,0))));//array to store the single center expansion
    vec2i Occj(TMO,vec1i(2,1)); //array to store the occupation in a single determinant
    vec1d rho(iRmax,0.); for (int iR=0; iR<iRmax; iR++) rho[iR]=rhoMin+iR*dh;
    vec1d R(iRmax,0.); mapping_rho_to_R(Ar,Br,Rmin,Rmax,R,rho,iRmax); //for (int iR=0; iR<iRmax; iR++) R[iR]=dh/2.+iR*dh; #Array to store radial grid
    
    vec1d EnergyMO(TMO,0.);//Save energies of the HF MOs

    //%%%%%%%% READING CI VECTOR OF BOUND STATES
    Load_CI(Occj,Sy.NMO,pathocc);
    
    //%%%%%%%% SINGLE CENTER EXPANSION FOR BOUND STATE MOs
    #include "Source_Main/load_mo.cpp"

    //%%%%%%%% READING ENERGIES OF HF MOs
    Load_BE_HF(Occj,Sy.NMO,EnergyMO,pathBE);
    
    //%%%%%%%% CALCULATING PHOTOIONIZATION CROSS SECTIONS
    system("mkdir Continuum_waves");
    system("mkdir Photoionization");
    ofstream fp_output;
    sname="Photoionization/cross_sections.txt";
    fp_output.open(sname.c_str());
    fp_output << " MO       Photoelectron(eV)       Dx(au)          Dy(au)          Dz(au)          D**2(Mb)" << endl;
    
    printf("\nAtom   Z      x(a.u.)          y(a.u.)          z(a.u.)\n");
    for (int i=0; i<Geometry.size(); i++) {
        printf("%i      %1.f    %2.4f            %2.4f            %2.4f\n",i+1, Geometry[i][1], Geometry[i][2], Geometry[i][3], Geometry[i][4]);
    }
    printf("\n\n");
    
    
    for(int ic=0; ic<TMO; ic++){
        vec1i ij; ij=Sy.ij(ic);
        int hspin=0;
        if(Occj[ic][0]==1 && Occj[ic][1]==1){
            //declaring continuum wave Pelm
            vec3x Pelm; Create_Pelm(Pelm,Lmax,mmax,iRmax);

            //declaring K-matrix
            cx_mat Kmat; Kmat.zeros(SpH(Lmax,mmax,mmax)+1,SpH(Lmax,mmax,mmax)+1); //K matrix, Rmat = - pi*Kmat
            
            //adding the hole in the occupation array
            Occj[ic][hspin]=0;
            
            // epsilon -> energy continuum wave
            // se -> spin continuum wave
            // (kc,sc) -> orbital and spin of the hole state
            int se=hspin; //Ionization of core-hole -> hspin (spin of core hole) = se (spin of continuum electron)
            double epsilon=omega+EnergyMO[ic];
            if (epsilon>0.) {
                printf("\nMolecular orbital %d.%d: %3.3f\n",ij[1]+1,ij[0]+1,epsilon);
                continuum_wave_1h(Pelm,epsilon,se,Kmat,iRj,R,RadialGrid_P,Lmax,mmax,Geometry,Occj,Plm);
                //Calculating photoionization cross sections
                vec2x dipoles(3,vec1x(SpH(Lmax,mmax,mmax)+1,complexd(0,0)));
                Electric_dipoles(Pelm,epsilon,se,ic,hspin,iRj,R,RadialGrid_P,Lmax,mmax,Plm,dipoles);
                double TotalXS_x=0.;
                double TotalXS_y=0.;
                double TotalXS_z=0.;
                for (int L=0; L<=Lmax; L++)
                {
                    int M1=( L<=mmax ? L : mmax);
                    for (int M=-M1; M<=M1; M++)
                    {
                        TotalXS_x+=norm(dipoles[0][SpH(L,M,mmax)]);
                        TotalXS_y+=norm(dipoles[1][SpH(L,M,mmax)]);
                        TotalXS_z+=norm(dipoles[2][SpH(L,M,mmax)]);
                    }
                }
                double temp=4.*pi*pi*omega/(3.*c);
                TotalXS_x*=temp;
                TotalXS_y*=temp;
                TotalXS_z*=temp;
                fp_output << ij[1] << "." << ij[0] << "       ";
                fp_output << epsilon*energy_au_eV << "          " << TotalXS_x << "          " << TotalXS_y << "          " << TotalXS_z << "          " << 2.*(TotalXS_x+TotalXS_y+TotalXS_z)*XSections_au_Mb << endl;
            }
            //adding the electron in the occupation array
            Occj[ic][hspin]=1;
        }//fi
    } // end ic -- loop over all MOs
    fp_output.close();

    
    
    
    
    
    
    
return 0;
}
