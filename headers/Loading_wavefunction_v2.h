/*// old version with likely library
#include "LikelyTricubic3_3.h"
#include "LikelyTricubic3_2.h"
#include "LikelyTricubic3.h"
 */
#include "Interpolation.h"
#include <gsl/gsl_sf_legendre.h>

//Class to count different symmetries
class Symmetry
{
public:
    Symmetry(string&);
    void TMO(int&);
    vec1i NMO;
    vec1i ij(int&);
    
private:
    int irrerep;
};

Symmetry::Symmetry(string& PointGroup)
{
    if (PointGroup=="D2h") {
        irrerep=8;
    }
    else if (PointGroup=="C2v" || PointGroup=="C2h" || PointGroup=="D2") {
        irrerep=4;
    }
    else if (PointGroup=="Cs" || PointGroup=="C2" || PointGroup=="Ci") {
        irrerep=2;
    }
    else {
        irrerep=1;
    }
    
    NMO.resize(irrerep,0);
}

void Symmetry::TMO(int& TMO)
{
    TMO=0;
    for (std::vector<int>::size_type i=0; i<NMO.size(); i++) {
        for (int j=0; j<NMO[i]; j++) {
            TMO++;
        }
    }
}

vec1i Symmetry::ij(int& k)
{
    int kk=0;
    vec1i ij(2,0);
    for (std::vector<int>::size_type i=0; i<NMO.size(); i++) {
        for (int j=0; j<NMO[i]; j++) {
            if(kk==k) {
                ij[0]=i;
                ij[1]=j;
            }
            kk++;
        }
    }
    return ij;
}



//Notation for saving spherical harmonic basis
int SpH(int l, int m, int mmax){
    int A=0;
    if (l<=mmax) {
        for (int j=0; j<l; j++) A+=(2*j+1);
        A--;
        if (m==0) A+=l+1;
        else A+=l+1+m;
    }
    else {
        for (int j=0; j<=mmax; j++) A+=(2*j+1);
        for (int j=mmax+1; j<l; j++) A+=(2*mmax+1);
        A--;
        if (m==0) A+=mmax+1;
        else A+=mmax+1+m;
    }
    return A;
}
//Function to load CI vectors and amplitudes
void Load_CI(vec2i& Occj,vec1i& NMO,string pathocc){
    ifstream fp_input;
    fp_input.open(pathocc.c_str());
    int TMO=0;
    for (std::vector<int>::size_type i=0; i<NMO.size(); i++) {
        for (int j=0; j<NMO[i]; j++) {
            string value;
            fp_input >> value;
            if(value=="2"){
                Occj[TMO][0]=1; Occj[TMO][1]=1;
            }
            else if(value=="a"){
                Occj[TMO][0]=1; Occj[TMO][1]=0;
            }
            else if(value=="b"){
                Occj[TMO][0]=0; Occj[TMO][1]=1;
            }
            else{
                Occj[TMO][0]=0; Occj[TMO][1]=0;
            }
            TMO++;
        }
    }
    fp_input.close();
}



//Function to load binding energies of a HF calculation
void Load_BE_HF(vec2i& Occj,vec1i& NMO,vec1d& EnergyMO,string pathBE){
    ifstream fp_input;
    fp_input.open(pathBE.c_str());
    int TMO=0;
    for (std::vector<int>::size_type i=0; i<NMO.size(); i++) {
        for (int j=0; j<NMO[i]; j++) {
            if(Occj[TMO][0]==1 && Occj[TMO][1]==1){
                double value;
                fp_input >> value >> value;
                EnergyMO[TMO]=value;
            }//fi
            TMO++;
        }
    }
    fp_input.close();
}


//Radial grid
double Rx(double Ar,double Br,double rho,double R,int mode){
    double fx;
    if (mode==0) fx=Ar*R+Br*log(R);
    else if (mode==1) fx=Ar*R+Br*log(R)-rho;
    return fx;
}
double dRx(double Ar,double Br,double R){
    double dfx=Ar+Br/R;
    return dfx;
}
double d2Rx(double Ar,double Br,double R){
    double dfx=-Br/(R*R);
    return dfx;
}
double d3Rx(double Ar,double Br,double R){
    double dfx=+2.*Br/(R*R*R);
    return dfx;
}

//Given the lower and upper limits of integration x1 and x2,
//this routine returns arrays x[0..n-1] and w[0..n-1] of length n,
//containing the abscissas and weights of the Gauss-Legendre n-point quadrature formula.
//Numerical Recipies 3rd edition
void gauleg(double x1,double x2, vec1d& x, vec1d& w)
{
    double EPS=1.0e-14; //EPS is the relative precision
    double z1,z,xm,xl,pp,p3,p2,p1;
    int n=x.size();
    int m=(n+1)/2; //The roots are symmetric in the interval, so we only have to find half of them.
    xm=0.5*(x2+x1);
    xl=0.5*(x2-x1);
    for (int i=0;i<m;i++) { //Loop over the desired roots.
        z=cos(3.141592654*(i+0.75)/(n+0.5));
        //Starting with this approximation to the ith root, we enter the main loop of refinement by Newton’s method.
        do {
            p1=1.0;
            p2=0.0;
            for (int j=0;j<n;j++) {//Loop up the recurrence relation to get the Legendre polynomial evaluated at z.
                p3=p2;
                p2=p1;
                p1=((2.0*j+1.0)*z*p2-j*p3)/(j+1);
            }
            //p1 is now the desired Legendre polynomial.
            //We next compute pp, its derivative, by a standard relation involving also p2, the polynomial of one lower order.
            pp=n*(z*p1-p2)/(z*z-1.0);
            z1=z;
            z=z1-p1/pp;
        } while (abs(z-z1) > EPS); //Newton's method
        x[i]=xm-xl*z;
        x[n-1-i]=xm+xl*z; //Scale the root to the desired interval, and put in its symmetric counterpart.
        w[i]=2.0*xl/((1.0-z*z)*pp*pp);
        w[n-1-i]=w[i]; //Compute the weight and its symmetric counterpart.
    }
}

double rtsafe(double Ar,double Br,double rho,double x1, double x2, double xacc) {
    //Using a combination of Newton-Raphson and bisection, return the root of a function bracketed between x1 and x2. The root will be refined until its accuracy is known within  ̇xacc. funcd is a user-supplied struct that returns the function value as a functor and the first derivative of the function at the point x as the function df (see text).
    int MAXIT=100; //Maximum allowed number of iterations.
    double xh,xl;
    double fl=Rx(Ar,Br,rho,x1,1);//funcd(x1);
    double fh=Rx(Ar,Br,rho,x2,1);//funcd(x2);
    if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
        throw("Root must be bracketed in rtsafe");
    if (fl == 0.0) return x1;
    if (fh == 0.0) return x2;
    if (fl < 0.0) { //Orient the search so that f .xl/ < 0.
        xl=x1;
        xh=x2;
    } else {
        xh=x1;
        xl=x2; }
    double rts=0.5*(x1+x2);  //Initialize the guess for root, the “stepsize before last,” and the last step.
    double dxold=abs(x2-x1);
    double dx=dxold;
    double f=Rx(Ar,Br,rho,rts,1);//funcd(rts);
    double df=dRx(Ar,Br,rts); //funcd.df(rts);
    for (int j=0;j<MAXIT;j++) { //Loop over allowed iterations.
        if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0) //Bisect if Newton out of range, or not decreasing fast enough.
            || (abs(2.0*f) > abs(dxold*df))) {
            dxold=dx;
            dx=0.5*(xh-xl);
            rts=xl+dx;
            if (xl == rts) return rts; //Change in root is negligible. Newton step acceptable. Take it.
        }
        else {
            dxold=dx;
            dx=f/df;
            double temp=rts;
            rts -= dx;
            if (temp == rts) return rts;
        }
        if (abs(dx) < xacc) return rts; //Convergence criterion.
        double f=Rx(Ar,Br,rho,rts,1);//funcd(rts);
        double df=dRx(Ar,Br,rts); //funcd.df(rts);
        //The one new function evaluation per iteration.
        if (f < 0.0) xl=rts; //Maintain the bracket on the root.
        else xh=rts;
    }
    throw("Maximum number of iterations exceeded in rtsafe");
}

void mapping_rho_to_R(double Ar,double Br,double x1,double x2, vec1d& R,vec1d& rho,int iRmax){
    //int iRmax=size(rho);
    for (int iR=0; iR<iRmax; iR++){
        R[iR]=rtsafe(Ar,Br,rho[iR],x1,x2,1.e-14);
    }
    ofstream fp_output_R; fp_output_R.open("print_R_rho.txt");
    for (int iR=0; iR<iRmax; iR++){
        fp_output_R.precision(5);
        fp_output_R << rho[iR] << " " << R[iR] << endl;
    }
    fp_output_R.close();
}


//FUNCTION TO READ A MO FROM A CUBE FILE
void loading_MO (ifstream& fp_input,int& Nx, int& Ny, int& Nz, vec1d& Origin,vec1d& spacings, vec2d& Geometry, vec3d& MO,string program)
{
    if (!fp_input.is_open())
    {
        cout << "error opening cube file " << endl;
        exit(1);
    }
    
    //printf("Reading cube file...\n");
    string trash;
    double value;
    for(int i=0; i<2; i++) getline(fp_input, trash);
    //reading number of atoms and reference origin
    int iatoms=0; fp_input >> iatoms;
    if(program=="Molpro") iatoms=-iatoms;
    printf("Number of atoms: %d\n",iatoms);
    for(int i=0; i<3; i++)
    {
        fp_input >> value;
        Origin[i]=value;
    }
    //reading Nx, Ny, and Nz + the spacings in each direction
    fp_input >> Nx >> spacings[0] >> value >> value;
    fp_input >> Ny >> value >> spacings[1] >> value;
    fp_input >> Nz >> value >> value >> spacings[2];
    
    printf("Origin cube file: (%3.6f,%3.6f,%3.6f)\n#NX: %d, #NY: %d, #NZ: %d \nSpacings, dx: %5.3f, dy: %5.3f, dz: %5.3f au\n",Origin[0],Origin[1],Origin[2],Nx,Ny,Nz,spacings[0],spacings[1],spacings[2]);

    Geometry.resize(iatoms,vec1d(5,0.));
    for(int i=0; i<iatoms; i++)
    {
        
        for(int j=0; j<5; j++)
        {
            fp_input >> value;
            Geometry[i][j]=value;
            //Geometry[i].push_back(double(value));
        }
    }
    
    if(program=="Molpro") fp_input >> value >> value;
    //cout << Geometry[0][0] << " " << Geometry[iatoms-1][3] << " " << Geometry[iatoms-1][4] << endl;
    //printf("Reading MO...\n");

    int count=0;
    //MO.reserve(Nx);
    MO.resize(Nx,vec2d(Ny,vec1d(Nz,0.)));
    for(int ix=0; ix<Nx; ix++)
    {
        //MO[ix].reserve(Ny);
        for(int iy=0; iy<Ny; iy++)
        {
            //MO[ix][iy].reserve(Nz);
            for(int iz=0; iz<Nz; iz++)
            {
                fp_input >> MO[ix][iy][iz];
                //if (count<5) { cout << ix << " " << iy << " " << iz << " " << MO[ix][iy][iz] << " " << count << endl; }
                count++;
            }//end iz
        }//end iy
    }//end ix
    
    //printf("Read MO...\n");
    fp_input.close();
}


//FUNCTION TO PERFORM SPHERICAL HARMONIC EXPANSION
void Do_Single_Center (int Nx,int Ny,int Nz,int Lmax,int mmax,int orderPolynomial1, int orderPolynomial2,int iRj, vec1d& R, vec1d& Origin,vec1d& spacings, vec3d& MO3D,vec2x& Plm,string filename)
{
    string s= "Performing spherical harmonic expansion: " + filename + "\n";
    printf(s.c_str());
    int iRmax=R.size();
    
    //%%%%%%%% Interpolating MO
    /* //using likely library
    vec1d MOToFit1D(Nx*Ny*Nz,0.);
    for (int ix=0; ix<Nx; ix++)
    {
        for (int iy=0; iy<Ny; iy++)
        {
            for (int iz=0; iz<Nz; iz++)
            {
                MOToFit1D[ix+Nx*(iy + Ny*iz)]=MO3D[ix][iy][iz];
            }
        }
    }
    likely::TriCubicInterpolator MO(MOToFit1D,spacings,Origin,Nx,Ny,Nz);
    */
    
    FittedDataREAL MO(MO3D,spacings,Origin);
    
    for (int ix=0; ix<Nx; ix++) {
        for (int iy=0; iy<Ny; iy++) {
            for (int iz=0; iz<Nz; iz++) {
                double x=Origin[0]+spacings[0]*ix;
                double y=Origin[1]+spacings[1]*iy;
                double z=Origin[2]+spacings[2]*iz;
                if(abs(MO(x,y,z)-MO3D[ix][iy][iz])>1.e-10){
                    printf("WARNING: Interpolated MO is not good: x %3.3f y %3.3f z %3.3f MO3D %2.10f MO %2.10f\n",x,y,z,MO3D[ix][iy][iz],MO(x,y,z));
                    //cout << "x " << x << " y " << y << " z " << z << " MO3D " << MO3D[ix][iy][iz] << " MO " << MO(x,y,z) << endl;
                }
            }
        }
    }
    
    
    //%%%%%%%% Calculation of Plm for the core orbitals
    vec1d root1(orderPolynomial1,0.);
    vec1d root2(orderPolynomial2,0.);
    vec1d weight1(orderPolynomial1,0.);
    vec1d weight2(orderPolynomial2,0.);
    gauleg(0,pi,root1,weight1); gauleg(0,2.*pi,root2,weight2);
    double *Ytemp = new double[gsl_sf_legendre_array_n(Lmax)];
    
    /*
    gsl_sf_legendre_array(GSL_SF_LEGENDRE_SPHARM,Lmax,cos(pi/2.),Ytemp);
    for (int l=0; l<=Lmax; l++)
        for (int m=-l; m<=l; m++){
            int Condon_Shortley=(m%2!=0 && m>0 ? -1 : 1);
            cout << "l " << l << " m " << m << " " << Condon_Shortley*Ytemp[gsl_sf_legendre_array_index(l,abs(m))] << endl;
        }
    */
    for (int iR=0; iR<iRmax; iR++) {
        if (iR<iRj)
        {
            for (int l=0; l<=Lmax; l++)
            {
                int M1=( l<=mmax ? l : mmax);
                for (int m=-M1; m<=M1; m++)
                {
                    for (int i=0; i<orderPolynomial1; i++)
                    {
                        double wei1=weight1[i];
                        double droot1=root1[i];
                        gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_SPHARM,Lmax,cos(root1[i]),-1,Ytemp);
                        complexd alpha2=0.;
                        for (int j=0; j<orderPolynomial2; j++)
                        {
                            //complexd temp=boost::math::spherical_harmonic(l,-m,root1[i],root2[j]);
                            double wei2=weight2[j];
                            double x=R[iR]*sin(root1[i])*cos(root2[j]);
                            double y=R[iR]*sin(root1[i])*sin(root2[j]);
                            double z=R[iR]*cos(root1[i]);
                            int Condon_Shortley=(m%2!=0 && m>0 ? -1 : 1);
                            complexd temp=(Condon_Shortley*Ytemp[gsl_sf_legendre_array_index(l,abs(m))])*exp(-c1*(m*root2[j]));
                            double sign=( m%2==0 ? 1. : -1.);
                            alpha2+=temp*(MO(x,y,z)*wei2*sign); //we multiply by complex conjugate Y_l^m = Y_l^-m
                        }
                        Plm[SpH(l,m,mmax)][iR]+=alpha2*(wei1*sin(droot1)*R[iR]);
                    }
                } //--- End m loop
            } //--- End l loop
        }
    }//--- End R loop
    delete [] Ytemp;
    
    //system("mkdir Plm_orbitals");
    string sname="Plm_orbitals/Plm_"+filename;
    ofstream fp_output; fp_output.open(sname.c_str());
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
        if (iR<iRj)
        {
            fp_output << R[iR] << " ";
            for (int l=0; l<=Lmax; l++)
            {
                int M1=( l<=mmax ? l : mmax);
                for (int m=-M1; m<=M1; m++)
                {
                    fp_output << Plm[SpH(l,m,mmax)][iR] << " ";
                }//--- End m loop
            }//--- End l loop
            fp_output << endl;
        }
    }//--- End R loop
    
    fp_output.close();
}

