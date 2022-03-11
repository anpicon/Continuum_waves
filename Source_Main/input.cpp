ifstream fp_input;

//%%%%%%%%%%%%%%%%%%% READING INPUT
if ( argc != 2 )
{
    cout<<"usage: "<< argv[0] <<" <filename>" << endl;
    exit(1);
}
else
{
    fp_input.open(argv[1]);
    if (!fp_input.is_open())
    {
        cout << "error opening file " << argv[1] << endl;
        exit(1);
    }
}
    
system("mkdir Plm_orbitals");
string sname;
string program; fp_input >> sname >> program; //Program to make electronic structure
double epsilon; fp_input >> sname >> epsilon; //Energy of the electron in the continuum; Ip-Eci[ek];
double Rmax; fp_input >> sname >> Rmax; //Max grid for ionization calculations
double Rj; fp_input >> sname >> Rj; //Max grid for MO expanded in spherical harmonics
int Lmax; fp_input >> sname >> Lmax; //Max L
int mmax; fp_input >> sname >> mmax; //Max M
int orderPolynomial1; fp_input >> sname >> orderPolynomial1; //for angular integration, max order-polynomial in theta
int orderPolynomial2; fp_input >> sname >> orderPolynomial2; //for angular integration, max order-polynomial in varphi
double dh; fp_input >> sname >> dh; //for defining radial grid
double Ar; fp_input >> sname >> Ar; //for defining radial grid
double Br; fp_input >> sname >> Br; //for defining radial grid
string TheoryLevel; fp_input >> sname >> TheoryLevel; //defining level calculation
string PointGroup; fp_input >> sname >> PointGroup; //defining point group
Symmetry Sy(PointGroup);
fp_input >> sname;
for (std::vector<int>::size_type i=0; i<(Sy.NMO).size(); i++) {
    fp_input >> Sy.NMO[i];
}
int hsym=0;
int hocc=0;
fp_input >> sname;
for (std::vector<int>::size_type i=0; i<(Sy.NMO).size(); i++) {
    int A1=0;
    fp_input >> A1;
    if (A1>Sy.NMO[i]) {
        printf("Error in input, hole orbital: %d and max. number orbital %d\n",A1,Sy.NMO[i]);
        exit(1);
    }
    if(A1!=0){
        hsym=i;
        hocc=A1;
    }
}
int hspin=0; fp_input >> sname >> sname;
if(sname=="a") hspin=0;
else if(sname=="b") hspin=1;
else {
    printf("Error in input, hole spin: %s\n",sname.c_str());
    exit(1);
}
double BE; fp_input >> sname >> BE; //Binding energy of the cation
string scube; fp_input >> scube >> scube; //filename of the cube files
string pathocc; fp_input >> pathocc >> pathocc; //filename of the occupation file
fp_input.close();
