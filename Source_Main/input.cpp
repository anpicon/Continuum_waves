ifstream fp_input;
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
//string program="Molpro";
string sname;
//double epsilon; fp_input >> sname >> epsilon; //Energy of the electron in the continuum; Ip-Eci[ek];
string program; fp_input >> sname >> program; //Name of the code used to calculate electronic structure
double omega; fp_input >> sname >> omega; //Photon energy
double Rmax; fp_input >> sname >> Rmax; //Max grid for ionization calculations
double Rj; fp_input >> sname >> Rj; //Max grid for MO expanded in spherical harmonics
int Lmax; fp_input >> sname >> Lmax; //Max L
int mmax; fp_input >> sname >> mmax; //Max M
int orderPolynomial1; fp_input >> sname >> orderPolynomial1; //for angular integration, max order-polynomial in theta
int orderPolynomial2; fp_input >> sname >> orderPolynomial2; //for angular integration, max order-polynomial in varphi
double dh; fp_input >> sname >> dh; //for defining radial grid
double Ar; fp_input >> sname >> Ar; //for defining radial grid
double Br; fp_input >> sname >> Br; //for defining radial grid
string TheoryLevel; fp_input >> sname >> TheoryLevel; //defining point group
string PointGroup; fp_input >> sname >> PointGroup; //defining point group
Symmetry Sy(PointGroup);
fp_input >> sname;
for (std::vector<int>::size_type i=0; i<(Sy.NMO).size(); i++) {
    fp_input >> Sy.NMO[i];
}
/*
int irrerep=0;
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
vec1i NMO(irrerep,0); fp_input >> sname;
for (std::vector<int>::size_type i=0; i<NMO.size(); i++) {
    fp_input >> NMO[i];
}
*/
string scube; fp_input >> scube >> scube; //filename of the cube files
string pathocc; fp_input >> pathocc >> pathocc; //filename of the occupation file
string pathBE; fp_input >> pathBE >> pathBE; //filename of the binding energy file
fp_input.close();
