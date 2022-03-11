#pragma omp parallel for schedule(static)
for(int k=0; k<TMO; k++){
    vec3d MO; //array to store the molecular orbital
    vec1d Origin; //array to save the origin
    vec2d Geometry_pr; //store geometry molecule
    vec1d CubeGrid_spacings(3,0.); //array to store the spacing of the cube file
    vec1d CubeGrid_origin(3,0.); //array to store the reference origin of the cube file
    string s;
    
    stringstream iname;
    vec1i ij; ij=Sy.ij(k);
    
    iname.seekp(0,ios::beg); iname << (ij[0]+1);
    int Nx,Ny,Nz; Nx=0; Ny=0; Nz=0;
    
    stringstream jname;
    jname.seekp(0,ios::beg); jname << (ij[1]+1);
    s= scube + "_orbital_" + jname.str() + "." + iname.str() + ".cube";
    printf("Reading: %s\n",s.c_str());
    //READING CUBE FILE
    ifstream f;
    f.open(s.c_str());
    loading_MO(f,Nx,Ny,Nz,CubeGrid_origin,CubeGrid_spacings,Geometry_pr,MO,program);
    if (k==0) {
        Geometry.resize(Geometry_pr.size(),vec1d(5,0.));
        for(int i=0; i<Geometry_pr.size(); i++)
        {
            for(int j=0; j<5; j++) Geometry[i][j]=Geometry_pr[i][j];
        }
    }
    f.close();
    
    //PERFORMING SINGLE CENTER EXPANSION
    s= "orbital_" + jname.str() + "." + iname.str() + ".cube";
    Do_Single_Center (Nx,Ny,Nz,Lmax,mmax,orderPolynomial1,orderPolynomial2,iRj,R,CubeGrid_origin,CubeGrid_spacings,MO,Plm[k],s);
}
