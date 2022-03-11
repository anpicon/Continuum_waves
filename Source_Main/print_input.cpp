printf("Program used: %s\n",program.c_str());
printf("Theory level: %s\n",TheoryLevel.c_str());
printf("Point group: %s\n",PointGroup.c_str());
printf("MOs: ");
for (std::vector<int>::size_type i=0; i<Sy.NMO.size(); i++) {
    printf("%d ",Sy.NMO[i]);
}
printf("\n");

/*
printf("Hole: ");
for (std::vector<int>::size_type i=0; i<NMO.size(); i++) {
    if (i==hsym) {
        printf("%d ",hocc);
    } else {
        printf("%d ",0);
    }
}
printf("\nSpin: %d\n",hspin);
printf("\nBinding energy after the hole (eV): %5.3f\n",BE*energy_au_eV);
*/

printf("\nPhoton energy: %3.3f au, %3.3f eV\n",omega,omega*energy_au_eV);
printf("One-center spatial grid:\n");
printf("Parameters for angular grid\n");
printf("Max L: %i\n",Lmax);
printf("Max M: %i\n",mmax);
printf("Polynomial order for angular integration, theta: %i\n",orderPolynomial1);
printf("Polynomial order for angular integration, varphi: %i\n\n",orderPolynomial2);
printf("Parameters for radial grid\n");
printf("Max radial distance: %5.2f au\n",Rmax);
printf("Max radial for bound MO: %5.2f au\n",Rj);
printf("dh (rho vector): %3.5f\n",dh);
printf("Ar: %5.3f au\n",Ar);
printf("Br: %5.3f au\n",Br);
