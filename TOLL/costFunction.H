
for(i=0; i<allcells; i++)
{
    dvg1[i]=dvgHea[i];
}

setCells(dvg1, cells_fluid, 0);

V=0;
for (i=0; i<allcells; i++)
{
    V=dvg1[i]+V;
}
reduce(V, sumOp<scalar>());
V=V/dvgnum;

dvgchange=0;
for (i=0; i<allcells; i++)
{
    dvgchange = dvgchange + Foam::mag(dvg0[i]-dvgHea[i]);
}
reduce(dvgchange, sumOp<scalar>());
dvgchange = dvgchange / nallcells;

DissPower=0;
for (i=0; i<nPOWERPatch; i++)
{
        DissPower = DissPower - sum(
        phi.boundaryField()[POWERPatchList[i]]*(p.boundaryField()[POWERPatchList[i]]
        + 0.5*magSqr(U.boundaryField()[POWERPatchList[i]]))); 
}
 
if(geo_dim==2) 
{
    DissPower = DissPower / len;
}
DissPowerbase = DissPower / D_normalization; 
reduce(DissPower, sumOp<scalar>()); 

jmixdegree=0;
for (i=0; i<nMIPatch; i++)
{
        jmixdegree =
            jmixdegree + sum(
            magSqr(C.boundaryField()[MIPatchList[i]]
            -Cv.boundaryField()[MIPatchList[i]])
            );
}

reduce(jmixdegree, sumOp<scalar>());
jmixdegree = jmixdegree / jnum;

if (Pstream::master())
{
    ofstream outfile("OptMessage.txt",std::ios::app);
    outfile<< iterTO << " " << jmixdegree << " " << DissPowerbase << " " 
        << dvgchange << " " <<  V << " " 
        << Uresidual0 << " " << Uaresidual0 << " "<< Ucresidual0 << " " 
        << loopnumU << " " << loopnumUa <<" "<< loopnumUc << " "
        << alphasolidmax.value() << " " << beta << " " << asqa << " "  
        << "\n";
    outfile.close();

    ofstream outfile6("Mixdegree.txt",std::ios::app);
    outfile6<<jmixdegree<<"\n";
    outfile6.close(); 

}

iterTO++;
