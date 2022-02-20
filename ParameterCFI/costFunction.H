
for(i=0;i<allcells;i++)
{
    dvg1[i]=dvgHea[i];
}

setCells(dvg1, cells_fluid, 0);

V=0;
for (i=0;i<allcells;i++)
{
    V=dvg1[i]+V;
}
reduce(V, sumOp<scalar>());
V=V/dvgnum;
 

dvgchange=0;
for (i=0;i<allcells;i++)
{
    dvgchange=dvgchange+Foam::mag(dvg0[i]-dvgHea[i]);
}
reduce(dvgchange, sumOp<scalar>());
dvgchange=dvgchange/nallcells;

//DissPower
DissPower=0;
for (i=0; i<nObjPatchpow; i++)
{
    DissPower =DissPower - sum(
    phi.boundaryField()[conPatchList[i]]*(p.boundaryField()[conPatchList[i]]
    + 0.5*magSqr(U.boundaryField()[conPatchList[i]])));
}
 
if(geo_dim==2) 
{
    DissPower=DissPower/len;
}

reduce(DissPower, sumOp<scalar>()); 
 
DissPowerbase=DissPower/D_normalization;
Info<<"constraint function (DissPowerbase): "<<DissPowerbase<<endl; 
  
jmixdegree=0;
for (iLoop=0; iLoop<nObjPatch; iLoop++)
{
        jmixdegree =
            jmixdegree + sum(
            magSqr(alphac.boundaryField()[objPatchList[iLoop]]
            -alphad.boundaryField()[objPatchList[iLoop]])
            );
}

reduce(jmixdegree, sumOp<scalar>());
jmixdegree=jmixdegree/jnum;
Info<<"Objective function : "<<jmixdegree<<endl;

if (Pstream::master())
{
ofstream outfile1("Voluse.txt",std::ios::app);
outfile1<<V<<"\n";
outfile1.close();

ofstream outfile2("DissPower.txt",std::ios::app);
outfile2<<DissPower<<"\n";
outfile2.close();
   
ofstream outfile3("DissPowerbase.txt",std::ios::app);
outfile3<<DissPowerbase<<"\n";
outfile3.close();

ofstream outfile4("J2.txt",std::ios::app);
outfile4<<J2<<"\n";
outfile4.close();

ofstream outfile5("Time.txt",std::ios::app);
outfile5<<runTime.elapsedClockTime()<<"\n";
outfile5.close(); 

ofstream outfile6("Mixdegree.txt",std::ios::app);
outfile6<<jmixdegree<<"\n";
outfile6.close(); 

ofstream outfile7("opt.txt",std::ios::app);
outfile7<<opt<<"\n";
outfile7.close();

ofstream outfile8("dvgchange.txt",std::ios::app);
outfile8<<dvgchange<<"\n";
outfile8.close();

ofstream outfile9("alphasolidmax.txt",std::ios::app);
outfile9<<alphasolidmax.value()<<"\n";
outfile9.close();

ofstream outfile10("asqa.txt",std::ios::app);
outfile10<<asqa<<"\n";
outfile10.close();

ofstream outfile11("resi1.txt",std::ios::app);
outfile11<<resi1<<"\n";
outfile11.close();
ofstream outfile12("resi2.txt",std::ios::app);
outfile12<<resi2<<"\n";
outfile12.close();
ofstream outfile13("resi3.txt",std::ios::app);
outfile13<<resi3<<"\n";
outfile13.close();

ofstream outfile14("fsensmax.txt",std::ios::app);
outfile14<<fsensmax<<"\n";
outfile14.close();

ofstream outfile15("gsensmax.txt",std::ios::app);
outfile15<<gsensmax<<"\n";
outfile15.close();

ofstream outfile16("beta.txt",std::ios::app);
outfile16<<beta<<"\n";
outfile16.close();
}
