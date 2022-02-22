// 
    std::remove("Voluse.txt");
    std::remove("DissPower.txt");
    std::remove("DissPowerbase.txt");
    std::remove("Time.txt");
    std::remove("Mixdegree.txt");
    std::remove("opt.txt");
    std::remove("dvgchange.txt");
    std::remove("jnum.txt");
    std::remove("J2.txt");
    std::remove("alphasolidmax.txt");
    std::remove("asqa.txt");
    std::remove("resi1.txt");
    std::remove("resi2.txt");
    std::remove("resi3.txt");
    std::remove("resif.txt");
    std::remove("resig.txt");
    std::remove("residvg.txt");
    std::remove("fsensmax.txt");
    std::remove("gsensmax.txt");
    std::remove("beta.txt");


    scalar C(0);
    scalar DissPower(0);

    scalar DissPowerbase (0);
    scalar V(0);
    scalar dvgnum(0);
    scalar dvgchange(0);
    scalar jmixdegree(0);
    scalar jnum(0);
    scalar J2(0);

    scalar resi1(0),resi2(0),resi3(0);

    scalar Uresidualsave(1.0),Uaresidualsave(1.0),Ucresidualsave(1.0);

    vector Uresidual,Uaresidual,Ucresidual;

    int i=0,j=0,k=0,jk,allcells,nallcells,mycells=0,opt=1;
    double fsensmax=0; double gsensmax=0;
    double beta=1;
    double eta5=0.5,del=10;
    int betaloop=0;
    scalar allcell=dvg.size();//allcell:length of dvg in the i-th core
    allcells=static_cast<int>(allcell); 
    List<scalar> cellsize(Pstream::nProcs());
    for(i=0;i<Pstream::nProcs();i++)
    {
       cellsize[i]=0;
    }

    cellsize[Pstream::myProcNo()]=allcells;
    reduce(cellsize, sumOp<List<scalar>>()); 
    reduce(allcell, sumOp<scalar>()); 
    nallcells=static_cast<int>(allcell);//nallcells:total length of dvg

    for(i=0;i<Pstream::myProcNo();i++)
    {
       mycells=mycells+cellsize[i];
    }

    double add1=0,add2=0,len=std::pow(mesh.V()[0],1.0/3);


    Info<< "nallcells" << nallcells << endl;
    Info<< "allcells" << allcells << endl;

    List<vectorField> AllcellList(Pstream::nProcs());
    vectorField localcell(mesh.C());
    AllcellList[Pstream::myProcNo()] = localcell;
    Pstream::gatherList(AllcellList);
    Pstream::scatterList(AllcellList);
    vectorField Allcellcenter 
    (
        ListListOps::combine<vectorField>
        (
            AllcellList,
            accessOp<vectorField>()
        )
    );  


//***********************************//Petsc
PetscInt n_dof,N_dof,m=1;//m:number of constraints
n_dof=allcells;
N_dof=nallcells;

PetscScalar Xmin=0; // Min. value of design variables
PetscScalar Xmax=1; // Max. value of design variables

Vec xmma,xold,xmin,xmax,xTilde,xPhys;
PetscScalar *xmma_array;
PetscScalar *xmin_array;


VecCreate(PETSC_COMM_WORLD,&xmma);
VecSetSizes(xmma,n_dof,N_dof);
VecSetFromOptions(xmma);
VecGetArray(xmma,&xmma_array);



PetscPrintf(PETSC_COMM_WORLD,"%i\n",xmma_array[0]);
PetscPrintf(PETSC_COMM_WORLD,"%i\n",xmma_array[1]);


VecDuplicate(xmma,&xTilde);


VecDuplicate(xmma,&xPhys);

for(i=0;i<allcells;i++)
{
  xmma_array[i]=dvg[i];
}
VecRestoreArray(xmma,&xmma_array);

VecDuplicate(xmma,&xold); //向量赋值
VecDuplicate(xmma,&xmin);
VecDuplicate(xmma,&xmax); 

PetscScalar aMMA[m];
PetscScalar cMMA[m];
PetscScalar dMMA[m];

for (PetscInt i=0;i<m;i++)
{
  aMMA[i]=0.0;
  dMMA[i]=1.0;
  cMMA[i]=1000.0;
}


MMA *mma;

mma=new MMA(N_dof,m,xmma, aMMA, cMMA, dMMA, mma_init, mma_dec, mma_inc);
//**********************************************//dfdx
Vec dfdx; // Sensitivities of objective
VecDuplicate(xmma,&dfdx); 
PetscScalar *dfdx_array;
VecGetArray(dfdx,&dfdx_array);
VecRestoreArray(dfdx,&dfdx_array); 
//*****************************************************//gx
PetscScalar *gx= new PetscScalar[m]; // Array with constraint values

//*****************************************************//dgdx
Vec *dgdx; // Sensitivities of constraints (vector array)
VecDuplicateVecs(xmma,m, &dgdx); 
PetscScalar **dgdx_array;
VecGetArrays(dgdx,m,&dgdx_array);
VecRestoreArrays(dgdx,m,&dgdx_array);

alphasolid = alphasolidmin+(alphasolidmax-alphasolidmin)*asqa*(1-dvgHea)/(asqa+dvgHea);


setCells(dvg2, cells_fluid, 0);
setCells(dvg2, cells_design, 1);

for (i=0;i<allcells;i++)
{
    dvgnum=dvg2[i]+dvgnum;
}
reduce(dvgnum, sumOp<scalar>());
Info<< "dvgnum" << dvgnum << endl;
    
jnum=0;
for (iLoop=0; iLoop<nObjPatch; iLoop++)
{
    jnum =jnum + sum(magSqr(alphad.boundaryField()[objPatchList[iLoop]]));
}
reduce(jnum,sumOp<scalar>());
Info<< "jnum" << jnum << endl;

for(i=0;i<allcells;i++)
{
    dvg0[i]=dvgHea[i];
}



