
scalar DissPower(0);
scalar DissPowerbase (0);
scalar V(0);
scalar dvgnum(0);
scalar dvgchange(0);
scalar jmixdegree(0);
scalar jnum(0);

// initialize the residual in the solver
scalar loopnumU(0),loopnumUa(0), loopnumUc(0);
scalar Uresidual0(10.0), Uaresidual0(10.0), Ucresidual0(10.0);
vector Uresidual = vector(10,10,10);
vector Uaresidual = vector(10,10,10);
vector Ucresidual = vector(10,10,10);

scalar cumulativeContErra = 0;
scalar cumulativeContErrc = 0;

int i=0, allcells, nallcells;
dimensionedScalar Rmin = Rmesh/(2*sqrt(3));
int iterTO=1; 
double fsensmax=0; 
double gsensmax=0;
double beta=1;
double eta5=0.5;
int betaloop=0;

scalar allcell=dvg.size();
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

double len=std::pow(mesh.V()[0],1.0/3);

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

Vec xmma,xold,xmin,xmax;
PetscScalar *xmma_array;

VecCreate(PETSC_COMM_WORLD,&xmma);
VecSetSizes(xmma,n_dof,N_dof);
VecSetFromOptions(xmma);
VecGetArray(xmma,&xmma_array);

PetscPrintf(PETSC_COMM_WORLD,"%i\n",xmma_array[0]);
PetscPrintf(PETSC_COMM_WORLD,"%i\n",xmma_array[1]);

for(i=0;i<allcells;i++)
{
    xmma_array[i]=dvg[i];
}
VecRestoreArray(xmma,&xmma_array);

VecDuplicate(xmma,&xold); // VecDuplicate: set the value of xold as xmma
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

// sum of dvg at cells_design
setCells(dvg2, cells_fluid, 0);
setCells(dvg2, cells_design, 1);
for (i=0;i<allcells;i++)
{
    dvgnum = dvg2[i] + dvgnum;
}
reduce(dvgnum, sumOp<scalar>());
Info<< "dvgnum" << dvgnum << endl;

// sum of Cv of all mesh at the outlet 
for (inum=0; inum<nMIPatch; inum++)
{
    jnum = jnum + sum(magSqr(Cv.boundaryField()[MIPatchList[inum]]));
}
reduce(jnum,sumOp<scalar>());
Info<< "jnum" << jnum << endl;

// change value of dvg
for(i=0;i<allcells;i++)
{
    dvg0[i]=dvgHea[i];
}
