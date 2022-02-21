
Info<<"sensitivity analysis"<<endl;

// alphadao=-(alphasolidmax-alphasolidmin)*asqa*(asqa+1)/((asqa+dvgHea)*(asqa+dvgHea));

fsens=-(alphasolidmax-alphasolidmin)*asqa*(asqa+1)/((asqa+dvgHea)*(asqa+dvgHea))*(U & Uc);
gsens=-(alphasolidmax-alphasolidmin)*asqa*(asqa+1)/((asqa+dvgHea)*(asqa+dvgHea))*(U & Ua);


fsensmax=gMax(mag(fsens.primitiveField()));
gsensmax=gMax(mag(gsens.primitiveField()));
fsensm=fsens/fsensmax;
gsensm=gsens/gsensmax;

setCells(fsens, cells_fluid, 0);
setCells(gsens, cells_fluid, 0);


if(beta<512&&opt>100&&betaloop>100&&dvgchange<1e-3)
{
    beta=2*beta;
    betaloop=0;
}

betaloop++;

for(i=0;i<allcells;i++)//ith node
{
    if(dvgPDE[i]<=eta5)
    {
        drho[i]=beta*Foam::exp(-beta*(1-dvgPDE[i]/eta5))+Foam::exp(-beta); 
    }
    else
    {
        drho[i]=beta*Foam::exp(-beta*(dvgPDE[i]-eta5)/(1-eta5))+Foam::exp(-beta); 
    }
}

fsensHea=fsensm*drho;
gsensHea=gsensm*drho;


fvScalarMatrix Eqnf (-fvm::laplacian(Rmin*Rmin,fsensPDE)+fvm::Sp(1,fsensPDE)-fsensHea);
Eqnf.solve();

fvScalarMatrix Eqng (-fvm::laplacian(Rmin*Rmin,gsensPDE)+fvm::Sp(1,gsensPDE)-gsensHea);
Eqng.solve();

gx[0]=(PetscScalar)(jmixdegree/J1-1.0);

VecGetArray(dfdx,&dfdx_array);
VecGetArrays(dgdx,m,&dgdx_array);

for(i=0;i<allcells;i++)
{
    dfdx_array[i]=fsensHea[i];
    dgdx_array[0][i]=gsensHea[i];
}
VecRestoreArray(dfdx,&dfdx_array);
VecRestoreArrays(dgdx,m,&dgdx_array);

mma->SetOuterMovelimit(Xmin,Xmax,movlim,xmma,xmin,xmax);

mma->Update(xmma,dfdx,gx,dgdx,xmin,xmax);

VecGetArray(xmma,&xmma_array);

// dvg old
for(i=0;i<allcells;i++)
{
    dvg0[i]=dvgHea[i];
}

for(i=0;i<allcells;i++)
{
    dvg[i]=xmma_array[i];
}

VecRestoreArray(xmma,&xmma_array);


fvScalarMatrix Eqndvg (-fvm::laplacian(Rmin*Rmin,dvgPDE)+fvm::Sp(1,dvgPDE)-dvg);
Eqndvg.solve();

for(i=0;i<allcells;i++)
{
   if(dvgPDE[i]<=eta5)
   {
     dvgHea[i]=eta5*(Foam::exp(-beta*(1.0-dvgPDE[i]/eta5))-(1.0-dvgPDE[i]/eta5)*Foam::exp(-beta)); 
   }
   else
   {
     dvgHea[i]=eta5+(1-eta5)*(1-Foam::exp(-beta*(dvgPDE[i]-eta5)/(1.0-eta5))+(dvgPDE[i]-eta5)*Foam::exp(-beta)/(1-eta5));
   }   
}

setCells(dvgHea, cells_fluid, 1);

alphasolidmax=alphasolidini*Foam::pow(alphaMaxdelta.value(),(opt/100));

alphasolidmax.value()=Foam::min(alphasolidmax.value(),alphaMax.value());

if(opt>=100)
{
    asqa=asqaini+(opt-100)*asqadelta;
}
asqa=Foam::min(asqa,asqamax);


alphasolid = alphasolidmin+(alphasolidmax-alphasolidmin)*asqa*(1-dvgHea)/(asqa+dvgHea);

opt++;
