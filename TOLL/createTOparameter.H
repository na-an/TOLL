IOdictionary optProperties
(
   IOobject
   (
       "optProperties",
       runTime.constant(),
       mesh,
       IOobject::MUST_READ,
       IOobject::NO_WRITE
   )
);

dimensionedScalar alphasolidini(optProperties.lookup("alphasolidini"));
dimensionedScalar alphaMaxdelta(optProperties.lookup("alphaMaxdelta"));
dimensionedScalar alphasolidmax(optProperties.lookup("alphasolidmax"));
dimensionedScalar alphaMax(optProperties.lookup("alphaMax"));

dimensionedScalar alphasolidmin(optProperties.lookup("alphasolidmin"));

dimensionedScalar Rmesh(optProperties.lookup("Rmesh"));

scalar V0(readScalar(optProperties.lookup("V0")));

scalar raa0(readScalar(optProperties.lookup("raa0")));
scalar mma_init(readScalar(optProperties.lookup("mma_init")));
scalar mma_dec(readScalar(optProperties.lookup("mma_dec")));
scalar mma_inc(readScalar(optProperties.lookup("mma_inc")));
scalar movlim(readScalar(optProperties.lookup("movlim")));

scalar D_normalization(readScalar(optProperties.lookup("D_normalization")));

scalar J1(readScalar(optProperties.lookup("J1")));

scalar geo_dim(readScalar(optProperties.lookup("geo_dim")));

scalar asqa(readScalar(optProperties.lookup("asqa")));
scalar asqaini(readScalar(optProperties.lookup("asqaini")));
scalar asqadelta(readScalar(optProperties.lookup("asqadelta")));
scalar asqamax(readScalar(optProperties.lookup("asqamax")));

scalar MaxUloop(readScalar(optProperties.lookup("MaxUloop")));
scalar MaxUaloop(readScalar(optProperties.lookup("MaxUaloop")));
scalar MaxUcloop(readScalar(optProperties.lookup("MaxUcloop")));
scalar Uresidualset(readScalar(optProperties.lookup("Uresidualset")));
scalar Ucresidualset(readScalar(optProperties.lookup("Ucresidualset")));
scalar Uaresidualset(readScalar(optProperties.lookup("Uaresidualset")));

word zoneName1 = "zone_fluid"; 
label zoneID1 = mesh.cellZones().findZoneID(zoneName1);
labelList cells_fluid = mesh.cellZones()[zoneID1];

word zoneName2 = "zone_design"; 
label zoneID2 = mesh.cellZones().findZoneID(zoneName2);
labelList cells_design = mesh.cellZones()[zoneID2];

const labelList& INLET1Cells = mesh.boundary()["INLET1"].faceCells();
const labelList& INLET2Cells = mesh.boundary()["INLET2"].faceCells();

volScalarField dvg
(
    IOobject
    (
        "dvg",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    scalar(V0),
    zeroGradientFvPatchScalarField::typeName 
);

setCells(dvg, cells_fluid, 1);

volScalarField dvgPDE
(
    IOobject
    (
        "dvgPDE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    scalar(V0),
    zeroGradientFvPatchScalarField::typeName 
);

volScalarField dvgHea
(
    IOobject
    (
        "dvgHea",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    scalar(V0),
    zeroGradientFvPatchScalarField::typeName 
);

volScalarField dvg0(dvg);
volScalarField dvg1(dvg);
volScalarField dvg2(dvg);
volScalarField drho(dvg);

volScalarField alphasolid
(
    IOobject
    (
        "alphasolid",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    alphasolidmin+(alphasolidmax-alphasolidmin)*asqa*(1-dvgHea)/(asqa+dvgHea)
);

List<scalarField> AllrList(Pstream::nProcs());
scalarField localr(dvg.internalField());
AllrList[Pstream::myProcNo()] = localr;
Pstream::gatherList(AllrList);
Pstream::scatterList(AllrList);
scalarField Alldvg
(
    ListListOps::combine<scalarField>
    (
        AllrList,
        accessOp<scalarField>()
    )
);

volScalarField fsens
(
    IOobject
    (
        "fsens",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (alphasolidmax-alphasolidmin)*asqa*(1-dvgHea)/(asqa+dvgHea)*(Ua & U),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField fsensm(fsens);

volScalarField fsensPDE
(
    IOobject
    (
        "fsensPDE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (alphasolidmax-alphasolidmin)*asqa*(1-dvgHea)/(asqa+dvgHea)*(Ua & U),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField fsensHea
(
    IOobject
    (
        "fsensHea",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (alphasolidmax-alphasolidmin)*asqa*(1-dvgHea)/(asqa+dvgHea)*(Ua & U),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField gsens
(
    IOobject
    (
        "gsens",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (alphasolidmax-alphasolidmin)*asqa*(1-dvgHea)/(asqa+dvgHea)*(Uc & U),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField gsensm(gsens);

volScalarField gsensPDE
(
    IOobject
    (
        "gsensPDE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (alphasolidmax-alphasolidmin)*asqa*(1-dvgHea)/(asqa+dvgHea)*(Uc & U),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField gsensHea
(
    IOobject
    (
        "gsensHea",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (alphasolidmax-alphasolidmin)*asqa*(1-dvgHea)/(asqa+dvgHea)*(Uc & U),
    zeroGradientFvPatchScalarField::typeName
);

dictionary power_dissipation=mesh.solutionDict().subDict("power_dissipation");
int nPOWERPatch=power_dissipation.lookupOrDefault<scalar>("number_power_dissipation_Patches",3);//number of patch 
wordList POWERPatchNames = power_dissipation.lookup("power_dissipation_patches_names");//function name
label POWERPatchList[nPOWERPatch];//patch list
int inum;
for (inum=0; inum<nPOWERPatch; inum++)
{
    POWERPatchList[inum] = mesh.boundaryMesh().findPatchID(POWERPatchNames[inum]);
}

dictionary mixing_index = mesh.solutionDict().subDict("mixing_index");
int nMIPatch = mixing_index.lookupOrDefault<scalar>("number_mixing_index_patches",0);	
wordList MIPatchNames = mixing_index.lookup("mixing_index_patches_names");
label MIPatchList[nMIPatch];
for (inum=0; inum<nMIPatch; inum++)
{
    MIPatchList[inum] =
         mesh.boundaryMesh().findPatchID(MIPatchNames[inum]);
}

scalar Cv_value = 0.5;
volScalarField Cv
(
    IOobject
    (
        "Cv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
	),
    mesh,
    Cv_value
);