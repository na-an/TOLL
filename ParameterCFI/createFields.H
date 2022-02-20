Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    fvc::flux(U)
);
label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, simple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

Info<< "\nReading g" << endl;
uniformDimensionedVectorField g
(
    IOobject
    (
        "g",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);                                   
// adjoint
Info<< "Reading field pa\n" << endl;
volScalarField pa
(
    IOobject
    (
        "pa",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Ua\n" << endl;
volVectorField Ua
(
    IOobject
    (
        "Ua",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

surfaceScalarField phia
(
    IOobject
    (
        "phia",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
	linearInterpolate(Ua) & mesh.Sf()
);

label paRefCell = 0;
scalar paRefValue = 0.0;
setRefCell
(
    pa,
    simple.dict(),
    paRefCell,
    paRefValue
);

mesh.setFluxRequired(pa.name());


Info<< "Reading field pc\n" << endl;
volScalarField pc
(
    IOobject
    (
        "pc",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Uc\n" << endl;
volVectorField Uc
(
    IOobject
    (
        "Uc",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

surfaceScalarField phic
(
    IOobject
    (
        "phic",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    linearInterpolate(Uc) & mesh.Sf()
);

label pcRefCell = 0;
scalar pcRefValue = 0.0;
setRefCell
(
    pc,
    simple.dict(),
    pcRefCell,
    pcRefValue
);

mesh.setFluxRequired(pc.name());

volScalarField alphac
(
    IOobject
    (
        "alphac",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField adalpha
(
    IOobject
    (
        "adalpha",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);

  
Info<< "Reading Material Properties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar changealpha("1", dimMass/(dimLength*dimTime*dimTime), 1.0);
dimensionedScalar changealphanan("1", dimMass/(dimLength*dimLength*dimLength*dimTime), 1.0);



dimensionedScalar Dab(transportProperties.lookup("Dab"));

dimensionedScalar rho1(transportProperties.lookup("rho1"));
dimensionedScalar rho2(transportProperties.lookup("rho2"));

dimensionedScalar nua(transportProperties.lookup("nua"));
dimensionedScalar lambdadan(transportProperties.lookup("lambdadan"));


dimensionedScalar alphasolidini(transportProperties.lookup("alphasolidini"));
dimensionedScalar alphaMaxdelta(transportProperties.lookup("alphaMaxdelta"));
dimensionedScalar alphasolidmax(transportProperties.lookup("alphasolidmax"));
dimensionedScalar alphaMax(transportProperties.lookup("alphaMax"));

dimensionedScalar alphasolidmin(transportProperties.lookup("alphasolidmin"));

dimensionedScalar Rmin(transportProperties.lookup("Rmin"));



scalar voluse(readScalar(transportProperties.lookup("voluse")));
scalar filter_Radius(readScalar(transportProperties.lookup("filter_Radius")));

scalar raa0(readScalar(transportProperties.lookup("raa0")));
scalar mma_init(readScalar(transportProperties.lookup("mma_init")));
scalar mma_dec(readScalar(transportProperties.lookup("mma_dec")));
scalar mma_inc(readScalar(transportProperties.lookup("mma_inc")));
scalar movlim(readScalar(transportProperties.lookup("movlim")));

scalar D_normalization(readScalar(transportProperties.lookup("D_normalization")));
scalar D0(readScalar(transportProperties.lookup("D0")));
scalar D1(readScalar(transportProperties.lookup("D1")));

scalar J0(readScalar(transportProperties.lookup("J0")));
scalar J1(readScalar(transportProperties.lookup("J1")));

scalar geo_dim(readScalar(transportProperties.lookup("geo_dim")));


scalar asqa(readScalar(transportProperties.lookup("asqa")));
scalar asqaini(readScalar(transportProperties.lookup("asqaini")));
scalar asqadelta(readScalar(transportProperties.lookup("asqadelta")));
scalar asqamax(readScalar(transportProperties.lookup("asqamax")));

volScalarField rho("rho", alphac*rho1 + (1-alphac)*rho2);

word zoneName1 = "zone_fluid"; 
label zoneID1 = mesh.cellZones().findZoneID(zoneName1);
labelList cells_fluid = mesh.cellZones()[zoneID1];

word zoneName2 = "zone_design"; 
label zoneID2 = mesh.cellZones().findZoneID(zoneName2);
labelList cells_design = mesh.cellZones()[zoneID2];


// design variable dvg


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
    scalar(voluse),
    zeroGradientFvPatchScalarField::typeName 
);

setCells(dvg, cells_fluid, 1);

//volScalarField dvgPDE(dvg);
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
    scalar(voluse),
    zeroGradientFvPatchScalarField::typeName 
);
//volScalarField dvgHea(dvg);
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
    scalar(voluse),
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

//Alldvg
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

volScalarField alphadao
(
    IOobject
    (
        "alphadao",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    -(alphasolidmax-alphasolidmin)*asqa*(asqa+1)/((asqa+dvgHea)*(asqa+dvgHea))
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
    (alphasolidmax-alphasolidmin)*asqa*(1-dvgHea)/(asqa+dvgHea)*(Uc & U),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField fsensm(fsens);
//volScalarField fsensPDE(fsens);
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
    (alphasolidmax-alphasolidmin)*asqa*(1-dvgHea)/(asqa+dvgHea)*(Uc & U),
    zeroGradientFvPatchScalarField::typeName
);
//volScalarField fsensHea(fsens);
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
    (alphasolidmax-alphasolidmin)*asqa*(1-dvgHea)/(asqa+dvgHea)*(Uc & U),
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
    (alphasolidmax-alphasolidmin)*asqa*(1-dvgHea)/(asqa+dvgHea)*(Ua & U),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField gsensm(gsens);
//volScalarField gsensPDE(gsens);
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
    (alphasolidmax-alphasolidmin)*asqa*(1-dvgHea)/(asqa+dvgHea)*(Ua & U),
    zeroGradientFvPatchScalarField::typeName
);
//volScalarField gsensHea(gsens);
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
    (alphasolidmax-alphasolidmin)*asqa*(1-dvgHea)/(asqa+dvgHea)*(Ua & U),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField UcU
(
    IOobject
    (
        "UcU",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    Uc & U
);

volScalarField UaU
(
    IOobject
    (
        "UaU",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    Ua & U
);


#include "createFvOptions.H"
//sensitivity and objFunction

dictionary optFunc = mesh.solutionDict().subDict("objectiveFunctionDict");
int nObjPatch = optFunc.lookupOrDefault<scalar>("numberObjectivePatches",0);	
wordList objPatchNames = optFunc.lookup("objectivePatchesNames");
label objPatchList[nObjPatch];
int iLoop;
for (iLoop=0; iLoop<nObjPatch; iLoop++)
{
    objPatchList[iLoop] =
         mesh.boundaryMesh().findPatchID(objPatchNames[iLoop]);
}

dictionary Fun1=mesh.solutionDict().subDict("power_dissipation");
int nObjPatchpow=Fun1.lookupOrDefault<scalar>("numberConstraintPatches",3);//number of patch 
wordList conPatchNames = Fun1.lookup("constraintPatchesNames");//function name
label conPatchList[nObjPatchpow];//patch list
int iLooppow;
for (iLooppow=0; iLooppow<nObjPatchpow; iLooppow++)
{
    conPatchList[iLooppow] = mesh.boundaryMesh().findPatchID(conPatchNames[iLooppow]);
}

Info<<"Initializing objective function calculation:"<<endl;
Info<<"Name of the patches for which the cost function will be calculated:"<<objPatchNames<<endl;

// Constants for flow uni at the outlet.
scalar alphalookup = optFunc.lookupOrDefault<scalar>("alphaduserdefnodim",scalar(0.5));

volScalarField alphad
(
    IOobject
    (
        "alphad",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
	),
    mesh,
    alphalookup
);

Info<< "Finish Reading field" << endl;
