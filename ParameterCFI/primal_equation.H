
resi1=0;
Uresidualsave=1.0;

while (Uresidualsave>1e-3&&resi1 <200)
{
    #include "alphaE.H"
    #include "UEqn.H"
    #include "pEqn.H"

    resi1++;
}
