
resi2=0;
Uaresidualsave=1.0;
while (Uaresidualsave>1e-3&&resi2<200)
{
    #include "adjointalphaE.H"
    #include "UaEqn.H"
    #include "paEqn.H"
    resi2++;
}
