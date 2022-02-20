
resi3=0;
Ucresidualsave=1.0;
while  (Ucresidualsave>1e-3&&resi3<200)
{
    #include "UcEqn.H"
    #include "pcEqn.H"    
    resi3++;
}
