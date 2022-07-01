loopnumU = 0;
Uresidual0 = 10.0;
Info<< "primal equation loop" << endl;
while (Uresidual0 > Uresidualset && loopnumU < MaxUloop)
{
    #include "CEqn.H"
    #include "UEqn.H"
    #include "pEqn.H"
    loopnumU++;
}