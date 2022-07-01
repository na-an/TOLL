loopnumUc = 0;
Ucresidual0 = 10.0;
Info<< "adjoint Disspower equation loop" << endl;
while  (Ucresidual0 > Ucresidualset && loopnumUc < MaxUcloop)
{
    #include "UcEqn.H"
    #include "pcEqn.H"
    loopnumUc++;
}
