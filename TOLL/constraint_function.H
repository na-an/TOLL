loopnumUa = 0;
Uaresidual0 = 10.0;
Info<< "adjoint MIX equation loop" << endl;
while (Uaresidual0 > Uaresidualset && loopnumUa < MaxUaloop)
{    
    
    #include "CaEqn.H"
    #include "UaEqn.H"
    #include "paEqn.H"

    loopnumUa++;
}
