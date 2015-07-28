/*##############################################################################################
#                                                                                              #
#     #############       #############       #############       ####                ####     #
#    #             #     #             #     #             #     #    #              #    #    #
#    #    #####    #     #    #########      #    #####    #     #    #              #    #    #
#    #    #   #    #     #    #              #    #   #    #     #    #              #    #    #
#    #    #####    #     #    #              #    #####    #     #    #              #    #    #
#    #             #     #    #########      #             #     #    #              #    #    #
#    #             #     #             #     #             #     #    #              #    #    #
#    #    #####    #      #########    #     #    #####    #     #    #              #    #    #
#    #    #   #    #              #    #     #    #   #    #     #    #              #    #    #
#    #    #   #    #      #########    #     #    #   #    #     #    #########      #    #    #
#    #    #   #    #     #             #     #    #   #    #     #             #     #    #    #
#     ####     ####       #############       ####     ####       #############       ####     #
#                                                                                              #
#   Department of Energy                                                                       #
#   Politecnico di Milano                                                                      #
#   Author: Stefano Rebughini <stefano.rebughini@polimi.it>                                    #
#                                                                                              #
################################################################################################
#                                                                                              #
#   License                                                                                    #
#                                                                                              #
#   This file is part of ASALI.                                                                #
#                                                                                              #
#   ASALI is free software: you can redistribute it and/or modify                              #
#   it under the terms of the GNU General Public License as published by                       #
#   the Free Software Foundation, either version 3 of the License, or                          #
#   (at your option) any later version.                                                        #
#                                                                                              #
#   ASALI is distributed in the hope that it will be useful,                                   #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of                             #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                              #
#   GNU General Public License for more details.                                               #
#                                                                                              #
#   You should have received a copy of the GNU General Public License                          #
#   along with ASALI. If not, see <http://www.gnu.org/licenses/>.                              #
#                                                                                              #
##############################################################################################*/

if ( input.getOdeSolver() == "BzzMath")
{
    #if ASALI_USE_BZZ == 1

    BzzVector yFbzz;
    BzzVector x0bzz;
    ChangeDimensions( dae->NumberOfEquations(), &yFbzz, true);
    ChangeDimensions( dae->NumberOfEquations(), &x0bzz, true);
    ChangeDimensions( dae->NumberOfEquations(), &yF, true);

    for (unsigned int i=1;i<=dae->NumberOfEquations();i++)
    {
        x0bzz[i] = yFode[i];
    }

    BzzVectorInt iAO;
    dae->AlgebraicEquations(iAO);

    BzzDaeObject o;
    o(x0bzz,0.,iAO,dae);

    o.SetTolAbs(1.e-12);
    o.SetTolRel(1.e-08);

    yFbzz = o(input.getReactorLength());

    FromBzzToOS(yFbzz,yF);

    #endif
}
else if ( input.getDaeSolver() == "Sundials" )
{
    #if ASALI_USE_SUNDIALS == 1

    OpenSMOKE::OpenSMOKEVectorDouble y0(dae->NumberOfEquations());
    ChangeDimensions(dae->NumberOfEquations(), &yF, true);

    for (unsigned int i=1;i<=NC;i++)
    {
        y0[i] = x0[i];
    }

    for (unsigned int i=1;i<=NC;i++)
    {
        y0[i+NC] = x0[i];
    }
 
    y0[dae->NumberOfEquations()] = input.getFeedPressure();

    OpenSMOKE::OpenSMOKEVectorBool   algebraic(dae->NumberOfEquations());
    dae->AlgebraicEquations(algebraic);

    OpenSMOKE::DAESystem_IDA_Template *daeSystemObject;
    daeSystemObject = OpenSMOKE::DAESystem_IDA_Template::GetInstance(dae->NumberOfEquations(),algebraic.GetHandle());
    daeSystemObject->SetDaeSystem(dae);

    OpenSMOKE::OpenSMOKE_IDA_Sundials<OpenSMOKE::DAESystem_IDA_Template>   o(daeSystemObject);
    o.SetDimensions(dae->NumberOfEquations());
    o.SetAbsoluteTolerance(input.getSpecieAbsTol());
    o.SetRelativeTolerance(input.getSpecieRelTol());
    o.SetMaximumNumberOfSteps(100000);
    o.SetAnalyticalJacobian(false);
    o.SetAlgebraicEquations(algebraic.GetHandle());
    o.SetInitialValues(0., x0.GetHandle(), residuals.GetHandle());

    // Solution
    o.Solve(input.getReactorLength());
    o.Solution(yF.GetHandle());
    #endif
}
