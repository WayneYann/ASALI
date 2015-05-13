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
    bvp->resize(N);
    bvp->setGrid(z);
    bvp->setDiscretizationScheme(discretizationScheme);

    BzzVector yFbzz;
    BzzVector x0bzz;
    ChangeDimensions( bvp->NumberOfEquations(), &yFbzz, true);
    ChangeDimensions( bvp->NumberOfEquations(), &x0bzz, true);
    ChangeDimensions( bvp->NumberOfEquations(), &yF, true);

    FromOSToBzz(x0,x0bzz);

    BzzVectorInt iAO;
    bvp->AlgebraicEquations(iAO);

    BzzDaeSparseObject o;
    o(x0bzz,0.,iAO,bvp,bvp->BlockDimension());

    if ( input.getConstraints() == true)
    {
        BzzVector yMax(NE);
        {
            unsigned counter = 1;
            for (unsigned int k=1;k<=N;k++)
            {
                for (unsigned int j=1;j<=NC;j++)
                    yMax[counter++] = 1.;
                for (unsigned int j=1;j<=NC;j++)
                    yMax[counter++] = 1.;
                for (unsigned int j=1;j<=SURF_NC;j++)
                    yMax[counter++] = 1.;
                yMax[counter++] = 5000.;
                yMax[counter++] = 5000.;
            }
        }
        BzzVector yMin(NE);
        {
            unsigned counter = 1;
            for (unsigned int k=1;k<=N;k++)
            {
                for (unsigned int j=1;j<=NC;j++)
                    yMin[counter++] = 0.;
                for (unsigned int j=1;j<=NC;j++)
                    yMin[counter++] = 0.;
                for (unsigned int j=1;j<=SURF_NC;j++)
                    yMin[counter++] = 0.;
                yMin[counter++] = 0.;
                yMin[counter++] = 0.;
            }
        }
        o.SetMinimumConstraints(&yMin);
        o.SetMaximumConstraints(&yMax);
    }

    {
        BzzVector tolAbs(NE);
        {
            unsigned counter = 1;
            for (unsigned int k=1;k<=N;k++)
            {
                for (unsigned int j=1;j<=NC;j++)
                    tolAbs[counter++] = input.getSpecieAbsTol();
                for (unsigned int j=1;j<=NC;j++)
                    tolAbs[counter++] = input.getSpecieAbsTol();
                for (unsigned int j=1;j<=SURF_NC;j++)
                    tolAbs[counter++] = input.getSpecieAbsTol();
                tolAbs[counter++] = input.getTemperatureAbsTol();
                tolAbs[counter++] = input.getTemperatureAbsTol();
            }
        }
        BzzVector tolRel(NE);
        {
            unsigned counter = 1;
            for (unsigned int k=1;k<=N;k++)
            {
                for (unsigned int j=1;j<=NC;j++)
                    tolRel[counter++] = input.getSpecieRelTol();
                for (unsigned int j=1;j<=NC;j++)
                    tolRel[counter++] = input.getSpecieRelTol();
                for (unsigned int j=1;j<=SURF_NC;j++)
                    tolRel[counter++] = input.getSpecieRelTol();
                tolRel[counter++] = input.getTemperatureRelTol();
                tolRel[counter++] = input.getTemperatureRelTol();
            }
        }

        o.SetTolAbs(tolAbs);
        o.SetTolRel(tolRel);
    }

    bvp->start();
    yFbzz = o(input.getIntegrationTime());

    FromBzzToOS(yFbzz,yF);

    bvp->end();
    #endif
}
else if ( input.getDaeSolver() == "Sundials" )
{
    #if ASALI_USE_SUNDIALS == 1
    bvp->resize(N);
    bvp->setGrid(z);

    OpenSMOKE::OpenSMOKEVectorDouble residuals(bvp->NumberOfEquations());
    OpenSMOKE::OpenSMOKEVectorBool   algebraic(bvp->NumberOfEquations());

    bvp->AlgebraicEquations(algebraic);
    ChangeDimensions( bvp->NumberOfEquations(), &yF, true);

    OpenSMOKE::DAESystem_IDA_Template *daeSystemObject;
    daeSystemObject = OpenSMOKE::DAESystem_IDA_Template::GetInstance(bvp->NumberOfEquations(),algebraic.GetHandle());
    daeSystemObject->SetDaeSystem(bvp);

    OpenSMOKE::OpenSMOKE_IDA_Sundials<OpenSMOKE::DAESystem_IDA_Template>   o(daeSystemObject);
    o.SetDimensions(bvp->NumberOfEquations());
    o.SetTridiagonalBlock(bvp->BlockDimension());
    o.SetAbsoluteTolerance(input.getSpecieAbsTol());
    o.SetRelativeTolerance(input.getSpecieRelTol());
    o.SetMaximumNumberOfSteps(100000);
    o.SetAnalyticalJacobian(false);
    o.SetAlgebraicEquations(algebraic.GetHandle());
    o.SetInitialValues(0., x0.GetHandle(), residuals.GetHandle());

    bvp->start();
    o.Solve(input.getIntegrationTime());
    o.Solution(yF.GetHandle(), residuals.GetHandle());

    bvp->end();
    #endif
}
