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
    ode->resize(N);
    ode->setGrid(z);
    ode->setDiscretizationScheme(discretizationScheme);

    BzzVector yFbzz;
    BzzVector x0bzz;
    ChangeDimensions( ode->NumberOfEquations(), &yFbzz, true);
    ChangeDimensions( ode->NumberOfEquations(), &x0bzz, true);
    ChangeDimensions( ode->NumberOfEquations(), &yF, true);

    FromOSToBzz(x0,x0bzz);

    BzzOdeSparseStiffObject o;
    o(x0bzz,0.,ode,ode->BlockDimension());

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
                tolRel[counter++] = input.getTemperatureRelTol();
            }
        }

        o.SetTolAbs(tolAbs);
        o.SetTolRel(tolRel);
    }

    ode->start();
    yFbzz = o(1e10);

    FromBzzToOS(yFbzz,yF);

    for (unsigned int i=1;i<=ode->NumberOfEquations();i++)
    {
        x0[i] = yF[i];
        yF[i] = 0.;
    }
    
    ode->end();
    #endif
}
else if ( input.getOdeSolver() == "Sundials" )
{
    #if ASALI_USE_SUNDIALS == 1
    ode->resize(N);
    ode->setGrid(z);
    ode->setDiscretizationScheme(discretizationScheme);

    ChangeDimensions( ode->NumberOfEquations(), &yF, true);

    OpenSMOKE::ODESystem_CVODE_Template *odeSystemObject;
    odeSystemObject = OpenSMOKE::ODESystem_CVODE_Template::GetInstance();
    odeSystemObject->SetOdeSystem(ode);

    OpenSMOKE::OpenSMOKE_CVODE_Sundials<OpenSMOKE::ODESystem_CVODE_Template>   o(odeSystemObject);
    o.SetDimensions(ode->NumberOfEquations());
    o.SetTridiagonalBlock(ode->BlockDimension());
    o.SetAbsoluteTolerance(input.getSpecieAbsTol());
    o.SetRelativeTolerance(input.getSpecieRelTol());
    o.SetMaximumNumberOfSteps(10000);
    o.SetAnalyticalJacobian(false);
    o.SetInitialValues(0., x0.GetHandle());

    ode->start();
    double tf = 1e10;
    o.Solve(tf);
    o.Solution(yF.GetHandle());
    
    for (unsigned int i=1;i<=ode->NumberOfEquations();i++)
    {
        x0[i] = yF[i];
        yF[i] = 0.;
    }
    
    ode->end();
    #endif
}
