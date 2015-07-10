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
    ChangeDimensions( ic->NumberOfEquations(), &yFbzz, true);
    ChangeDimensions( ic->NumberOfEquations(), &x0bzz, true);
    ChangeDimensions( ic->NumberOfEquations(), &yF, true);

    FromOSToBzz(x0,x0bzz);

    BzzOdeStiffObject o;
    o(x0bzz,0.,ic);

    if ( input.getConstraints() == true)
    {
        BzzVector yMax(ic->NumberOfEquations());
        {
            unsigned counter = 1;
            {
                for (unsigned int j=1;j<=NC;j++)
                    yMax[counter++] = 1.;
                yMax[counter++] = 5000.;
                for (unsigned int j=1;j<=SURF_NP;j++)
                    yMax[counter++] = 5000.;
                for (unsigned int j=1;j<=SURF_NC;j++)
                    yMax[counter++] = 1.;
                yMax[counter++] = 5000.;
            }
        }
        BzzVector yMin(ic->NumberOfEquations());
        {
            unsigned counter = 1;
            {
                for (unsigned int j=1;j<=NC;j++)
                    yMax[counter++] = 0.;
                yMax[counter++] = 0.;
                for (unsigned int j=1;j<=SURF_NP;j++)
                    yMax[counter++] = 0.;
                for (unsigned int j=1;j<=SURF_NC;j++)
                    yMax[counter++] = 1.;
                yMax[counter++] = 0.;
            }
        }

        o.SetMinimumConstraints(&yMin);
        o.SetMaximumConstraints(&yMax);
    }

    {
        BzzVector tolAbs(ic->NumberOfEquations());
        {
            unsigned counter = 1;
            {
                for (unsigned int j=1;j<=NC;j++)
                    tolAbs[counter++] = input.getSpecieAbsTol();;
                tolAbs[counter++] = input.getSpecieAbsTol();;
                for (unsigned int j=1;j<=SURF_NP;j++)
                    tolAbs[counter++] = input.getSpecieAbsTol();;
                for (unsigned int j=1;j<=SURF_NC;j++)
                    tolAbs[counter++] = input.getSpecieAbsTol();
                tolAbs[counter++] = input.getTemperatureAbsTol();
            }
        }
        BzzVector tolRel(ic->NumberOfEquations());
        {
            unsigned counter = 1;
            {
                for (unsigned int j=1;j<=NC;j++)
                    tolRel[counter++] = input.getSpecieRelTol();;
                tolRel[counter++] = input.getSpecieRelTol();;
                for (unsigned int j=1;j<=SURF_NP;j++)
                    tolRel[counter++] = input.getSpecieRelTol();;
                for (unsigned int j=1;j<=SURF_NC;j++)
                    tolRel[counter++] = input.getSpecieRelTol();
                tolRel[counter++] = input.getTemperatureRelTol();
            }
        }

       o.SetTolAbs(tolAbs);
       o.SetTolRel(tolRel);
    }

    ic->start();
    yFbzz = o(0.5);

    FromBzzToOS(yFbzz,yF);

	for(unsigned int i=1;i<=ic->NumberOfEquations();i++)
	{
		x0[i] = yF[i];
		yF[i] = 0.;
	}

    ic->end();
    #endif
}
else if ( input.getOdeSolver() == "Sundials" )
{
    #if ASALI_USE_SUNDIALS == 1
    ChangeDimensions( ic->NumberOfEquations(), &yF, true);

    OpenSMOKE::ODESystem_CVODE_Template *icSystemObject;
    icSystemObject = OpenSMOKE::ODESystem_CVODE_Template::GetInstance();
    icSystemObject->SetOdeSystem(ic);

    OpenSMOKE::OpenSMOKE_CVODE_Sundials<OpenSMOKE::ODESystem_CVODE_Template>   o(icSystemObject);
    o.SetDimensions(ic->NumberOfEquations());
    o.SetAbsoluteTolerance(input.getSpecieAbsTol());
    o.SetRelativeTolerance(input.getSpecieRelTol());
    o.SetMaximumNumberOfSteps(10000);
    o.SetAnalyticalJacobian(false);
    o.SetInitialValues(0., x0.GetHandle());

    ic->start();
    o.Solve(0.5);
    o.Solution(yF.GetHandle());

	for(unsigned int i=1;i<=ic->NumberOfEquations();i++)
	{
		x0[i] = yF[i];
		yF[i] = 0.;
	}

    ic->end();
    #endif
}
