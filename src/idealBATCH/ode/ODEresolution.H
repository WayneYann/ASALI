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
    ChangeDimensions( ode->NumberOfEquations(), &yFbzz, true);
    ChangeDimensions( ode->NumberOfEquations(), &x0bzz, true);
    ChangeDimensions( ode->NumberOfEquations(), &yF, true);

    FromOSToBzz(x0,x0bzz);

    BzzOdeStiffObject o;
    o(x0bzz,0.,ode);

    if ( input.getConstraints() == true)
    {
        BzzVector yMax(ode->NumberOfEquations());
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
        BzzVector yMin(ode->NumberOfEquations());
        {
            unsigned counter = 1;
            {
                for (unsigned int j=1;j<=NC;j++)
                    yMin[counter++] = 0.;
                yMin[counter++] = 0.;
                for (unsigned int j=1;j<=SURF_NP;j++)
                    yMin[counter++] = 0.;
                for (unsigned int j=1;j<=SURF_NC;j++)
                    yMin[counter++] = 0.;
                yMin[counter++] = 0.;
            }
        }

        o.SetMinimumConstraints(&yMin);
        o.SetMaximumConstraints(&yMax);
    }

    {
        BzzVector tolAbs(ode->NumberOfEquations());
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
        BzzVector tolRel(ode->NumberOfEquations());
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

    ode->start();
    yFbzz = o(input.getIntegrationTime());

    FromBzzToOS(yFbzz,yF);

    ode->end();
    #endif
}
else if ( input.getOdeSolver() == "Sundials" )
{
    #if ASALI_USE_SUNDIALS == 1
    ChangeDimensions( ode->NumberOfEquations(), &yF, true);

    OpenSMOKE::ODESystem_CVODE_Template *odeSystemObject;
    odeSystemObject = OpenSMOKE::ODESystem_CVODE_Template::GetInstance();
    odeSystemObject->SetOdeSystem(ode);

    OpenSMOKE::OpenSMOKE_CVODE_Sundials<OpenSMOKE::ODESystem_CVODE_Template>   o(odeSystemObject);
    o.SetDimensions(ode->NumberOfEquations());
    o.SetAbsoluteTolerance(input.getSpecieAbsTol());
    o.SetRelativeTolerance(input.getSpecieRelTol());
    o.SetMaximumNumberOfSteps(10000);
    o.SetAnalyticalJacobian(false);
    o.SetInitialValues(0., x0.GetHandle());

    ode->start();
    o.Solve(input.getIntegrationTime());
    o.Solution(yF.GetHandle());
    
    ode->end();
    #endif
}
else if ( input.getOdeSolver() == "OpenSMOKE" )
{
	Eigen::VectorXd x0eigen(ode->NumberOfEquations());
	Eigen::VectorXd yFeigen(ode->NumberOfEquations());

	memcpy(x0eigen.data(), x0.GetHandle(), x0.Size()*sizeof(double));

    ChangeDimensions( ode->NumberOfEquations(), &yF, true);

	typedef OdeSMOKE::KernelDense<OpenSMOKE::ODESystem_ODE_OpenSMOKE> denseOde;
	typedef OdeSMOKE::MethodGear<denseOde> methodGearOde;
	OdeSMOKE::MultiValueSolver<methodGearOde> o;

	o.SetOdeSystem(ode);
	o.SetInitialConditions(0., x0eigen);
    o.SetAbsoluteTolerances(input.getSpecieAbsTol());
    o.SetRelativeTolerances(input.getSpecieRelTol());

    if ( input.getConstraints() == true)
    {
        Eigen::VectorXd yMax(ode->NumberOfEquations());
        {
            unsigned counter = 0;
            {
                for (unsigned int j=0;j<NC;j++)
                    yMax(counter++) = 1.;
                yMax(counter++) = 5000.;
                for (unsigned int j=0;j<SURF_NP;j++)
                    yMax(counter++) = 5000.;
                for (unsigned int j=0;j<SURF_NC;j++)
                    yMax(counter++) = 1.;
                yMax(counter++) = 5000.;
            }
        }
        Eigen::VectorXd  yMin(ode->NumberOfEquations());
        {
            unsigned counter = 0;
            {
                for (unsigned int j=0;j<NC;j++)
                    yMin(counter++) = 0.;
                yMin(counter++) = 0.;
                for (unsigned int j=0;j<SURF_NP;j++)
                    yMin(counter++) = 0.;
                for (unsigned int j=0;j<SURF_NC;j++)
                    yMin(counter++) = 0.;
                yMin(counter++) = 0.;
            }
        }

        o.SetMinimumValues(yMin);
        o.SetMaximumValues(yMax);
    }

    ode->start();
    OdeSMOKE::OdeStatus status = o.Solve(input.getIntegrationTime());
    o.Solution(yFeigen);
    
    memcpy(yF.GetHandle(), yFeigen.data(), yFeigen.size()*sizeof(double));
    
    ode->end();
}
