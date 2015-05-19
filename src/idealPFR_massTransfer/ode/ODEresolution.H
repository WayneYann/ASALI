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

	for (unsigned int i=1;i<=NC;i++)
	{
		x0bzz[i] = x0[i];
	}
	
	x0bzz[ode->NumberOfEquations()] = input.getFeedPressure();

	BzzOdeStiffObject o;
	o(x0bzz,0.,ode);

	o.SetTolAbs(1.e-12);
	o.SetTolRel(1.e-08);

	yFbzz = o(input.getReactorLength());

	FromBzzToOS(yFbzz,yF);

	#endif
}
else if ( input.getOdeSolver() == "Sundials" )
{
	#if ASALI_USE_SUNDIALS == 1

	OpenSMOKE::OpenSMOKEVectorDouble y0(ode->NumberOfEquations());
	ChangeDimensions(ode->NumberOfEquations(), &yF, true);

	for (unsigned int i=1;i<=NC;i++)
	{
		y0[i] = x0[i];
	}
	
	y0[ode->NumberOfEquations()] = input.getFeedPressure();

	OpenSMOKE::ODESystem_CVODE_Template *odeSystemObject;
	odeSystemObject = OpenSMOKE::ODESystem_CVODE_Template::GetInstance();
	odeSystemObject->SetOdeSystem(ode);

	OpenSMOKE::OpenSMOKE_CVODE_Sundials<OpenSMOKE::ODESystem_CVODE_Template>   o(odeSystemObject);
	o.SetDimensions(ode->NumberOfEquations());
	o.SetAbsoluteTolerance(1e-12);
	o.SetRelativeTolerance(1e-08);
	o.SetMaximumNumberOfSteps(20000);
	o.SetAnalyticalJacobian(false);
	o.SetInitialValues(0, y0.GetHandle());

	// Solution
	o.Solve(input.getReactorLength());
	o.Solution(yF.GetHandle());
	#endif
}
