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

namespace ASALI
{
class ODESystem
#if ASALI_USE_BZZ == 1
 : public BzzOdeSystemObject
#endif
{

public:

	ODESystem(OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&          thermodynamicsMap, 
			  OpenSMOKE::KineticsMap_CHEMKIN<double>&                kineticsMap,
			  OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN<double>&  thermodynamicsSurfaceMap, 
			  OpenSMOKE::KineticsMap_Surface_CHEMKIN_Lumped<double>& kineticsSurfaceMap,
			  OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>&     transportMap);

	#include "vector.h"

	void setHomogeneusReactions(const bool flag)   { homogeneusReactions_ = flag; }

	void setHeterogeneusReactions(const bool flag) { heterogeneusReactions_ = flag; }


	void setReactorGeometry(const double alfa, const double av, const double G, const double Dh);

	void setPressure(const double p);

	void setTemperature (const double T);

	void setInert(const std::string inert);

	unsigned int NumberOfEquations()	const {return NE_;};

	void start();
	void end();

	unsigned int Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy);

	unsigned int Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);

	#if ASALI_USE_BZZ == 1
	virtual void GetSystemFunctions(BzzVector &y, double t, BzzVector &dy);
	virtual void ObjectBzzPrint(void);
	virtual ~ODESystem(){};
	#endif


private:

	double MWbulk_;
	double MWwall_;
	double cTotBulk_;
	double cTotWall_;
	double rhoBulk_;
	double rhoWall_;
	double etaMix_;
	double p_;
	double T_;
	double L_;
	double Dh_;
	double alfa_;
	double G_;
	double SD_;
	double av_;
	double t_;

	unsigned int NC_;
	unsigned int SURF_NP_;
	unsigned int SURF_NC_;
	unsigned int NE_;

	unsigned int iteration;

	std::string inert_;

	bool homogeneusReactions_ ;
	bool heterogeneusReactions_;

	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&			thermodynamicsMap_;				//!< thermodynamic map
	OpenSMOKE::KineticsMap_CHEMKIN<double>&					kineticsMap_;					//!< kinetic map
	OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN<double>&	thermodynamicsSurfaceMap_;		//!< thermodynamic map
	OpenSMOKE::KineticsMap_Surface_CHEMKIN<double>&			kineticsSurfaceMap_;			//!< kinetic map
	OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>& 		transportMap_;					//!< transport map

	OpenSMOKE::OpenSMOKEVectorDouble omegaBulk_;
	OpenSMOKE::OpenSMOKEVectorDouble omegaWall_;

	OpenSMOKE::OpenSMOKEVectorDouble teta_;

	OpenSMOKE::OpenSMOKEVectorDouble RfromGas_;
	OpenSMOKE::OpenSMOKEVectorDouble RfromSurface_;
	OpenSMOKE::OpenSMOKEVectorDouble Rsurface_;

	OpenSMOKE::OpenSMOKEVectorDouble diffG_;
	OpenSMOKE::OpenSMOKEVectorDouble Kmat_;
	OpenSMOKE::OpenSMOKEVectorDouble jSolid_;

	OpenSMOKE::OpenSMOKEVectorDouble xBulk_;
	OpenSMOKE::OpenSMOKEVectorDouble xWall_;
	OpenSMOKE::OpenSMOKEVectorDouble cBulk_;
	OpenSMOKE::OpenSMOKEVectorDouble cWall_;
	OpenSMOKE::OpenSMOKEVectorDouble RsurfacePhases_;

	OpenSMOKE::OpenSMOKEVectorDouble dyOS_;
	OpenSMOKE::OpenSMOKEVectorDouble  yOS_;

	void MassTransferCoefficient(const double z);

};


ODESystem::ODESystem(	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&          thermodynamicsMap, 
						OpenSMOKE::KineticsMap_CHEMKIN<double>&                kineticsMap,
						OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN<double>&  thermodynamicsSurfaceMap, 
						OpenSMOKE::KineticsMap_Surface_CHEMKIN_Lumped<double>& kineticsSurfaceMap,
						OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>&     transportMap):
	thermodynamicsMap_(thermodynamicsMap), 
	kineticsMap_(kineticsMap),
	thermodynamicsSurfaceMap_(thermodynamicsSurfaceMap), 
	kineticsSurfaceMap_(kineticsSurfaceMap),
	transportMap_(transportMap)
	{
		MWbulk_    = 0.;
		MWwall_    = 0.;
		cTotBulk_  = 0.;
		cTotWall_  = 0.;
		rhoBulk_   = 0.;
		rhoWall_   = 0.;
		etaMix_    = 0.;
		p_         = 0.;
		T_         = 0.;
		Dh_        = 0.;
		alfa_      = 0.;
		G_         = 0.;
		SD_        = 0.;
		av_        = 0.;
		t_         = 0.;

		NC_        = 0;
		SURF_NP_   = 0;
		SURF_NC_   = 0;
		NE_        = 0;

		iteration  = 0;

		NC_      = thermodynamicsMap_.NumberOfSpecies();
		SURF_NC_ = thermodynamicsSurfaceMap_.number_of_site_species();
		SURF_NP_ = thermodynamicsSurfaceMap_.number_of_site_phases(0);
		NE_      = NC_ + NC_ + SURF_NC_ ;

		ChangeDimensions(NC_, &omegaBulk_, true);
		ChangeDimensions(NC_, &omegaWall_, true);
		ChangeDimensions(SURF_NC_, &teta_, true);

		ChangeDimensions(NC_, &RfromGas_,      true);
		ChangeDimensions(NC_, &RfromSurface_,  true);
		ChangeDimensions(SURF_NC_, &Rsurface_, true);

		ChangeDimensions(NC_, &xBulk_,  true);
		ChangeDimensions(NC_, &xWall_,  true);
		ChangeDimensions(NC_, &cBulk_,  true);
		ChangeDimensions(NC_, &cWall_,  true);
		ChangeDimensions(NC_, &diffG_,  true);
		ChangeDimensions(NC_, &Kmat_,   true);
		ChangeDimensions(NC_, &jSolid_, true);

		ChangeDimensions(NE_, &dyOS_, true);
		ChangeDimensions(NE_, &yOS_,  true);
	}

void ODESystem::setInert(const std::string inert)
{
	inert_		= inert;
}

void ODESystem::setReactorGeometry(const double alfa, const double av, const double G, const double Dh)
{
	alfa_			= alfa;
	G_				= G;
	av_				= av;
	Dh_             = Dh;
}

void ODESystem::setPressure(const double p)
{
	p_ = p;
}

void ODESystem::setTemperature (const double T)
{
	T_ = T;
}

void ODESystem::MassTransferCoefficient(const double z)
{
	double Re = G_*Dh_/(etaMix_);
	OpenSMOKE::OpenSMOKEVectorDouble Sc(NC_);
	OpenSMOKE::OpenSMOKEVectorDouble Sh(NC_);
	double zNew;
	double zStar;
	for (unsigned int i=1;i<=NC_;i++)
	{
		Sc[i] = etaMix_/(rhoBulk_*diffG_[i]);
		zNew  = 1e-06;
		zStar = zNew/(Dh_*Re*Sc[i]);
		zStar = fabs(zStar);
		Sh[i] = 5.21 + 6.874*pow((1000.*zStar),-0.35)*exp(-71.2*zStar);
		Kmat_[i] = Sh[i]*diffG_[i]/Dh_;
	}
}

#if ASALI_USE_BZZ == 1
void ODESystem::GetSystemFunctions(BzzVector &y, double t, BzzVector &dy)
{
	FromBzzToOS(y,yOS_);
	#include "ODEequations.H"
	ObjectBzzPrint();
	FromOSToBzz(dyOS_,dy);
}

void ODESystem::ObjectBzzPrint(void)
{
	#include "printIntegration.H"
}
#endif

unsigned int ODESystem::Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
{
	ChangeDimensions(NE_, &yOS_,  true);
	ChangeDimensions(NE_, &dy,    true);
	ChangeDimensions(NE_, &dyOS_, true);

	for (unsigned int i=1;i<=NE_;i++)
		yOS_[i] = y[i];

	#include "ODEequations.H"

	for (unsigned int i=1;i<=NE_;i++)
		dy[i] = dyOS_[i];
	
	return 0;
}

unsigned int ODESystem::Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
{
	t_ = t;
	#include "printIntegration.H"
	return 0;
}

void ODESystem::start()
{
	std::cout << "\n######################################" << std::endl;
	std::cout << "# ODE system:                 START  #" << std::endl;
	std::cout << "######################################\n" << std::endl;
}

void ODESystem::end()
{
	std::cout << "\n######################################" << std::endl;
	std::cout << "# ODE system:                 END    #" << std::endl;
	std::cout << "######################################\n" << std::endl;
}

}
