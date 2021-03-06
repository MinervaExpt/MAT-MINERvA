#ifndef PARTICLERESPONSECONSTANTS_CXX
#define PARTICLERESPONSECONSTANTS_CXX 1
#include "utilities/ParticleResponseDefaults.h"

using namespace PlotUtils;
double PartRespDefaults::GetDefaultPartRespFracUnc( std::string particle )
{
  if( particle.find("proton")!=std::string::npos || particle.compare("p")==0 )       return PartRespDefaults::protonResponse;
  if( particle.find("low_proton")!=std::string::npos)   return PartRespDefaults::protonResponseLowEkin ; 
  if( particle.find("mid_proton")!=std::string::npos )   return PartRespDefaults::protonResponseMidEkin ; 
  if( particle.find("high_proton")!=std::string::npos )  return PartRespDefaults::protonResponseHighEkin; 
  if( particle.find("low_neutron")!=std::string::npos || particle.compare("lown")==0)  return PartRespDefaults::neutronResponseLowEkin ; 
  if( particle.find("mid_neutron")!=std::string::npos || particle.compare("midn")==0)  return PartRespDefaults::neutronResponseMidEkin ; 
  if( particle.find("high_neutron")!=std::string::npos || particle.compare("highn")==0) return PartRespDefaults::neutronResponseHighEkin; 
  if( particle.find("meson")!=std::string::npos )        return PartRespDefaults::mesonResponse          ;
  if( particle.find("em")!=std::string::npos )           return PartRespDefaults::electromagneticResponse;
  if( particle.find("xtalk")!=std::string::npos )        return PartRespDefaults::xtalkResponse  ;
  if( particle.find("other")!=std::string::npos )        return PartRespDefaults::otherParticleResponse  ;
  return 0;
};
#endif
