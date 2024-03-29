#ifndef RESPONSESYSTEMATICS_CXX
#define RESPONSESYSTEMATICS_CXX

#include "utilities/ParticleResponseDefaults.h"
#include "universes/ResponseSystematics.h"
#include <iostream>

using namespace PlotUtils;

//Including a boolean to give user ability to include/exclude neutron systematics...Do not turn on unless you know what you are doing
//Any trace of Neutron Systematics  will/should vanish from here once it is confirmed that
//it is implemented properly on Hadron Reweight Systematics.

namespace PlotUtils
{
  template <class T>
  std::vector<T*>GetResponseSystematics(typename T::config_t chain, int sigma,
				      bool NEUTRON=false, bool use_new_part_resp = false, bool PROTON = false){
    std::cout << "Recoil Systematics created with neutron " << NEUTRON << std::endl;
    std::cout << "Recoil Systematics created with proton "  << PROTON << std::endl;
                                                                       
    std::vector<T*>ret;                                                
                                                                       
    std::vector<std::string> response_systematics;                     
      if(PROTON)
      {
        response_systematics.push_back("low_proton");
        response_systematics.push_back("mid_proton");
        response_systematics.push_back("high_proton");
      }
      else response_systematics.push_back("proton");
      response_systematics.push_back("meson");
      response_systematics.push_back("em");
      //response_systematics.push_back("xtalk"); //DO NOT USE
      response_systematics.push_back("other");
      if(NEUTRON){
        response_systematics.push_back("low_neutron");
        response_systematics.push_back("mid_neutron");
        response_systematics.push_back("high_neutron");
      }
 
    for(std::vector<std::string>::const_iterator syst = response_systematics.begin(); 
            syst != response_systematics.end(); ++syst){
      ret.push_back( new PlotUtils::ResponseUniverse<T>(chain, sigma, *syst, use_new_part_resp) );
    }
    
    return ret;
  }

  //Overloaded version used to pass a specific recoil name for the response branches
  template <class T>
  std::vector<T*>GetResponseSystematics(typename T::config_t chain, int sigma, std::string name_tag, 
					bool NEUTRON=false, bool use_new_part_resp = false, bool PROTON = false){
    std::cout << "Recoil Systematics created with neutron " << NEUTRON << std::endl;
    std::cout << "Recoil Systematics created with proton "  << PROTON << std::endl;
                                                                       
    std::vector<T*>ret;                                                
                                                                       
    std::vector<std::string> response_systematics;                     
      if(PROTON)
      {
        response_systematics.push_back("low_proton");
        response_systematics.push_back("mid_proton");
        response_systematics.push_back("high_proton");
      }
      else response_systematics.push_back("proton");
      response_systematics.push_back("meson");
      response_systematics.push_back("em");
      //response_systematics.push_back("xtalk"); //DO NOT USE
      response_systematics.push_back("other");
      if(NEUTRON){
        response_systematics.push_back("low_neutron");
        response_systematics.push_back("mid_neutron");
        response_systematics.push_back("high_neutron");
      }
 
    for(std::vector<std::string>::const_iterator syst = response_systematics.begin(); 
            syst != response_systematics.end(); ++syst){
      ret.push_back( new PlotUtils::ResponseUniverse<T>(chain, sigma, name_tag, *syst, use_new_part_resp) );
    }
    
    return ret;
  }

  //Overloaded version used to pass a specific recoil name for the response branches && vector of subdetector names
  /*template <class T>
  std::vector<T*>GetResponseSystematics(typename T::config_t chain, int sigma, std::string name_tag, std::vector<std::string> subdetectors, 
					bool NEUTRON=false, bool use_new_part_resp = false, bool PROTON = false){
    std::cout << "Recoil Systematics created with neutron " << NEUTRON << std::endl;
    std::cout << "Recoil Systematics created with proton "  << PROTON << std::endl;
                                                                       
    std::vector<T*>ret;

    std::vector<std::string> response_systematics;                     
    for(auto subdet : subdetectors){                                             
                                                                       
      if(PROTON)
        {
          response_systematics.push_back(Form("%s_low_proton", subdet.c_str()));
          response_systematics.push_back(Form("%s_mid_proton", subdet.c_str()));
          response_systematics.push_back(Form("%s_high_proton", subdet.c_str()));
        }
        else response_systematics.push_back(Form("%s_proton", subdet.c_str()));
        response_systematics.push_back(Form("%s_meson", subdet.c_str()));
        response_systematics.push_back(Form("%s_em", subdet.c_str()));
        //response_systematics.push_back("xtalk"); //DO NOT USE
        response_systematics.push_back(Form("%s_other", subdet.c_str()));
        if(NEUTRON){
          response_systematics.push_back(Form("%s_low_neutron", subdet.c_str()));
          response_systematics.push_back(Form("%s_mid_neutron", subdet.c_str()));
          response_systematics.push_back(Form("%s_high_neutron", subdet.c_str()));
        }
    
    }
 
    for(std::vector<std::string>::const_iterator syst = response_systematics.begin(); 
            syst != response_systematics.end(); ++syst){
      ret.push_back( new PlotUtils::ResponseUniverse<T>(chain, sigma, name_tag, *syst, use_new_part_resp) );
    }
    
    return ret;
  }*/

  template <class T>
  std::map< std::string, std::vector<T*> > GetResponseSystematicsMap(typename T::config_t chain, bool NEUTRON=false, bool use_new_part_resp = false, bool PROTON = false) {
    std::map< std::string, std::vector<T*> > ret;

    std::cout << "Lecacy error recoil branches" << std::endl;

    std::vector<double> sigmas;
    sigmas.push_back(-1.);
    sigmas.push_back(+1.);

    std::vector<std::string> response_systematics;
      if(PROTON)
      {
        response_systematics.push_back("low_proton");
        response_systematics.push_back("mid_proton");
        response_systematics.push_back("high_proton");
      }
      else response_systematics.push_back("proton");
      response_systematics.push_back("meson");
      response_systematics.push_back("em");
      //response_systematics.push_back("xtalk"); //DO NOT USE
      response_systematics.push_back("other");
      if(NEUTRON){
        response_systematics.push_back("low_neutron");
        response_systematics.push_back("mid_neutron");
        response_systematics.push_back("high_neutron");
      }

    for(std::vector<double>::const_iterator sigma = sigmas.begin(); 
            sigma != sigmas.end(); ++sigma) {
      for(std::vector<std::string>::const_iterator syst = response_systematics.begin(); 
              syst != response_systematics.end(); ++syst){
        ret[*syst].push_back(new PlotUtils::ResponseUniverse<T>(chain, *sigma, *syst, use_new_part_resp));
      }
    }

    return ret;
  }

  //Overloaded version used to pass a specific recoil name for the response branches
  template <class T>
  std::map< std::string, std::vector<T*> > GetResponseSystematicsMap(typename T::config_t chain, std::string name_tag, bool NEUTRON=false, bool use_new_part_resp = false, bool PROTON = false) {
    std::map< std::string, std::vector<T*> > ret;

    std::cout << "Lecacy error recoil branches" << std::endl;

    std::vector<double> sigmas;
    sigmas.push_back(-1.);
    sigmas.push_back(+1.);

    std::vector<std::string> response_systematics;
      if(PROTON)
      {
        response_systematics.push_back("low_proton");
        response_systematics.push_back("mid_proton");
        response_systematics.push_back("high_proton");
      }
      else response_systematics.push_back("proton");
      response_systematics.push_back("meson");
      response_systematics.push_back("em");
      //response_systematics.push_back("xtalk"); //DO NOT USE
      response_systematics.push_back("other");
      if(NEUTRON){
        response_systematics.push_back("low_neutron");
        response_systematics.push_back("mid_neutron");
        response_systematics.push_back("high_neutron");
      }

    for(std::vector<double>::const_iterator sigma = sigmas.begin(); 
            sigma != sigmas.end(); ++sigma) {
      for(std::vector<std::string>::const_iterator syst = response_systematics.begin(); 
              syst != response_systematics.end(); ++syst){
        ret[*syst].push_back(new PlotUtils::ResponseUniverse<T>(chain, *sigma, name_tag, *syst, use_new_part_resp));
      }
    }

    return ret;
  }

  template <class T>
  std::map< std::string, std::vector<T*> > GetResponseSystematicsMap(typename T::config_t chain, bool ID, bool OD, std::string name_tag, bool NEUTRON=false, bool use_new_part_resp = false, bool PROTON = false, bool nucl = false, bool tracker = false, bool ecal = false, bool hcal = false ) {
    std::map< std::string, std::vector<T*> > ret;

    std::cout << "New error recoil branches" << std::endl;

    std::vector<double> sigmas;
    sigmas.push_back(-1.);
    sigmas.push_back(+1.);

    std::vector<std::string> response_systematics;
      if(PROTON)
      {
        response_systematics.push_back("low_proton");
        response_systematics.push_back("mid_proton");
        response_systematics.push_back("high_proton");
      }
      else response_systematics.push_back("proton");
      response_systematics.push_back("meson");
      response_systematics.push_back("em");
      //response_systematics.push_back("xtalk"); //DO NOT USE
      response_systematics.push_back("other");
      if(NEUTRON){
        response_systematics.push_back("low_neutron");
        response_systematics.push_back("mid_neutron");
        response_systematics.push_back("high_neutron");
      }

    for(std::vector<double>::const_iterator sigma = sigmas.begin(); 
            sigma != sigmas.end(); ++sigma) {
      for(std::vector<std::string>::const_iterator syst = response_systematics.begin(); 
              syst != response_systematics.end(); ++syst){
        ret[*syst].push_back(new PlotUtils::ResponseUniverse<T>(chain, *sigma,  ID, OD, name_tag, *syst, use_new_part_resp, nucl, tracker, ecal, hcal, true));
      }
    }

    return ret;
  }
}

  //Overloaded version used to pass a specific recoil name for the response branches && vector of subdetector names
  /*template <class T>
  std::map< std::string, std::vector<T*> > GetResponseSystematicsMap(typename T::config_t chain, std::string name_tag, std::vector<std::string> subdetectors, bool NEUTRON=false, bool use_new_part_resp = false, bool PROTON = false) {
    std::map< std::string, std::vector<T*> > ret;

    std::vector<double> sigmas;
    sigmas.push_back(-1.);
    sigmas.push_back(+1.);

    std::vector<std::string> response_systematics;                     
    for(auto subdet : subdetectors){                                             


        if(PROTON)
        {
          response_systematics.push_back(Form("%s_low_proton", subdet.c_str()));
          response_systematics.push_back(Form("%s_mid_proton", subdet.c_str()));
          response_systematics.push_back(Form("%s_high_proton", subdet.c_str()));
        }
        else response_systematics.push_back(Form("%s_proton", subdet.c_str()));
        response_systematics.push_back(Form("%s_meson", subdet.c_str()));
        response_systematics.push_back(Form("%s_em", subdet.c_str()));
        //response_systematics.push_back("xtalk"); //DO NOT USE
        response_systematics.push_back(Form("%s_other", subdet.c_str()));
        if(NEUTRON){
          response_systematics.push_back(Form("%s_low_neutron", subdet.c_str()));
          response_systematics.push_back(Form("%s_mid_neutron", subdet.c_str()));
          response_systematics.push_back(Form("%s_high_neutron", subdet.c_str()));
        } 
    }

    for(std::vector<double>::const_iterator sigma = sigmas.begin(); 
            sigma != sigmas.end(); ++sigma) {
      for(std::vector<std::string>::const_iterator syst = response_systematics.begin(); 
              syst != response_systematics.end(); ++syst){
        ret[*syst].push_back(new PlotUtils::ResponseUniverse<T>(chain, *sigma, name_tag, *syst, use_new_part_resp, true));
      }
    }

    return ret;
  }
}*/
  
//=================================================================================
// Response 
//=================================================================================
// Constructor
template <class T>
ResponseUniverse<T>::ResponseUniverse(typename T::config_t chw, 
				      double nsigma, std::string response_name, bool use_new_part_resp)
  : T(chw, nsigma), 
    m_name(response_name), m_branch_name("part_response_recoil_" + response_name + "_id_err"),
    m_use_new_part_resp(use_new_part_resp), 
    m_particle_passive_name("part_response_recoil_passive_" + response_name + "_id"),
    m_container_particle_passive_name("part_response_container_recoil_passive_"+response_name+"_id")
{
  m_frac_unc = PartRespDefaults::GetDefaultPartRespFracUnc( response_name ); 
}

//Constructor to allow for a specific recoil name being added
template <class T>
ResponseUniverse<T>::ResponseUniverse(typename T::config_t chw, 
				      double nsigma, std::string name_tag, std::string response_name, bool use_new_part_resp)
  : T(chw, nsigma), 
    m_name(response_name), m_branch_name( name_tag.length() == 0 ? "part_response_recoil_" + response_name + "_id_err": "part_response_recoil_" + name_tag + "_" + response_name + "_id_err"),
    m_use_new_part_resp(use_new_part_resp), 
    m_particle_passive_name("part_response_recoil_passive_" + response_name + "_id"),
    m_container_particle_passive_name("part_response_container_recoil_passive_"+response_name+"_id")
{
  m_frac_unc = PartRespDefaults::GetDefaultPartRespFracUnc( response_name ); 
}

//Constructor to allow for a specific recoil name being added
/*template <class T>
ResponseUniverse<T>::ResponseUniverse(typename T::config_t chw, 
				      double nsigma, std::string name_tag, std::string response_name, bool use_new_part_resp, bool p4)
  : T(chw, nsigma), 
    m_name(response_name), m_branch_name( name_tag.length() == 0 ? "part_response_recoil_" + response_name + "_err": "part_response_recoil_" + name_tag + "_" + response_name + "_err"),
    m_use_new_part_resp(use_new_part_resp), 
    m_particle_passive_name("part_response_recoil_passive_" + response_name + "_id"),
    m_container_particle_passive_name("part_response_container_recoil_passive_"+response_name+"_id")
{
  m_frac_unc = PartRespDefaults::GetDefaultPartRespFracUnc( response_name ); 
}*/

template <class T>
ResponseUniverse<T>::ResponseUniverse(typename T::config_t chw, 
				      double nsigma, bool ID, bool OD, std::string name_tag, std::string response_name, bool use_new_part_resp, bool nucl, bool tracker, bool ecal, bool hcal, bool p4)
  : T(chw, nsigma), 
    m_name(response_name),
    m_nametag(name_tag),
    m_ID(ID), m_OD(OD), m_nucl(ID == true ? false : nucl), m_tracker(ID == true ? false : tracker), m_ecal(ID == true ? false : ecal), m_hcal(ID == true ? false : hcal),
    m_p4(p4),
    m_use_new_part_resp(use_new_part_resp), 
    m_particle_passive_name("part_response_recoil_passive_" + response_name + "_id"),
    m_container_particle_passive_name("part_response_container_recoil_passive_"+response_name+"_id")
{
  m_frac_unc = PartRespDefaults::GetDefaultPartRespFracUnc( response_name ); 
}


template <class T>
double ResponseUniverse<T>::GetRecoilShift() const {
  double shift_val = 0;
  if( m_use_new_part_resp )
  {
    //ID only
    //Total recoil
    //  part_response_total_recoil_passive_id/od
    //  part_response_container_total_recoil_passive_id/od
    //Particle recoil
    //  part_response_recoil_passive_<response_name>_id/od
    //  part_response_container_recoil_passive_<response_name>_id/od
    double total_passive_recoil = T::GetDouble("part_response_total_recoil_passive_id");
    double part_passive_recoil  = T::GetDouble(m_particle_passive_name.c_str()); 
     
    for( std::vector<int>::const_iterator idx = T::m_non_cal_idx.begin(); idx != T::m_non_cal_idx.end(); ++idx  )
    {
      total_passive_recoil -= T::GetVecElem("part_response_container_total_recoil_passive_id", *idx);
      part_passive_recoil  -= T::GetVecElem( m_container_particle_passive_name.c_str(), *idx);
    }
      
    shift_val = T::m_nsigma*( (part_passive_recoil/total_passive_recoil)*m_frac_unc )*T::GetCalRecoilEnergy();
  }
  else
  { //std::cout<< "Using err branches" << std::endl;
    //std::cout<< m_branch_name.c_str() << std::endl;
    //double recoil_E_shift = T::GetDouble(m_branch_name.c_str());
    double recoil_E_shift = 0;
    if(m_p4){
      if(m_ID){
        std::string nucl = "part_response_recoil_" + m_nametag + "_nucl_" + m_name +"_err";
        std::string tracker = "part_response_recoil_" + m_nametag + "_tracker_" + m_name +"_err";
        std::string ecal = "part_response_recoil_" + m_nametag + "_ecal_" + m_name +"_err";
        std::string hcal = "part_response_recoil_" + m_nametag + "_hcal_" + m_name +"_err";

        recoil_E_shift += T::GetDouble(nucl.c_str()) + T::GetDouble(tracker.c_str()) + T::GetDouble(ecal.c_str()) + T::GetDouble(hcal.c_str());
        //std::cout<< "id" << std::endl;
      }
      if(m_OD){
        std::string od = "part_response_recoil_" + m_nametag + "_" + m_name + "_od_err";
        recoil_E_shift += T::GetDouble(od.c_str());
        //std::cout<< "od" << std::endl;
      }
      if (m_nucl){
        std::string nucl = "part_response_recoil_" + m_nametag + "_nucl_" + m_name +"_err";
        recoil_E_shift += T::GetDouble(nucl.c_str());
        //std::cout<< "nucl"<<  std::endl;
      }
      if (m_tracker){
        std::string tracker = "part_response_recoil_" + m_nametag + "_tracker_" + m_name +"_err";
        recoil_E_shift += T::GetDouble(tracker.c_str());
        //std::cout<< "tracker"<<  std::endl;
      }
      if (m_ecal){
        std::string ecal = "part_response_recoil_" + m_nametag + "_ecal_" + m_name +"_err";
        recoil_E_shift += T::GetDouble(ecal.c_str());
        //std::cout<< "ecal" << std::endl;
      }
      if (m_hcal){
        std::string hcal = "part_response_recoil_" + m_nametag + "_hcal_" + m_name +"_err";
        recoil_E_shift += T::GetDouble(hcal.c_str());
        //std::cout<< "hcal" << std::endl;
      }
    }
    else{
      std::cout << "Legacy error recoil branches" << std::endl;
      recoil_E_shift += T::GetDouble(m_branch_name.c_str());
    }

    // Guard against implementation of default branch value (which is -1)
    // That is, if the branch is filled with its default, treat the shift as 0
    recoil_E_shift = (recoil_E_shift>0)? recoil_E_shift : 0;
    shift_val = T::m_nsigma*recoil_E_shift;
  }
  return shift_val;
}

//template <class T>
//double ResponseUniverse<T>::GetRecoilEnergy() const {
//  double shift_val = GetRecoilShift();
//  return shift_val+T::GetRecoilEnergy();
//}

template <class T>
double ResponseUniverse<T>::GetCalRecoilEnergy() const {
  double shift_val = GetRecoilShift();
  return shift_val+T::GetCalRecoilEnergy();
}

template <class T>
std::string ResponseUniverse<T>::ShortName() const {
  return "response_" + m_name;
};

template <class T>
std::string ResponseUniverse<T>::LatexName() const {
  return "Response " + m_name;
};


#endif // RESPONSESYSTEMATICS_CXX
