// -*- C++ -*-
#include "Rivet/Analysis.hh"
//#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/HeavyHadrons.hh"
//! c++ headers
#include <iostream>     
#include <algorithm>    
#include <vector>

namespace Rivet {

  //! Analysis to get non-photonic electron (i.e. from bottom or charm hadron decays) spectra in 200 GeV p+p collisions 
  class STAR_2012_I1111571: public Analysis {
  public:    
    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(STAR_2012_I1111571);

    /// Book histograms and initialise projections before the run
    void init()
    {
      
      // Initialise and register projections
      //const ChargedFinalState cfs(Cuts::abseta < 1);
      const HeavyHadrons hh(Cuts::abseta < 1);
      declare(hh, "hh");
      // Book histograms
      book(_h_c_cbar , 1, 1, 1);
      
    }
    /// Perform the per-event analysis
    void analyze(const Event& event)
    {
      
      const FinalState& part = applyProjection<FinalState>(event, "hh");

      for (const Particle& charm : part.particles()) 
      {

	      if(!PID::hasCharm(charm.pid())) continue;   

      	const double pT = charm.pT() / GeV;
      	_h_c_cbar->fill(pT, (1./(2.0 * M_PI * 2.0))/pT); //2*pi*pT*deta - dpT should be done automatically
	      
      }
    } // end event ana   
    /// Normalise histograms etc., after the run
    void finalize() {
      //scale(_h_c_cbar,1./sumOfWeights());
      scale(_h_c_cbar,crossSection()/sumOfWeights()/millibarn);
      MSG_DEBUG("crossSection()     = " << crossSection());
      MSG_DEBUG("sumOfWeights()     = " << sumOfWeights());
    }
    
  private:
    Histo1DPtr _h_c_cbar;
  };

  DECLARE_RIVET_PLUGIN(STAR_2012_I1111571);

}
