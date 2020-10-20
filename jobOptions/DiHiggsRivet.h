#ifndef TRUTHRIVETTOOLS_DIHIGGSRIVET_H
#define TRUTHRIVETTOOLS_DIHIGGSRIVET_H

/// DiHiggs namespace
namespace HH {

  /// Error code: whether the classification was successful or failed
  enum ErrorCode {
    UNDEFINED=-99,
    SUCCESS = 0,               ///< successful classification
    PRODMODE_DEFINED = 1,      ///< production mode not defined
    MOMENTUM_CONSERVATION = 2, ///< failed momentum conservation
    DIHIGGS_IDENTIFICATION = 3,  ///< failed to identify DiHiggs Sys
  };

  /// Higgs production modes, corresponding to input sample
  enum ProdMode {
    UNKNOWN = 0,
    GGF = 1, VBF = 2
  };
  

#ifdef ROOT_TLorentzVector

    typedef TLorentzVector TLV;
    typedef std::vector<TLV> TLVs;
    
    template <class vec4>
      TLV MakeTLV(vec4 const p) { return TLV(p.px(),p.py(),p.pz(),p.E()); }
    
    template <class Vvec4>
      inline TLVs MakeTLVs(Vvec4 const &rivet_jets){ 
      TLVs jets; for ( auto jet:rivet_jets ) jets.push_back(MakeTLV(jet)); 
      return jets; 
    }
    
    // Structure holding information about the current event:
    // Four-momenta and event classification according to the
    // Higgs Template Cross Section
    struct DiHiggsSys {
      // Higgs production mode
      HH::ProdMode prodMode;
      // The Higgs boson
      TLV H1, H2;
      // The Higgs boson decay products
      TLV P4DECAY_H1, P4DECAT_H2;
  
      // Jets are built ignoring Higgs decay products and leptons from V decays
      // jets with pT > 25 GeV and 30 GeV
      TLVs jets25, jets30;

      HH::ErrorCode errorCode;
    };

#endif

} // namespace STXS


#ifdef RIVET_Particle_HH
//#ifdef HIGGSTRUTHCLASSIFIER_HIGGSTRUTHCLASSIFIER_CC
//#include "Rivet/Particle.hh"
namespace Rivet {

  /// @struct HiggsClassification
  /// @brief Structure holding information about the current event:
  ///        Four-momenta and event classification according to the
  ///        Higgs Template Cross Section
  struct DiHiggsSys {
    /// Higgs production mode
    HH::ProdMode prodMode;
    /// The Higgs boson
    Rivet::Particle H1, H2;
    /// Vector boson produced in asscoiation with the Higgs
    Rivet::FourMomentum P4DECAY_H1, P4DECAY_H2;
    /// Jets built ignoring Higgs decay products and leptons from V decays, pT thresholds at 25 GeV and 30 GeV
    Rivet::Jets jets25, jets30;
    /// Error code: Whether classification was succesful or some error occured
    HH::ErrorCode errorCode;
  };
} // namespace Rivet
#endif
#endif
