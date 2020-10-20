#ifndef TRUTHRIVETTOOLS_DIHIGGSRIVET_CC
#define TRUTHRIVETTOOLS_DIHIGGSRIVET_CC

// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Projections/FastJets.hh"

// Definition of the StatusCode and Category enums
#include "DiHiggsRivet.h"

namespace Rivet {
  
  class DiHiggsRivet : public Analysis {
  public:
    // Constructor
    DiHiggsRivet()
      : Analysis("DiHiggsRivet"),
    m_ProdMode(HH::UNKNOWN) {}

  public:

    /// @name Utility methods
    /// Methods to identify the Higgs boson and
    /// associated vector boson and to build jets
    /// @{

    /// follow a "propagating" particle and return its last instance
    Particle getLastInstance(Particle ptcl) {
      if ( ptcl.genParticle()->end_vertex() ) {
	if ( !hasChild(ptcl.genParticle(),ptcl.pdgId()) ) return ptcl;
	else return getLastInstance(ptcl.children()[0]);
      }
      return ptcl;
    }
    
    /// @brief Whether particle p originate from any of the ptcls
    bool originateFrom(const Particle& p, const Particles& ptcls ) {
      const GenVertex* prodVtx = p.genParticle()->production_vertex();
      if (prodVtx == nullptr) return false;
      // for each ancestor, check if it matches any of the input particles
      for (auto ancestor:particles(prodVtx, HepMC::ancestors)){ 
	for ( auto part:ptcls ) 
	  if ( ancestor==part.genParticle() ) return true;
      }
      // if we get here, no ancetor matched any input particle
      return false; 
    }
    
    /// @brief Whether particle p originates from p2
    bool originateFrom(const Particle& p, const Particle& p2 ) { 
      Particles ptcls = {p2}; return originateFrom(p,ptcls);
    }
    
    /// @brief Checks whether the input particle has a child with a given PDGID 
    bool hasChild(const GenParticle *ptcl, int pdgID) {
      for (auto child:Particle(*ptcl).children())
        if (child.pdgId()==pdgID) return true;
      return false;
    }
    
    /// @brief Checks whether the input particle has a parent with a given PDGID 
    bool hasParent(const GenParticle *ptcl, int pdgID) {
      for (auto parent:particles(ptcl->production_vertex(),HepMC::parents))
        if (parent->pdg_id()==pdgID) return true;
      return false;
    }

    /// @brief Return true if particle decays to quarks
    bool quarkDecay(const Particle &p) {
      for (auto child:p.children())
        if (PID::isQuark(child.pdgId())) return true;
      return false;
    }

    /// @brief Return true if particle decays to charged leptons.
    bool ChLeptonDecay(const Particle &p) {
      for (auto child:p.children())
        if (PID::isChLepton(child.pdgId())) return true;
      return false;
    }

    /// @brief Return bin index of x given the provided bin edges. 0=first bin, -1=underflow bin.
    int getBin(double x, std::vector<double> bins) {
      if (bins.size()==0||x<bins[0]) return -1; // should not happen!
      for (size_t i=1;i<bins.size();++i)
        if (x<bins[i]) return i-1;
      return bins.size()-1;
    }
    
    /// @brief VBF topology selection
    /// 0 = fail loose selection: m_jj > 350 GeV
    /// 1 pass loose, but fail additional cut pT(Hjj)<25. 2 pass pT(Hjj)>25 selection
    /// 3 pass tight (m_jj>700 GeV), but fail additional cut pT(Hjj)<25. 4 pass pT(Hjj)>25 selection
    int vbfTopology(const Jets &jets, const Particle &higgs) {
      if (jets.size()<2) return 0;
      const FourMomentum &j1=jets[0].momentum(), &j2=jets[1].momentum();
      double mjj = (j1+j2).mass();
      if(mjj>350 && mjj<=700) return (j1+j2+higgs.momentum()).pt()<25 ? 1 : 2;
      else if(mjj>700) return (j1+j2+higgs.momentum()).pt()<25 ? 3 : 4;
      else return 0;
    }

    /// @brief Whether the Higgs is produced in association with a vector boson (VH)
    //bool isVH(HH::ProdMode p) { return p==HH::WH || p==HH::QQ2ZH || p==HH::GG2ZH; }
    
    /// @name Default Rivet analysis methods and steering methods
    /// @{

    /// @brief Sets the Higgs production mode
    void setProdMode( HH::ProdMode prodMode ){ m_ProdMode = prodMode; }

    /// @brief default Rivet Analysis::init method
    /// Booking of histograms, initializing Rivet projection
    /// Extracts Higgs production mode from shell variable if not set manually using setProdMode
    void init() {
      printf("==============================================================\n");
      printf("========     HH Production Initialization     =========\n");
      printf("==============================================================\n");
      // check that the production mode has been set
      // if running in standalone Rivet the production mode is set through an env variable
      if (m_ProdMode==HH::UNKNOWN) {
	char *pm_env = getenv("PRODMODE");
	string pm(pm_env==nullptr?"":pm_env);
        if      ( pm == "GGF"   ) m_ProdMode = HH::GGF; 
        else if ( pm == "VBF"   ) m_ProdMode = HH::VBF;
        else {
          MSG_WARNING("No PRODMODE shell variable found. Needed when running Rivet stand-alone.");
        }
      }



      // Projections for final state particles
      const FinalState FS; 
      addProjection(FS,"FS");
	
      // initialize the histograms with for each of the stages
      initializeHistos();
            
      m_sumw = 0.0;
      m_event = 0;
      printf("==============================================================\n");
      printf("========                   prod mode %d              =========\n",m_ProdMode);
      printf("========          Sucessful Initialization           =========\n");
      printf("==============================================================\n");
    }

   // Selecting DiHiggs system
    DiHiggsSys selectHH(const Event& event, const HH::ProdMode prodMode ) {

	if (m_ProdMode==HH::UNKNOWN) m_ProdMode = prodMode;
	DiHiggsSys HHSys;
	 
	HHSys.prodMode = prodMode;
	if(prodMode == HH::UNKNOWN){ HHSys.errorCode = HH::UNDEFINED; MSG_WARNING(" Error Production Mode not Defined " ); return HHSys;}
	
	// find the hard scatter vertex and Higgs 
    
        GenVertex* HSVTX = nullptr;
	int NH=0;
 	std::vector<const HepMC::GenParticle* > HiggsV; 
        for (auto ptcl : Rivet::particles(event.genEvent()) ) {
	
		if( not PID::isHiggs(ptcl->pdg_id()) ) continue;
		if(!hasChild(ptcl, PID::HIGGS)) // select only final state higgs 
		{
			HiggsV.push_back(ptcl);
			NH++; 
		}

        }
	
	if (NH != 2 or HiggsV.size() != 2) { HHSys.errorCode = HH::DIHIGGS_IDENTIFICATION; MSG_WARNING(" Error More/Less Higgs are founded, Found only Higgs " << HiggsV.size()); return HHSys; }

	if (HSVTX == nullptr and HiggsV[0]->production_vertex()) {
		// check if both Higgs has the same production vertex
		//if(HiggsV[0]->production_vertex() != HiggsV[1]->production_vertex()) { MSG_WARNING(" Error Higgs doesn't have the same production mode "); return HHSys;}
			
		// set hard scatter vertex to Higgs production vertex
		HSVTX = HiggsV[0]->production_vertex();
	}

	// set the two Higgs to HHSys, before order higgs by pT
	
	HHSys.H1 = Particle((HiggsV[0]->momentum().perp() > HiggsV[1]->momentum().perp()) ? HiggsV[0] : HiggsV[1]); 
	HHSys.H2 = Particle((HiggsV[0]->momentum().perp() < HiggsV[1]->momentum().perp()) ? HiggsV[0] : HiggsV[1]); 

	// building jets
	// final state particles - stable
	const ParticleVector FS = applyProjection<FinalState>(event, "FS").particles();
	Particles hadrons;
	FourMomentum sum(0,0,0,0), sumH1(0,0,0,0), sumH2(0,0,0,0);

	
	for ( auto p : FS ) {
		
		sum += p.momentum();
		
		// compute leading Higgs decay 4p
		if(originateFrom(p, HHSys.H1)) { sumH1 += p.momentum(); continue;} 		
		
		// sub-leading 
		if(originateFrom(p, HHSys.H2)) { sumH2 += p.momentum(); continue;} 
		
		// but default we shouldn't have lepton final state from V decays. 
		hadrons += p;
	}	
    
	HHSys.P4DECAY_H1 = sumH1;
	HHSys.P4DECAY_H2 = sumH2;

	FinalState FS_tmp; 
	FastJets jets(FS_tmp, FastJets::ANTIKT, 0.4);
	jets.calc(hadrons);

	HHSys.jets25 = jets.jetsByPt( Cuts::pT > 25.0 );
	HHSys.jets30 = jets.jetsByPt( Cuts::pT > 30.0 );
	
	if( sum.pt() > 0.1 ) { HHSys.errorCode = HH::MOMENTUM_CONSERVATION ; MSG_WARNING(" No momentum conservation up to 0.1 MeV : sum %f  " << sum.pt() ); return HHSys; }
	
	return HHSys;
    }
    // Perform the per-event analysis
    void analyze(const Event& event) {

	// get the classification
	DiHiggsSys HHSys = selectHH(event, m_ProdMode);

	const double weight = event.weight();
	m_sumw += weight;
	m_event++;
	if(m_event%10==0) std::cout<<" Number of events "<<m_event<<std::endl; 

	FourMomentum hh = HHSys.H1.momentum() + HHSys.H2.momentum();
	
	HH_m->fill(hh.mass(), weight);
	HH_pT->fill(  hh.pt(), weight); 
	
 }
   

    void finalize() {
      std::cout << "=========================================================" << std::endl;
      std::cout << " SumOfWeights  " << Analysis::sumOfWeights() << std::endl;
      std::cout << "Cross Section " << Analysis::crossSection() << std::endl;
      std::cout << " NumberOfEvents " << Analysis::numEvents() << std::endl;
      double sf = m_sumw>0?1.0/m_sumw:1.0;
  
      for( auto hist : {HH_m, HH_pT}) {
      	Analysis::normalize(hist, Analysis::crossSection());
      }
    }
    
    /*
     *  initialize histograms
     */

    void initializeHistos(){

	HH_m = bookHisto1D("HH_m",150,250,1000);
	HH_pT = bookHisto1D("HH_pT",200,0,1000);
      
    }
    /// @}

    /*
     *    initialize private members used in the classification procedure
     */
    
  private:
    double m_sumw;
    int m_event;
    HH::ProdMode m_ProdMode;
    std::map<HH::ErrorCode,size_t> m_errorCount;
   
    // histos 
    Histo1DPtr HH_m;
    Histo1DPtr HH_pT;
  };
  DECLARE_RIVET_PLUGIN(DiHiggsRivet);
  //#endif

}

#endif
    
