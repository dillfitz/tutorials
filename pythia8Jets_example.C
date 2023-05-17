/// Jet hadronization studies with Pythia8 + FASTJET for B+ in jet correlation studies ///
/// Author(s): Dillon Fitzgerald, Nicole Kuchta, Julia Marchese, Jose Marco Arias (Markiito) ///

#include "TSystem.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include "fastjet/ClusterSequence.hh"
#include <iostream>
#include "TVector3.h"
#include "TLorentzVector.h"
#include <vector>

using namespace fastjet;
using namespace std;

class MyInfo: public 
PseudoJet::UserInfoBase {
  public:
    MyInfo(int id) : _pdg_id(id){}
    int pdg_id() const {return _pdg_id;}
    int _pdg_id;
    };

int verbosity = 0;
void pythia8Jets_example(Int_t nev  = 1000 , Int_t ndeb = 1)
{
	double ptHatMin = 17.0;
	//double ptHatMax = 50.0; //
	TString ptHatMin_str = "PhaseSpace:pTHatMin = ";
	//TString ptHatMax_str = "PhaseSpace:pTHatMax = "; //
	ptHatMin_str += ptHatMin;
	//ptHatMax_str += ptHatMax;//

	TString outFileName = "pythia8Jets_example.root";
	//outFileName += ptHatMin;
	//outFileName += "_";
	//outFileName += ptHatMax;
	//outFileName += ".root";

	TFile *outFile = new TFile(outFileName,"RECREATE");
	// Load libraries //
  gSystem->Load("libEG");
  gSystem->Load("libEGPythia8");

	// Setup histograms // 
	TH1D *ptHat_hist = new TH1D("ptHat", ";p_{T} [GeV/c];", 100, 0.0, 100.0);

	// Set up jet tree //
  int  nConstituents, eventNum;
  double jPt, jPx, jPy, jPz, jE; 

  TTree* jettree;
  jettree = new TTree("jets","A tree with jet info");
  jettree->Branch("nConstituents", &nConstituents, "nConstituents/I"); 
  jettree->Branch("jPt" , &jPt,  "jPt/D" );
  jettree->Branch("jPx" , &jPx, "jPx/D");
  jettree->Branch("jPy" , &jPy, "jPy/D");
  jettree->Branch("jPz" , &jPz, "jPz/D");
  jettree->Branch("jE" , &jE, "jE/D");
  jettree->Branch("eventNum" , &eventNum,  "eventNum/I");

  //Set up B-hadron Tree
  int bPDG;
  double bPx, bPy, bPz, bE;
  double JPsiPx, JPsiPy, JPsiPz, JPsiE; 
  double KaonPx, KaonPy, KaonPz, KaonE; 
  double Mu1Px, Mu1Py, Mu1Pz, Mu1E;
  double Mu2Px, Mu2Py, Mu2Pz, Mu2E;


  // Set up jet constituent tree //
  double cPt, cJt, cR, cz;
  int cPid;

  TTree* contree;
  contree = new TTree("cons", "A tree with jet's constituent information");
  contree->Branch("cPt" , &cPt,  "cPt/D" ); 
  contree->Branch("cJt", &cJt, "cJt/D" );      
  contree->Branch("cR", &cR, "cR/D" );         
  contree->Branch( "cz" , &cz, "cz/D" );
  contree->Branch( "cPid" , &cPid, "cPid/I" );  

  // choose a jet definition //
  double R = 0.5;
  JetDefinition jet_def(antikt_algorithm, R);
  vector<PseudoJet> parts;

	// Array of particles //
  TClonesArray* particles = new TClonesArray("TParticle", 1000000);  
	// Create pythia8 object //
  TPythia8* pythia8 = new TPythia8();

  #if PYTHIA_VERSION_INTEGER == 8235
    // Pythia 8.235 is known to cause crashes: //
    printf("ABORTING PYTHIA8 TUTORIAL!\n");
   	printf("The version of Pythia you use is known to case crashes due to memory errors.\n");
   	printf("They have been reported to the authors; the Pythia versions 8.1... are known to work.\n");
   	return;
  #endif

  // Configure
  pythia8->ReadString("HardQCD:all = on"); 
  pythia8->ReadString("Random:setSeed = on");
  
  // use a reproducible seed: always the same results for the tutorial.
  pythia8->ReadString("Random:seed = 42"); //make this random per event
  
  // Here is the pT hat cut... //
  pythia8->ReadString(ptHatMin_str);
  //pythia8->ReadString(ptHatMax_str);  //
  

	// Initialize
	// RHIC energy //
  //pythia8->Initialize(2212 /* p */, 2212 /* p */, 130 /* GeV */);

	// LHC energy //
  pythia8->Initialize(2212 /* p */, 2212 /* p */, 13000. /* GeV */);                        //13.7 TeV Max achieved thie July!
  
	// EIC beams and energies //
  //pythia8->Initialize(11 /* e ^-*/, 2212 /* p */, 18. /* GeV */, 275. /* GeV */);


  // Event loop
  for (Int_t iev = 0; iev < nev; iev++) 
	{
    eventNum = iev;
    pythia8->GenerateEvent();
    //cout << "pythia hard process... " << pythia8->Pythia8()->info.name() << " " << pythia8->Pythia8()->info.code() << endl;
    if (verbosity > 0)
      pythia8->EventListing();

    pythia8->ImportParticles(particles,"All");
    Int_t np = particles->GetEntriesFast();
    parts.clear();
       
		// Particle loop (1)
    for (Int_t ip = 0; ip < np; ip++) 
    {

      TParticle* part = (TParticle*) particles->At(ip);
      Int_t ist = part->GetStatusCode();
	    Int_t pdg = part->GetPdgCode();
      
	    // LHCb acceptance cut... //
	    //if (part->Eta() < 2. || part->Eta() > 5.)  {continue;}

      // Positive codes are final particles.
      if (ist <= 0) continue;
      
      Float_t eta = part->Eta();    //pseudoripidity of a particular particle 
      Float_t pt  = part->Pt();
      // Only cluster particles into jets that are above 200 MeV
      if (pt < 0.2) continue;
      
      // Do not cluster neutrinos into jets... (i.e. skip them!) //
      if (abs(pdg) == 12 || abs(pdg) == 14 || abs(pdg) == 16) continue;
      
      // append particle 4 vectors into a vector using fastjet::PseudoJet object for jet clustering
      parts.push_back(PseudoJet(part->Px(), part->Py(), part->Pz(), part->Energy() ));
      parts.back().set_user_info(new MyInfo(pdg));
      
    } // end particle loop (1)

    // Cluster the jets //
    ClusterSequence cs(parts, jet_def);
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

    // Jet loop //
    for (unsigned i = 0; i < jets.size(); i++)
    {
      if (jets[i].pt() < (ptHatMin - 5.0)) {continue;}
      // LHCb acceptance cut for jets... Is it redundant to have both? I think we would experimentally... //
	    //if (jets[i].eta() < 2 + R || jets[i].eta() > 4.5 - R) {continue;}
        
      PseudoJet jet = jets[i];
	    jPt = jet.pt(); 
	    vector<PseudoJet> constituents = jet.constituents();
	    nConstituents = constituents.size();

      TVector3 jet_3vec( jet.px(), jet.py(), jet.pz() );
      
      jPx = jet.px();
      jPy = jet.py();
      jPz = jet.pz();
      jE  = jet.E();
	  
      //ptHat->Fill(ptHat_hist)

      // Constituent loop //
      for (unsigned j = 0; j < constituents.size(); j++)
      {
        PseudoJet con = constituents[j];
        //cout << con.user_info<MyInfo>().pdg_id() << endl;
        
        cPid = con.user_info<MyInfo>().pdg_id();

		    TVector3 cons_3vec( con.px() , con.py(), con.pz()); 
        cz = ( jet_3vec * cons_3vec ) / ( jet.modp2() );
	      cJt = (jet_3vec.Cross(cons_3vec).Mag())/(jet_3vec.Mag());
        cR = sqrt( pow( jet.phi() - con.phi() , 2) + pow(jet.pseudorapidity() - con.pseudorapidity(), 2));

        cPt = con.pt();       
        contree->Fill();
        
      }
      // end constituent loop //
      jettree->Fill();

    }
    // end jet loop //
  }
  // end event loop //
  pythia8->PrintStatistics();
  outFile->Write();
  outFile->Close();
}
