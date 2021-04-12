/*
Simple macro showing how to match generated jets with reconstructed jets by minimizing delta R.

root -l Matching.C'("sims/NC_DIS_EIC.root")'
*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif

//------------------------------------------------------------------------------

void Matching(const char *inputFile)
{
  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchEvent = treeReader->UseBranch("Event");
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  
  // Book histograms
  
  TH1 *histGenJetEta = new TH1F("genjet_eta","GenJet Eta Distribution; Eta; Count", 100, -4, 5);
  TH1 *histMatchedGenJetEta = new TH1F("matched_genjet_eta","Matched GenJet Eta Distribution; Eta; Count", 100, -4, 5);
  TH1 *histMatchEffEta = new TH1F("match_eff_eta","Matching Efficiency vs Eta; Eta; Matching Efficiency (Matched GenJets / Total GenJets", 100, -4, 5);
  TH1 *histElectronEta = new TH1F("electron_eta","Electron Eta; Eta; Count", 100, -4, 5);
  int jetEntries = 0;
  int currentJet = 0;

  // Loop over all events
  Jet *genjet;
  Jet *jet;
  Jet *matchedJet;
  Jet *matchedGenJet;
  Jet *unjet;
  Jet *ungenjet;
  GenParticle *particle;
  for(Int_t entry = 0; entry < numberOfEntries; ++entry){
    // Load selected branches with data from specified event
	treeReader->ReadEntry(entry);
	
	int Jets = branchJet->GetEntriesFast();
   	int GenJets = branchGenJet->GetEntriesFast();
	int GenJetEntries[GenJets];
	int ElectronIndex = -1;
	float MaxElectronEnergy = 0;
	float ElectronEta = 0;
	for(int i = 0; i < GenJets; i++){
		GenJetEntries[i] = 0;
		genjet = (Jet*) branchGenJet->At(i);
		for(int j = 0; j < genjet->Particles.GetEntriesFast(); j++){
			particle = (GenParticle*) genjet->Particles.At(j);
			if(particle->PID == 11 && particle->E > MaxElectronEnergy){
				MaxElectronEnergy = particle->E;
				ElectronEta = particle->Eta;
				ElectronIndex = i;
			}
		}
	}
	if(ElectronIndex >= 0){
		GenJetEntries[ElectronIndex] = -1;
		histElectronEta->Fill(ElectronEta);
	}

	int JetEntries[Jets];
	for(int i = 0; i < Jets; i++){
                JetEntries[i] = 0;
        }

	float minDeltaR = 1;
	int jetIndex = 999;
	int genjetIndex = 999;
	int match = 1;

	//Find matches while there are pairs left to match
	while(Jets!= 0 && GenJets != 0){
		jetIndex = -1;
		genjetIndex = -1;
		minDeltaR = 1;
		//Loop over all unpaired reconstructed jets
		for(int i = 0; i < branchJet->GetEntriesFast(); i++){
			if(JetEntries[i] == 0){
				//Loop over all unpaired generated jets
				for(int j = 0; j < branchGenJet->GetEntriesFast(); j++){
					if(GenJetEntries[j] == 0){
						jet = (Jet*) branchJet->At(i);
        					auto jetMomentum = jet->P4();
						genjet = (Jet*) branchGenJet->At(j);
               					auto genJetMomentum = genjet->P4();

						//Find the closest match between all unpaired jets in the event
						if(genJetMomentum.DeltaR(jetMomentum) < minDeltaR){
							minDeltaR = genJetMomentum.DeltaR(jetMomentum);
							jetIndex = i;
							genjetIndex = j;
						}					
					}
				}
			}
		}
		GenJets--;
		Jets--;
		if(genjetIndex >= 0 && jetIndex >= 0){
			GenJetEntries[genjetIndex] = match;
        		JetEntries[jetIndex] = match;
			match++;
		}	
    	}
 
    float genjetEnergy = 0;
    float jetEnergy = 0;
    float deltaR = 0;
    
    //Find the indexes of all matched pairs
    for(int i = 0; i < sizeof(JetEntries) / sizeof(JetEntries[0]); i++){
	for(int j = 0; j < sizeof(GenJetEntries) / sizeof(GenJetEntries[0]); j++){
		if(JetEntries[i] == GenJetEntries[j] && JetEntries[i] > 0 && GenJetEntries[j] > 0){
			
			matchedJet = (Jet*) branchJet->At(i);
                        auto jetMomentum = matchedJet->P4();
                        matchedGenJet = (Jet*) branchGenJet->At(j);
                        auto genJetMomentum = matchedGenJet->P4();
			
			deltaR = genJetMomentum.DeltaR(jetMomentum);
			
			//Do analysis on matched jets here
			genjetEnergy = matchedGenJet->PT * TMath::CosH(matchedGenJet->Eta);			

			//Fill histograms
			histMatchEffEta->Fill(matchedGenJet->Eta);
			histGenJetEta->Fill(matchedGenJet->Eta);
			histMatchedGenJetEta->Fill(matchedGenJet->Eta);
		}
	}
    }

    //Find the indexes of unmatched genjets
    for(int i = 0; i < sizeof(GenJetEntries) / sizeof(GenJetEntries[0]); i++){
	if(GenJetEntries[i] == 0){
		ungenjet = (Jet*) branchGenJet->At(i);

		//Do analysis on unmatched genjets here
		auto genJetMomentum = ungenjet->P4();
		genjetEnergy = ungenjet->PT * TMath::CosH(ungenjet->Eta);	

		//Fill histograms
		histGenJetEta->Fill(ungenjet->Eta);
	}
    }
    
    //Find the indexes of unmatched reconstructed jets
    for(int i = 0; i < sizeof(JetEntries) / sizeof(JetEntries[0]); i++){
        if(JetEntries[i] == 0){
                unjet = (Jet*) branchJet->At(i);

                //Do analysis on unmatched reconstructed jets here
                auto genJetMomentum = ungenjet->P4();
                jetEnergy = unjet->PT * TMath::CosH(unjet->Eta);

                //Fill histograms
        }
    } 
  }

  // Show resulting histograms  
  auto c=new TCanvas("Canvas","Canvas",1500,700);
  c->Divide(2,2);
  c->SetLeftMargin(0.15);
  for(int i = 1; i <= 4; i++){
   c->cd(i)->SetLeftMargin(0.15);
   c->cd(i)->SetRightMargin(0.15);
  }
  c->cd(1);
  histGenJetEta->Draw();
  c->cd(2);
  histMatchedGenJetEta->Draw();
  c->cd(3);
  histMatchEffEta->Divide(histGenJetEta);
  histMatchEffEta->Draw(); 
  c->cd(4);
  histElectronEta->Draw();
}
