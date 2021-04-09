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
  
  // Book histograms
  
  TH1 *histDeltaR = new TH1F("delta_R","Delta R; Delta R; Count", 100, 0, 1);
  TH2 *histEDeltaR = new TH2F("energy_delta_R","GenJet Energy vs Delta R; Genjet Energy; Delta R", 100, 0, 200, 100, 0, 1);

  TH1 *histUnGenJets = new TH1F("unmatched_genjet_energy","Unmatched GenJet Energy; Energy; Count", 100, 0, 200);
  TH1 *histUnJets = new TH1F("unmatched_jet_energy","Unmatched Reconstructed Jet Energy; Energy; Count", 100, 0, 200);

  int jetEntries = 0;
  int currentJet = 0;

  // Loop over all events
  Jet *genjet;
  Jet *jet;
  Jet *matchedJet;
  Jet *matchedGenJet;
  Jet *unjet;
  Jet *ungenjet;

  for(Int_t entry = 0; entry < numberOfEntries; ++entry){
    // Load selected branches with data from specified event
	treeReader->ReadEntry(entry);
	
	int Jets = branchJet->GetEntriesFast();
   	int GenJets = branchGenJet->GetEntriesFast();
	int GenJetEntries[GenJets];

	for(int i = 0; i < GenJets; i++){
		GenJetEntries[i] = 0;
	}

	int JetEntries[Jets];

	for(int i = 0; i < Jets; i++){
		JetEntries[i] = 0;
	}

	float minDeltaR = 999;
	int jetIndex = 999;
	int genjetIndex = 999;
	int match = 1;

	//Find matches while there are pairs left to match
	while(Jets!= 0 && GenJets != 0){
		jetIndex = 999;
		genjetIndex = 999;
		minDeltaR = 999;
		//Loop over all unpaired reconstructed jets
		for(int i = 0; i < branchJet->GetEntriesFast(); i++){
			if(JetEntries[i] < 1){
				//Loop over all unpaired generated jets
				for(int j = 0; j < branchGenJet->GetEntriesFast(); j++){
					if(GenJetEntries[j] < 1){
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
	GenJetEntries[genjetIndex] = match;
        JetEntries[jetIndex] = match;
	match++;	
    }
 
    float genjetEnergy = 0;
    float jetEnergy = 0;
    float deltaR = 0;
    
    //Find the indexes of all matched pairs
    for(int i = 0; i < sizeof(JetEntries) / sizeof(JetEntries[0]); i++){
	for(int j = 0; j < sizeof(GenJetEntries) / sizeof(GenJetEntries[0]); j++){
		if(JetEntries[i] == GenJetEntries[j] && JetEntries[i] != 0 && GenJetEntries[j] != 0){
			
			matchedJet = (Jet*) branchJet->At(i);
                        auto jetMomentum = matchedJet->P4();
                        matchedGenJet = (Jet*) branchGenJet->At(j);
                        auto genJetMomentum = matchedGenJet->P4();
			
			deltaR = genJetMomentum.DeltaR(jetMomentum);
			
			//Apply matching cuts here
			if(deltaR <= 1){

				//Do analysis on matched jets here
				genjetEnergy = matchedGenJet->PT * TMath::CosH(matchedGenJet->Eta);			

				//Fill histograms
				histEDeltaR->Fill(genjetEnergy, deltaR);
				histDeltaR->Fill(deltaR);
			}
			else{
				JetEntries[i] = 0;
				GenJetEntries[j] = 0;
			}
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
		histUnGenJets->Fill(genjetEnergy);
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
                histUnJets->Fill(jetEnergy);
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
  c->cd(1)->SetLogy();
  histDeltaR->Draw();
  c->cd(2);
  histDeltaR->Draw();
  histEDeltaR->Draw("Colz");
  c->cd(3);
  c->cd(3)->SetLogy();
  histUnGenJets->Draw();
  c->cd(4);
  c->cd(4)->SetLogy();
  histUnJets->Draw();
 
}
