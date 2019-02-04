// **********************
//
// author: Alessandra Cappati
//         04/02/2019
//
// run with:
//        root -l -b -q 'read_tree.C("input_file_path_and_name.root", "output_directory_for_plots")'
// e.g.
//        root -l -b -q 'read_tree.C("../L1MuPhase2Ntuple_output.root","outputPlots")'
//
// **********************

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <iomanip>

#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TSystem.h"
#include "TTree.h"
#include "TCanvas.h"

using namespace std;


void read_tree(TString inputFileName = "../L1MuPhase2Ntuple_output.root", TString plotOutputPath = "outputPlots"){


  // define variables
  short int l1mu_Nmuons   = 0;
  vector<float> *l1mu_et  = 0;  // l1mu pT
  vector<float> *l1mu_eta = 0;
  vector<float> *l1mu_phi = 0;


  // open input tree
  TFile* inputFile = new TFile(inputFileName);
  TTree* inputTree = (TTree*)inputFile->Get("L1MuGlobalNtupleMaker/mytree");

  
  // associate branch address
  inputTree->SetBranchAddress("l1mu_Nmuons", &l1mu_Nmuons);
  inputTree->SetBranchAddress("l1mu_et",     &l1mu_et);  
  inputTree->SetBranchAddress("l1mu_eta",    &l1mu_eta);
  inputTree->SetBranchAddress("l1mu_phi",    &l1mu_phi);


  // define histos to fill
  TH1F* hist_l1mu_et  = new TH1F("hist_l1mu_et",  "hist_l1mu_et",  100,  0.,  200.);
  TH1F* hist_l1mu_eta = new TH1F("hist_l1mu_eta", "hist_l1mu_eta", 100, -3.,  3.);
  TH1F* hist_l1mu_phi = new TH1F("hist_l1mu_phi", "hist_l1mu_phi", 100, -3.2, 3.2);


  // --- loop over tree entries 
  Long64_t entries = inputTree->GetEntries();
  for(Long64_t z=0; z<entries; z++){

    inputTree->GetEntry(z);


    // fill histos 
    for(int i = 0; i<l1mu_Nmuons; i++){
      hist_l1mu_et->Fill(l1mu_et->at(i)); 
      hist_l1mu_eta->Fill(l1mu_eta->at(i)); 
      hist_l1mu_phi->Fill(l1mu_phi->at(i)); 
    }
    

  }// end loop over tree entries



  // create output directory for plots
  gSystem->Exec(("mkdir -p "+plotOutputPath));


  // draw histos
  TCanvas* c_l1mu_et = new TCanvas("c_l1mu_et", "c_l1mu_et");
  c_l1mu_et->cd();
  hist_l1mu_et->Draw("hist");
  c_l1mu_et->SaveAs((plotOutputPath + "/" + c_l1mu_et->GetName() + ".png"));

  TCanvas* c_l1mu_eta = new TCanvas("c_l1mu_eta", "c_l1mu_eta");
  c_l1mu_eta->cd();
  hist_l1mu_eta->Draw("hist");
  c_l1mu_eta->SaveAs((plotOutputPath + "/" + c_l1mu_eta->GetName() + ".png"));

  TCanvas* c_l1mu_phi = new TCanvas("c_l1mu_phi", "c_l1mu_phi");
  c_l1mu_phi->cd();
  hist_l1mu_phi->Draw("hist");
  c_l1mu_phi->SaveAs((plotOutputPath + "/" + c_l1mu_phi->GetName() + ".png"));


  // delete pointers
  delete hist_l1mu_et;  
  delete hist_l1mu_eta; 
  delete hist_l1mu_phi; 

  delete c_l1mu_et;
  delete c_l1mu_eta;
  delete c_l1mu_phi;

  delete inputTree;
  delete inputFile;
  
}
