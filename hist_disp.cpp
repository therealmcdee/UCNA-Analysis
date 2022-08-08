#include "TFile.h"
#include "TTree.h"


struct events{
	unsigned int nasic;
	unsigned int t0;
	ULong64_t ntime;
	Double_t rtime;
	unsigned int nevnt;
	unsigned int lg[32];
	unsigned int hg[32];
	unsigned int hit[32];
};



void hist_gain(char asic[1], char subj[2], int file){
  TFile *fdata;
  fdata = new TFile(Form("data_files/dt_data_%d.root",file),"READ");

  TTree* tree = (TTree*) fdata->Get("dt55");
  
  events hug;

  tree->SetBranchAddress(Form("asic%s",asic), &hug.nasic);

  TH1D *hg[32];

  for (int i = 0; i< 32; i++) {
	  hg[i] = new TH1D(Form("hg[%d]",i),Form("Chn %d;ADC;Counts",i),1024,0,16384);
  };
  for (int i= 0; i<tree->GetEntries();i++) {
	  tree->GetEntry(i);
	  for (int j = 0; j<32; j++) {
		  if (strcmp(subj,"hg") == 0) {
		  	hg[j]->Fill(hug.hg[j]);
		  } else if (strcmp(subj,"lg") == 0) {
			  hg[j]->Fill(hug.lg[j]);
		  };
	  };
  };

  TCanvas *canvas = new TCanvas("canvas","pads");
  canvas->Divide(4,8);
  for (int i = 1; i<33; i++){
	  //if (i==0) {
	  //	  hg[i]->Draw();
	  //} else {
	  //	  hg[i]->SetLineColor(i+10);
	  //	  hg[i]->Draw("same");
	  //};
	  canvas->cd(i);
	  hg[i-1]->Draw();
  }
}

//------------------------------------------------------------//

void evt_analysis(int file) {
	events asic[4];
	TFile *fdata;
	fdata = new TFile(Form("data_files/dt_data_%d.root",file),"READ");
	TTree *t = (TTree*) fdata->Get("dt55");
	char names[5] = "ABCD";

	for (int j = 0; j<4; j++) {
		t->SetBranchAddress(Form("asic%c",names[j]),&asic[j].nasic);
	};

	TH1D *sipm[4];
	TH1D *mult = new TH1D("mult","Multiplicity of SiPM Triggers",128,0,128);
	for (int i = 0; i<4; i++) {
		sipm[i] = new TH1D(Form("asic%c",names[i]),Form("multiplicity%c",names[i]),32,0,32);
		}

	//Print(Form("TOTAL EVENTS = %lld",t->GetEntries()));
	for (int j = 0; j<t->GetEntries();j++) {
		int multiplicity = 0;
		t->GetEntry(j);
		for (int k = 0; k<4; k++) {
			int mult_indi = 0;
			for (int w = 0; w<32; w++) {
				if (asic[k].hg[w] > 960) {
					multiplicity = multiplicity + 1;
					mult_indi = mult_indi + 1;
				}
			}
			sipm[k]->Fill(mult_indi);
		}	
		mult->Fill(multiplicity);
	}
		
	TCanvas *c1 = new TCanvas("c1","multi");
	mult->Draw();
	TCanvas *c2 = new TCanvas("c2","individual");
	c2->Divide(2,2);
	for (int i = 0; i<4; i++) {
		c2->cd(i+1);
		sipm[i]->Draw();
	}
}
