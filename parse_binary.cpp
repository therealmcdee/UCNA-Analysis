#include <TTree.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct events{
  unsigned int nasic;
  unsigned int t0;
  ULong64_t    ntime;
  Double_t     rtime;
  unsigned int nevnt;
  unsigned int lg[32];
  unsigned int hg[32];
  unsigned int hit[32];
};


void zero_evnt(events ev){
  ev.nasic = 0;
  ev.t0    = 0;
  ev.ntime = 0;
  ev.rtime = 0;
  ev.nevnt = 0;

  for( int i = 0 ; i < 32 ; i++){
    ev.lg[i] = 0;
    ev.hg[i] = 0;
    ev.hit[i]= 0;
  }

}



void doit(int n){

  FILE *fdata;
  fdata = fopen(Form("data_files/%d.data",n),"r");

  // buff to hold data from the binary file
  unsigned int buff;
  int  ecnt   = 1;
  bool header = false;
  int  nline  = 0;
  
  unsigned int nasic;
  ULong64_t    ltime,mtime;
  unsigned int last_t0=0;

  // structure for the asic data
  events asic_event[4];
  // Output ROOT File and Tree
  TFile *fOut  = new TFile(Form("data_files/dt_data_%d.root",n),"RECREATE");
  TTree *tData = new TTree("dt55","Data for DT5550w");
  
  // Define the branches
  tData->Branch("asicA",&asic_event[0],"nasic/i:t0:ntime/l:rtime/D:nevnt/i:lg[32]:hg[32]:hit[32]");
  tData->Branch("asicB",&asic_event[1],"nasic/i:t0:ntime/l:rtime/D:nevnt/i:lg[32]:hg[32]:hit[32]");
  tData->Branch("asicC",&asic_event[2],"nasic/i:t0:ntime/l:rtime/D:nevnt/i:lg[32]:hg[32]:hit[32]");
  tData->Branch("asicD",&asic_event[3],"nasic/i:t0:ntime/l:rtime/D:nevnt/i:lg[32]:hg[32]:hit[32]");
  
  if ( fdata != NULL){ 
    std::cout << "File is open \n";
    // start with the structure zeroed out
    for (int i = 0 ; i < 4 ; i++)
      {zero_evnt(asic_event[i]);}
	
    while(!feof(fdata)){
      // read a 32-bit chunk of data
      size_t read = fread(&buff,sizeof(buff),1,fdata);

      // check for the data packet header
      if ( ((buff >> 4 )) == 0x8000000 && !(header)){
	std::cout << "Found a header" << std::endl;
	std::cout << "At event : " << ecnt << std::endl;
	ecnt++;
	nline++;
	header = true;
	nasic = (buff & 0xf);
	asic_event[nasic].nasic = nasic;
	std::cout << "ASIC ID :" << nasic << std::endl;
      } else if ( buff == 0xc0000000){
	// Check for the footer
	std::cout << "Found a footer" << std::endl;
	nline = 0;
	header=false;
      } else if (header && nline > 0){
	// If we've found a header and incrimented the line counter
	// then start parseing the event data.
	switch (nline){
	case 1 :
	  // Get the t0 counter
	  asic_event[nasic].t0 = buff;
	  nline++;
	  std::cout << "T0 counter : " << asic_event[nasic].t0 << std::endl;
	  
	  if ( asic_event[nasic].t0 - last_t0 < 10 ) {
	    std::cout << "Found a multi-asic event" << std::endl;
	    // if the difference in the counter from the last event is
	    // small then cluster the events.  Typically this is 1 cycle
	    // off for each asic that triggered.  Setting a threshold of
	    // 10 just to be sure.
	    // !!!! Should look at this behavior at exetremely high rates
	  } else {
	    // if the next event occurs after the cluster time then
	    // fill the tree and re-zero the data structure.
	    tData->Fill();
	    std::cout << "New Event\n";
	    for (int i = 0 ; i < 4 ; i++)
	      {zero_evnt(asic_event[i]);}
	  }
	  last_t0 = asic_event[nasic].t0;
	  break;
	case 2:
	  // Store the "most significant bits" of the timestamp.
	  mtime = buff;
	  nline++;
	  break;
	case 3:
	  // store the least significant bits of the timestamp.
	  ltime = buff;
	  // shift the MSB over 32 bits and stitch together the full time stamp.
	  asic_event[nasic].ntime = ((mtime << 32) | ltime);
	  asic_event[nasic].rtime = asic_event[nasic].ntime * 8e-10;
	  nline++;
	  break;
	case 4:
	  asic_event[nasic].nevnt = buff;
	  nline++;
	  std::cout << "Event number : " << asic_event[nasic].nevnt << std::endl;
	  break;
	default:
	  if(nline > 4 && nline < 37){
	    std::cout << " Nline is : " << nline-5 << std::endl;
	    asic_event[nasic].hg[nline-5] = (buff & 0x3fff);
	    asic_event[nasic].lg[nline-5] = ((buff >> 14) & 0x3fff);
	    asic_event[nasic].hit[nline-5] =  ((buff >> 28) & 0x1);
	    nline++;
	  }
	}
      }
    }
  } else {
    std::cout << "File not opened!\n";
  }
  
  fclose(fdata);
  fOut->Write("",TObject::kOverwrite);
  fOut->Close();

}
