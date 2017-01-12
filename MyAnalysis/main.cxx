#include <fstream>
#include <sstream>      
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>

#include "TApplication.h"
#include "TFile.h"
#include "TString.h"
#include "TChain.h"

#include "MyAnalysis.h"
#include "Pax.h"
#include <json/json.h>

using namespace std;

int main(int argc, char **argv)
{
  TFile *file;
  TChain *chain;

  if(argc==3)
  {
    ifstream fin(argv[1]);
    if(!fin.is_open()){
      cerr << "file open error:" << argv[1] << endl;
      exit(1);
    }

    TString fname;
    int itr=0;
    while(fin >> fname){
      cout << "Input filename["<<itr<<"] = " << fname << endl;
      if(itr==0){
         file = new TFile(fname);
         chain = (TChain*)file->Get("tree");
      }
      if(itr>0)chain->Add(fname);
      ++itr;
    }
  }
  else
  {
    cerr << "Number of input parameters is wrong!!" << endl;
    exit(1);
  }

  string pax_metadata = file->Get("pax_metadata")->GetTitle();
  ofstream fout;
  fout.open("dummy.json");
  fout<<pax_metadata<<endl;
  fout.close();
  system("bash data/json.sh");

  MyAnalysis *ana = new MyAnalysis();
  ana->SetInputChain(chain);
  ana->SetOutputFileName(argv[2]);
  ana->Loop();
  
  return 0;
}
