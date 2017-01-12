#include <iostream>
#include <iomanip>
#include <map>
#include <algorithm>
#include <functional>
#include <vector>

#include "MyAnalysis.h"
#include "Pax.h"

using namespace std;

void MyAnalysis::get_leaky_PMTs(vector<int>& leaky_pmts, vector<double>& ap_rate)
{
   TFile *fin = new TFile("data/AP_Run_3907.root", "read");
   TTree *t = (TTree*)fin->Get("t1");
 
   Int_t PMT;
   Double_t totalAP;

   t->SetBranchAddress("PMT", &PMT);
   t->SetBranchAddress("totalAP", &totalAP);

   const Int_t Nevt = t->GetEntries();

   for(Int_t i=0; i<Nevt; i++)
   {
      t->GetEntry(i);

      if(totalAP > 0.5)
      {
         leaky_pmts.push_back(PMT);
         ap_rate.push_back(totalAP);
      }
   }
}

float MyAnalysis::get_area_fraction_leaky_PMTs(Double_t *area_per_channel, Float_t area, vector<int> leaky_pmts)
{
   Double_t area_in_leaky_pmts = 0.;

   for(int i=0; i<leaky_pmts.size(); i++)
   {
      area_in_leaky_pmts += area_per_channel[leaky_pmts.at(i)]; 
   }

   float area_fraction = (float)area_in_leaky_pmts/area;

   return area_fraction;
}

float MyAnalysis::get_hits_fraction_leaky_PMTs(Short_t *hits_per_channel, Int_t hits, vector<int> leaky_pmts)
{
   Short_t hits_in_leaky_pmts = 0.;

   for(int i=0; i<leaky_pmts.size(); i++)
   {
      hits_in_leaky_pmts += hits_per_channel[leaky_pmts.at(i)];
   }

   float hits_fraction = (float)hits_in_leaky_pmts/(float)hits;

   return hits_fraction;
}
