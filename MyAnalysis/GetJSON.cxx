#include <iostream>
#include <iomanip>
#include <map>
#include <algorithm>
#include <functional>
#include <vector>

#include "TFile.h"
#include "TTree.h"

#include "MyAnalysis.h"
#include "Pax.h"
#include <json/json.h>

using namespace std;

TTree* MyAnalysis::GetJSONInformation(Json::Value json, vector<int> leaky_pmt_id)
{
   TTree* tree = new TTree("json","json");

   Float_t e_lifetime;
   Float_t drift_velocity;
   vector<Float_t> pmt_gain;
   vector<Float_t> pmt_qe;

   tree->Branch("e_lifetime", &e_lifetime, "e_lifetime/F");
   tree->Branch("drift_velocity", &drift_velocity, "drift_velocity/F");
   tree->Branch("pmt_gain", &pmt_gain);
   tree->Branch("pmt_qe", &pmt_qe);

   e_lifetime     = json["configuration"]["DEFAULT"]["electron_lifetime_liquid"].asFloat();
   drift_velocity = json["configuration"]["DEFAULT"]["drift_velocity_liquid"].asFloat();

   cout<<e_lifetime<<" "<<drift_velocity<<endl;

   for(int i=0; i<260; i++)
   {
      pmt_gain.push_back(json["configuration"]["DEFAULT"]["gains"][i].asFloat());
   }
   for(int i=0; i<259; i++)
   {
      pmt_qe.push_back(json["configuration"]["DEFAULT"]["quantum_efficiencies"][i].asFloat());
   }
  
   tree->Fill();
  	
   return tree; 
}
