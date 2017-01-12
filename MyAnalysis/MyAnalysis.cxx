#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <algorithm>
#include <functional>

#include "MyAnalysis.h"
#include "Pax.h"
#include "Dict.cc"
#include <json/json.h>

using namespace std;

typedef pair<float, int> sort_t;

void MyAnalysis::Loop()
{
   initialize();

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();

   for (Long64_t jentry=0; jentry<nentries; jentry++) 
   {
      if(jentry%1000==0) cout << "Event No." << jentry << " / " << nentries << endl;
   
      fChain->GetEntry(jentry);  

      clear();
      execute();
   }

   finalize();
}

void MyAnalysis::initialize(){

   gROOT->ProcessLine("#include <vector>");

   /// Get Leaky PMT information ///
   get_leaky_PMTs(leaky_pmts,ap_rate);

   for(int i=0; i<leaky_pmts.size(); i++)
   {
      if(i==0)cout<<"Leaky PMTs = "<<leaky_pmts.at(i);
      if(i!=0)cout<<", "<<leaky_pmts.at(i);
   }
   cout<<""<<endl;

   /// Get Pax Information ///
   events = new Event();
   fChain->SetBranchAddress("events", &events);

   /// Definition of output TFile ///
   m_file = new TFile(OutputFileName,"recreate");
   cout<<"Output filename = "<<OutputFileName<<endl;

   /// Get JSON Information and SaveAs TTree /// 
   ifstream pax_json("pax_metadata.json", ifstream::binary);
   bool parsingSuccessful = reader.parse(pax_json, json, false);
   json_tree = GetJSONInformation(json,leaky_pmts);

   /// Definition of output TTree ///
   out_tree = new TTree("tree", "tree");

   /// General information ///
   out_tree->Branch("Filename", &Filename);
   out_tree->Branch("EventId", &EventId, "EventId/I");
   out_tree->Branch("time", &time, "time/L");

   /// Interactions ///
   out_tree->Branch("n_int", &n_int, "n_int/I");
   out_tree->Branch("int_s1_id", &int_s1_id);
   out_tree->Branch("int_s2_id", &int_s2_id);
   out_tree->Branch("int_s1_area_correction", &int_s1_area_correction);
   out_tree->Branch("int_s1_pattern_fit", &int_s1_pattern_fit);
   out_tree->Branch("int_s1_saturation_correction", &int_s1_saturation_correction);
   out_tree->Branch("int_s1_spatial_correction", &int_s1_spatial_correction);
   out_tree->Branch("int_s2_area_correction", &int_s2_area_correction);
   out_tree->Branch("int_s2_lifetime_correction", &int_s2_lifetime_correction);
   out_tree->Branch("int_x", &int_x);
   out_tree->Branch("int_y", &int_y);
   out_tree->Branch("int_z", &int_z);
   out_tree->Branch("int_xy_posrec_algorithm", &int_xy_posrec_algorithm);
   out_tree->Branch("int_xy_posrec_goodness_of_fit", &int_xy_posrec_goodness_of_fit);
   out_tree->Branch("int_xy_posrec_ndf", &int_xy_posrec_ndf);
   out_tree->Branch("int_drift_time", &int_drift_time);

   /// S1 ///
   out_tree->Branch("n_s1", &n_s1, "n_s1/I");
   out_tree->Branch("s1_id", &s1_id);
   out_tree->Branch("s1_area", &s1_area);
   out_tree->Branch("s1_area_fraction_top", &s1_area_fraction_top);
   out_tree->Branch("s1_area_fraction_leaky_pmts", &s1_area_fraction_leaky_pmts);
   out_tree->Branch("s1_height", &s1_height);
   out_tree->Branch("s1_top_hitpattern_spread", &s1_top_hitpattern_spread);
   out_tree->Branch("s1_bottom_hitpattern_spread", &s1_bottom_hitpattern_spread);
   out_tree->Branch("s1_center_time", &s1_center_time);
   out_tree->Branch("s1_hit_time_mean", &s1_hit_time_mean);
   out_tree->Branch("s1_hit_time_std", &s1_hit_time_std);
   out_tree->Branch("s1_n_hits", &s1_n_hits);
   out_tree->Branch("s1_hits_fraction_top", &s1_hits_fraction_top);
   out_tree->Branch("s1_hits_fraction_leaky_pmts", &s1_hits_fraction_leaky_pmts);
   out_tree->Branch("s1_n_contributing_channels", &s1_n_contributing_channels);
   out_tree->Branch("s1_n_contributing_channels_top", &s1_n_contributing_channels_top);
   out_tree->Branch("s1_n_saturated_channels", &s1_n_saturated_channels);
   out_tree->Branch("s1_range_area_decile_10p", &s1_range_area_decile_10p);
   out_tree->Branch("s1_range_area_decile_50p", &s1_range_area_decile_50p);
   out_tree->Branch("s1_n_noise_pulses", &s1_n_noise_pulses);
   out_tree->Branch("s1_x", &s1_x);
   out_tree->Branch("s1_y", &s1_y);
   out_tree->Branch("s1_goodness_of_fit", &s1_goodness_of_fit);
   out_tree->Branch("s1_ndf", &s1_ndf);
   out_tree->Branch("s1_algorithm", &s1_algorithm);

   /// S2 ///
   out_tree->Branch("n_s2", &n_s2, "n_s2/I");
   out_tree->Branch("s2_id", &s2_id);
   out_tree->Branch("s2_area", &s2_area);
   out_tree->Branch("s2_area_fraction_top", &s2_area_fraction_top);
   out_tree->Branch("s2_area_fraction_leaky_pmts", &s2_area_fraction_leaky_pmts);
   out_tree->Branch("s2_height", &s2_height);
   out_tree->Branch("s2_top_hitpattern_spread", &s2_top_hitpattern_spread);
   out_tree->Branch("s2_bottom_hitpattern_spread", &s2_bottom_hitpattern_spread);
   out_tree->Branch("s2_center_time", &s2_center_time);
   out_tree->Branch("s2_hit_time_mean", &s2_hit_time_mean);
   out_tree->Branch("s2_hit_time_std", &s2_hit_time_std);
   out_tree->Branch("s2_n_hits", &s2_n_hits);
   out_tree->Branch("s2_hits_fraction_top", &s2_hits_fraction_top);
   out_tree->Branch("s2_hits_fraction_leaky_pmts", &s2_hits_fraction_leaky_pmts);
   out_tree->Branch("s2_n_contributing_channels", &s2_n_contributing_channels);
   out_tree->Branch("s2_n_contributing_channels_top", &s2_n_contributing_channels_top);
   out_tree->Branch("s2_n_saturated_channels", &s2_n_saturated_channels);
   out_tree->Branch("s2_range_area_decile_10p", &s2_range_area_decile_10p);
   out_tree->Branch("s2_range_area_decile_50p", &s2_range_area_decile_50p);
   out_tree->Branch("s2_n_noise_pulses", &s2_n_noise_pulses);
   out_tree->Branch("s2_x", &s2_x);
   out_tree->Branch("s2_y", &s2_y);
   out_tree->Branch("s2_goodness_of_fit", &s2_goodness_of_fit);
   out_tree->Branch("s2_ndf", &s2_ndf);
   out_tree->Branch("s2_algorithm", &s2_algorithm);

   /// Other Peaks (maybe noise or afterpulse) /// 
   out_tree->Branch("n_other_peaks", &n_other_peaks, "n_other_peaks/I");
   out_tree->Branch("other_peaks_area", &other_peaks_area);
   out_tree->Branch("other_peaks_area_fraction_top", &other_peaks_area_fraction_top);
   out_tree->Branch("other_peaks_area_fraction_leaky_pmts", &other_peaks_area_fraction_leaky_pmts);
   out_tree->Branch("other_peaks_height", &other_peaks_height);
   out_tree->Branch("other_peaks_top_hitpattern_spread", &other_peaks_top_hitpattern_spread);
   out_tree->Branch("other_peaks_bottom_hitpattern_spread", &other_peaks_bottom_hitpattern_spread);
   out_tree->Branch("other_peaks_center_time", &other_peaks_center_time);
   out_tree->Branch("other_peaks_hit_time_mean", &other_peaks_hit_time_mean);
   out_tree->Branch("other_peaks_hit_time_std", &other_peaks_hit_time_std);
   out_tree->Branch("other_peaks_n_hits", &other_peaks_n_hits);
   out_tree->Branch("other_peaks_hits_fraction_top", &other_peaks_hits_fraction_top);
   out_tree->Branch("other_peaks_hits_fraction_leaky_pmts", &other_peaks_hits_fraction_leaky_pmts);
   out_tree->Branch("other_peaks_n_contributing_channels", &other_peaks_n_contributing_channels);
   out_tree->Branch("other_peaks_n_contributing_channels_top", &other_peaks_n_contributing_channels_top);
   out_tree->Branch("other_peaks_n_saturated_channels", &other_peaks_n_saturated_channels);
   out_tree->Branch("other_peaks_range_area_decile_10p", &other_peaks_range_area_decile_10p);
   out_tree->Branch("other_peaks_range_area_decile_50p", &other_peaks_range_area_decile_50p);
   out_tree->Branch("other_peaks_n_noise_pulses", &other_peaks_n_noise_pulses);
   out_tree->Branch("other_peaks_x", &other_peaks_x);
   out_tree->Branch("other_peaks_y", &other_peaks_y);
   out_tree->Branch("other_peaks_goodness_of_fit", &other_peaks_goodness_of_fit);
   out_tree->Branch("other_peaks_ndf", &other_peaks_ndf);
   out_tree->Branch("other_peaks_algorithm", &other_peaks_algorithm);

}

void MyAnalysis::execute(){

   /// Get Pax Information ///
   vector <Interaction> interactions = events->interactions;
   vector <Peak> peaks = events->peaks;

   /// Event Information ///
   EventId  = events->event_number;
   Filename = events->dataset_name; 
   time     = events->start_time;

   Float_t s1=0.;
   Float_t s2=0.;

   vector<int> s1_index;
   vector<int> s2_index;
   vector<float> s1_s2_index;
 
   /// Interaction Information ///
   for(int i=0; i<interactions.size(); i++)
   {  
       ++n_int;
       int_drift_time.push_back(interactions.at(i).drift_time);
       int_s1_id.push_back(interactions.at(i).s1);
       int_s2_id.push_back(interactions.at(i).s2);
       int_s1_area_correction.push_back(interactions.at(i).s1_area_correction);
       int_s1_pattern_fit.push_back(interactions.at(i).s1_pattern_fit);
       int_s1_saturation_correction.push_back(interactions.at(i).s1_saturation_correction);
       int_s1_spatial_correction.push_back(interactions.at(i).s1_spatial_correction);
       int_s2_area_correction.push_back(interactions.at(i).s2_area_correction);
       int_s2_lifetime_correction.push_back(interactions.at(i).s2_area_correction);
       int_x.push_back(interactions.at(i).x);
       int_y.push_back(interactions.at(i).y);
       int_z.push_back(interactions.at(i).z);
       int_xy_posrec_algorithm.push_back(interactions.at(i).xy_posrec_algorithm);
       int_xy_posrec_goodness_of_fit.push_back(interactions.at(i).xy_posrec_goodness_of_fit);
       int_xy_posrec_ndf.push_back(interactions.at(i).xy_posrec_ndf);

       /// For avoiding double count of S1 ///
       bool IsUsedS1=false;

       for(unsigned int j=0; j<s1_index.size(); j++)
       {
          if(interactions.at(i).s1==s1_index.at(j))
          {
             IsUsedS1=true;
          }
       }

       /// For avoiding double count of S2 ///
       bool IsUsedS2=false;

       for(unsigned int j=0; j<s2_index.size(); j++)
       {
          if(interactions.at(i).s2==s2_index.at(j))
          {
             IsUsedS2=true;
          }
       }

       /// S1 Information ///
       if(IsUsedS1==false)
       {
	  ++n_s1;

	  s1_index.push_back(interactions.at(i).s1);
	  s1_s2_index.push_back(interactions.at(i).s1);

          int id_s1 = interactions.at(i).s1;
          s1 = peaks.at(id_s1).area;

          s1_id.push_back(id_s1);
          s1_area.push_back(s1);
          s1_area_fraction_top.push_back(peaks.at(id_s1).area_fraction_top);
          s1_height.push_back(peaks.at(id_s1).height);
          s1_top_hitpattern_spread.push_back(peaks.at(id_s1).top_hitpattern_spread);
          s1_bottom_hitpattern_spread.push_back(peaks.at(id_s1).bottom_hitpattern_spread);
          s1_center_time.push_back(peaks.at(id_s1).center_time);
          s1_hit_time_mean.push_back(peaks.at(id_s1).hit_time_mean);
          s1_hit_time_std.push_back(peaks.at(id_s1).hit_time_std);
          s1_n_hits.push_back(peaks.at(id_s1).n_hits);
          s1_hits_fraction_top.push_back(peaks.at(id_s1).hits_fraction_top);
          s1_n_contributing_channels.push_back(peaks.at(id_s1).n_contributing_channels);
          s1_n_contributing_channels_top.push_back(peaks.at(id_s1).n_contributing_channels_top);
          s1_n_saturated_channels.push_back(peaks.at(id_s1).n_saturated_channels);
          s1_range_area_decile_10p.push_back(peaks.at(id_s1).range_area_decile[1]);
          s1_range_area_decile_50p.push_back(peaks.at(id_s1).range_area_decile[5]);
          s1_n_noise_pulses.push_back(peaks.at(id_s1).n_noise_pulses);
          s1_area_fraction_leaky_pmts.push_back(get_area_fraction_leaky_PMTs(peaks.at(id_s1).area_per_channel, s1, leaky_pmts));
          s1_hits_fraction_leaky_pmts.push_back(get_hits_fraction_leaky_PMTs(peaks.at(id_s1).hits_per_channel, peaks.at(id_s1).n_hits, leaky_pmts));

          std::vector <ReconstructedPosition> rec_pos = peaks.at(id_s1).reconstructed_positions;

          vector<Float_t> x;
          vector<Float_t> y;
          vector<Float_t> goodness_of_fit;
          vector<Float_t> ndf;
          vector<TString> algorithm;

          for(int j=0; j<rec_pos.size(); j++)
          {
             x.push_back(rec_pos.at(j).x);
             y.push_back(rec_pos.at(j).y);
             goodness_of_fit.push_back(rec_pos.at(j).goodness_of_fit);
             ndf.push_back(rec_pos.at(j).ndf);
             algorithm.push_back(rec_pos.at(j).algorithm);
          }

          s1_x.push_back(x);
          s1_y.push_back(y);
          s1_goodness_of_fit.push_back(goodness_of_fit);
          s1_ndf.push_back(ndf);
          s1_algorithm.push_back(algorithm);

       }

       /// S2 Information ///
       if(IsUsedS2==false)
       {
	  ++n_s2;

	  s2_index.push_back(interactions.at(i).s2);
	  s1_s2_index.push_back(interactions.at(i).s2);

          int id_s2 = interactions.at(i).s2;
          s2 = peaks.at(id_s2).area;

          s2_id.push_back(id_s2);
          s2_area.push_back(s2);
          s2_area_fraction_top.push_back(peaks.at(id_s2).area_fraction_top);
          s2_height.push_back(peaks.at(id_s2).height);
          s2_top_hitpattern_spread.push_back(peaks.at(id_s2).top_hitpattern_spread);
          s2_bottom_hitpattern_spread.push_back(peaks.at(id_s2).bottom_hitpattern_spread);
          s2_center_time.push_back(peaks.at(id_s2).center_time);
          s2_hit_time_mean.push_back(peaks.at(id_s2).hit_time_mean);
          s2_hit_time_std.push_back(peaks.at(id_s2).hit_time_std);
          s2_n_hits.push_back(peaks.at(id_s2).n_hits);
          s2_hits_fraction_top.push_back(peaks.at(id_s2).hits_fraction_top);
          s2_n_contributing_channels.push_back(peaks.at(id_s2).n_contributing_channels);
          s2_n_contributing_channels_top.push_back(peaks.at(id_s2).n_contributing_channels_top);
          s2_n_saturated_channels.push_back(peaks.at(id_s2).n_saturated_channels);
          s2_range_area_decile_10p.push_back(peaks.at(id_s2).range_area_decile[1]);
          s2_range_area_decile_50p.push_back(peaks.at(id_s2).range_area_decile[5]);
          s2_n_noise_pulses.push_back(peaks.at(id_s2).n_noise_pulses);
          s2_area_fraction_leaky_pmts.push_back(get_area_fraction_leaky_PMTs(peaks.at(id_s2).area_per_channel, s2, leaky_pmts));
          s2_hits_fraction_leaky_pmts.push_back(get_hits_fraction_leaky_PMTs(peaks.at(id_s2).hits_per_channel, peaks.at(id_s2).n_hits, leaky_pmts));

          std::vector <ReconstructedPosition> rec_pos = peaks.at(id_s2).reconstructed_positions;

          vector<Float_t> x;
          vector<Float_t> y;
          vector<Float_t> goodness_of_fit;
          vector<Float_t> ndf;
          vector<TString> algorithm;

          for(int j=0; j<rec_pos.size(); j++)
          {
             x.push_back(rec_pos.at(j).x);
             y.push_back(rec_pos.at(j).y);
             goodness_of_fit.push_back(rec_pos.at(j).goodness_of_fit);
             ndf.push_back(rec_pos.at(j).ndf);
             algorithm.push_back(rec_pos.at(j).algorithm);
          }
    
          s2_x.push_back(x);
          s2_y.push_back(y);
          s2_goodness_of_fit.push_back(goodness_of_fit);
          s2_ndf.push_back(ndf);
          s2_algorithm.push_back(algorithm);
       }
    
   }   
   /// Other Peak Information (sorted by energy) ///
   vector<int> op_id;

   for(unsigned int i=0; i<peaks.size(); i++)
   {
      bool IsOtherPeaks=true;

      for(unsigned int j=0; j<s1_s2_index.size(); j++)
      {
         if(i==s1_s2_index.at(j))
         {
            IsOtherPeaks=false;
         }
      }

      if(IsOtherPeaks==true)
      {
         op_id.push_back(i);
      }
   }

   vector<sort_t> sort_peaks(op_id.size()); 

   for(unsigned int i=0; i<op_id.size(); i++)
   {
      sort_peaks[i].first  = peaks.at(op_id.at(i)).area;
      sort_peaks[i].second = op_id.at(i);

   }
   sort(sort_peaks.begin(), sort_peaks.end(),  greater<pair<float, int> >());

   /// Sorted by energy (in p.e.) ///
   for(unsigned int i=0; i<sort_peaks.size(); i++)
   {
       int id = sort_peaks[i].second; 

       if(peaks.at(id).area>1.0)
       {
          ++n_other_peaks;
          other_peaks_area.push_back(peaks.at(id).area);
          other_peaks_area_fraction_top.push_back(peaks.at(id).area_fraction_top);
          other_peaks_height.push_back(peaks.at(id).height);
          other_peaks_top_hitpattern_spread.push_back(peaks.at(id).top_hitpattern_spread);
          other_peaks_bottom_hitpattern_spread.push_back(peaks.at(id).bottom_hitpattern_spread);
          other_peaks_center_time.push_back(peaks.at(id).center_time);
          other_peaks_hit_time_mean.push_back(peaks.at(id).hit_time_mean);
          other_peaks_hit_time_std.push_back(peaks.at(id).hit_time_std);
          other_peaks_n_hits.push_back(peaks.at(id).n_hits);
          other_peaks_hits_fraction_top.push_back(peaks.at(id).hits_fraction_top);
          other_peaks_n_contributing_channels.push_back(peaks.at(id).n_contributing_channels);
          other_peaks_n_contributing_channels_top.push_back(peaks.at(id).n_contributing_channels_top);
          other_peaks_n_saturated_channels.push_back(peaks.at(id).n_saturated_channels);
          other_peaks_range_area_decile_10p.push_back(peaks.at(id).range_area_decile[1]);
          other_peaks_range_area_decile_50p.push_back(peaks.at(id).range_area_decile[5]);
          other_peaks_n_noise_pulses.push_back(peaks.at(id).n_noise_pulses);
          other_peaks_area_fraction_leaky_pmts.push_back(get_area_fraction_leaky_PMTs(peaks.at(id).area_per_channel, peaks.at(id).area  , leaky_pmts));
          other_peaks_hits_fraction_leaky_pmts.push_back(get_hits_fraction_leaky_PMTs(peaks.at(id).hits_per_channel, peaks.at(id).n_hits, leaky_pmts));

          std::vector <ReconstructedPosition> rec_pos = peaks.at(id).reconstructed_positions;

          vector<Float_t> x;
          vector<Float_t> y;
          vector<Float_t> goodness_of_fit;
          vector<Float_t> ndf;
          vector<TString> algorithm;

          for(int j=0; j<rec_pos.size(); j++)
          {
             x.push_back(rec_pos.at(j).x);
             y.push_back(rec_pos.at(j).y);
             goodness_of_fit.push_back(rec_pos.at(j).goodness_of_fit);
             ndf.push_back(rec_pos.at(j).ndf);
             algorithm.push_back(rec_pos.at(j).algorithm);
          }
          
          other_peaks_x.push_back(x);
          other_peaks_y.push_back(y);
          other_peaks_goodness_of_fit.push_back(goodness_of_fit);
          other_peaks_ndf.push_back(ndf);
          other_peaks_algorithm.push_back(algorithm);
      }
   }    

   /// Fill all variables to TTree event by event ///
   out_tree->Fill();

   /// clear vectors ///
   if(interactions.size()>0)interactions.clear();
   if(peaks.size()>0)peaks.clear();
   if(s1_index.size()>0)s1_index.clear();
   if(s2_index.size()>0)s2_index.clear();
   if(s1_s2_index.size()>0)s1_index.clear();
   if(op_id.size()>0)op_id.clear();
   if(sort_peaks.size()>0)sort_peaks.clear();
}

void MyAnalysis::clear(){

   /// Interactions ///
   if(n_int>0)
   {
      int_s1_id.clear();
      int_s2_id.clear();
      int_s1_area_correction.clear();
      int_s1_pattern_fit.clear();
      int_s1_saturation_correction.clear();
      int_s1_spatial_correction.clear();
      int_s2_area_correction.clear();
      int_s2_lifetime_correction.clear();
      int_x.clear();
      int_y.clear();
      int_z.clear();
      int_xy_posrec_algorithm.clear();
      int_xy_posrec_goodness_of_fit.clear();
      int_xy_posrec_ndf.clear();
      int_drift_time.clear();
   }

   /// S1 ///
   if(n_s1>0)
   {
      s1_id.clear();
      s1_area.clear();
      s1_area_fraction_top.clear();
      s1_area_fraction_leaky_pmts.clear();
      s1_height.clear();
      s1_top_hitpattern_spread.clear();
      s1_bottom_hitpattern_spread.clear();
      s1_center_time.clear();
      s1_hit_time_mean.clear();
      s1_hit_time_std.clear();
      s1_n_hits.clear();
      s1_hits_fraction_top.clear();
      s1_hits_fraction_leaky_pmts.clear();
      s1_n_contributing_channels.clear();
      s1_n_contributing_channels_top.clear();
      s1_n_saturated_channels.clear();
      s1_range_area_decile_10p.clear();
      s1_range_area_decile_50p.clear();
      s1_n_noise_pulses.clear();
      s1_x.clear();
      s1_y.clear();
      s1_goodness_of_fit.clear();
      s1_ndf.clear();
      s1_algorithm.clear();
   }

   /// S2 ///
   if(n_s2>0)
   {
      s2_id.clear();
      s2_area.clear();
      s2_area_fraction_top.clear();
      s2_area_fraction_leaky_pmts.clear();
      s2_height.clear();
      s2_top_hitpattern_spread.clear();
      s2_bottom_hitpattern_spread.clear();
      s2_center_time.clear();
      s2_hit_time_mean.clear();
      s2_hit_time_std.clear();
      s2_n_hits.clear();
      s2_hits_fraction_top.clear();
      s2_hits_fraction_leaky_pmts.clear();
      s2_n_contributing_channels.clear();
      s2_n_contributing_channels_top.clear();
      s2_n_saturated_channels.clear();
      s2_range_area_decile_10p.clear();
      s2_range_area_decile_50p.clear();
      s2_n_noise_pulses.clear();
      s2_x.clear();
      s2_y.clear();
      s2_goodness_of_fit.clear();
      s2_ndf.clear();
      s2_algorithm.clear();
   }

   /// Other Peaks ///
   if(n_other_peaks>0)
   {
      other_peaks_area.clear();
      other_peaks_area_fraction_top.clear();
      other_peaks_area_fraction_leaky_pmts.clear();
      other_peaks_height.clear();
      other_peaks_top_hitpattern_spread.clear();
      other_peaks_bottom_hitpattern_spread.clear();
      other_peaks_center_time.clear();
      other_peaks_hit_time_mean.clear();
      other_peaks_hit_time_std.clear();
      other_peaks_n_hits.clear();
      other_peaks_hits_fraction_top.clear();
      other_peaks_hits_fraction_leaky_pmts.clear();
      other_peaks_n_contributing_channels.clear();
      other_peaks_n_contributing_channels_top.clear();
      other_peaks_n_saturated_channels.clear();
      other_peaks_range_area_decile_10p.clear();
      other_peaks_range_area_decile_50p.clear();
      other_peaks_n_noise_pulses.clear();
      other_peaks_x.clear();
      other_peaks_y.clear();
      other_peaks_goodness_of_fit.clear();
      other_peaks_ndf.clear();
      other_peaks_algorithm.clear();
   }
 
   n_int         = 0;
   n_s1          = 0;
   n_s2          = 0;
   n_other_peaks = 0;
}

void MyAnalysis::finalize(){

   m_file->Write();
   m_file->Close();
   cout<<"Now Finalizing My Analysis"<<endl;
}
