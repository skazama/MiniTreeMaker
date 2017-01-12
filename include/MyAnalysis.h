#ifndef MyAnalysis_h
#define MyAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TObject.h>
#include <vector>

#include <Pax.h>
#include <json/json.h>

using namespace std;

class MyAnalysis 
{

public:

   void Loop();
   void initialize();
   void execute();
   void finalize();
   void clear();

   void SetOutputFileName(TString f){OutputFileName = f;}
   void SetInputChain(TChain* chain){fChain = chain;}
 
   void get_leaky_PMTs(vector<int>& leaky_pmts, vector<double>& ap_rate);
   float get_area_fraction_leaky_PMTs(Double_t *area_per_channel, Float_t area, vector<int> leaky_pmts);
   float get_hits_fraction_leaky_PMTs(Short_t *hits_per_channel, Int_t hits, vector<int> leaky_pmts);
   TTree* GetJSONInformation(Json::Value json, vector<int> leaky_pmt_id);

private:

   TFile *m_file;
   TChain *fChain;
   TString OutputFileName;
   TTree *out_tree;
   TTree *json_tree;

   /// Pax Information ///
   Event* events;

   /// JSON Information ///
   Json::Value json;
   Json::Reader reader;
  
   /// Leaky PMT Information ///
   vector<int> leaky_pmts;
   vector<double> ap_rate;

   /// General information ///
   TString Filename;
   Int_t EventId;
   Long64_t time;

   /// Interactions ///
   Int_t n_int;
   vector<Int_t> int_s1_id;
   vector<Int_t> int_s2_id;
   vector<Float_t> int_s1_area_correction;
   vector<Float_t> int_s1_pattern_fit;
   vector<Float_t> int_s1_saturation_correction;
   vector<Float_t> int_s1_spatial_correction;
   vector<Float_t> int_s2_area_correction;
   vector<Float_t> int_s2_lifetime_correction;
   vector<Float_t> int_x;
   vector<Float_t> int_y;
   vector<Float_t> int_z;
   vector<TString> int_xy_posrec_algorithm;
   vector<Float_t> int_xy_posrec_goodness_of_fit;
   vector<Float_t> int_xy_posrec_ndf;
   vector<Float_t> int_drift_time;

   /// S1 ///
   Int_t n_s1;
   vector<Int_t> s1_id;
   vector<Float_t> s1_area;
   vector<Float_t> s1_area_fraction_top;
   vector<Float_t> s1_area_fraction_leaky_pmts;
   vector<Float_t> s1_height;
   vector<Float_t> s1_top_hitpattern_spread;
   vector<Float_t> s1_bottom_hitpattern_spread;
   vector<Float_t> s1_center_time;
   vector<Float_t> s1_hit_time_mean;
   vector<Float_t> s1_hit_time_std;
   vector<Int_t> s1_n_hits;
   vector<Float_t> s1_hits_fraction_top;
   vector<Float_t> s1_hits_fraction_leaky_pmts;
   vector<Int_t> s1_n_contributing_channels;
   vector<Int_t> s1_n_contributing_channels_top;
   vector<Int_t> s1_n_saturated_channels;
   vector<Double_t> s1_range_area_decile_10p;
   vector<Double_t> s1_range_area_decile_50p;
   vector<Int_t> s1_n_noise_pulses;
   vector<vector<Float_t> > s1_x;
   vector<vector<Float_t> > s1_y;
   vector<vector<Float_t> > s1_goodness_of_fit;
   vector<vector<Float_t> > s1_ndf;
   vector<vector<TString> > s1_algorithm;

   /// S2 ///
   Int_t n_s2;
   vector<Int_t> s2_id;
   vector<Float_t> s2_area;
   vector<Float_t> s2_area_fraction_top;
   vector<Float_t> s2_area_fraction_leaky_pmts;
   vector<Float_t> s2_height;
   vector<Float_t> s2_top_hitpattern_spread;
   vector<Float_t> s2_bottom_hitpattern_spread;
   vector<Float_t> s2_center_time;
   vector<Float_t> s2_hit_time_mean;
   vector<Float_t> s2_hit_time_std;
   vector<Int_t> s2_n_hits;
   vector<Float_t> s2_hits_fraction_top;
   vector<Float_t> s2_hits_fraction_leaky_pmts;
   vector<Int_t> s2_n_contributing_channels;
   vector<Int_t> s2_n_contributing_channels_top;
   vector<Int_t> s2_n_saturated_channels;
   vector<Double_t> s2_range_area_decile_10p;
   vector<Double_t> s2_range_area_decile_50p;
   vector<Int_t> s2_n_noise_pulses;
   vector<vector<Float_t> > s2_x;
   vector<vector<Float_t> > s2_y;
   vector<vector<Float_t> > s2_goodness_of_fit;
   vector<vector<Float_t> > s2_ndf;
   vector<vector<TString> > s2_algorithm;

   /// Other Peaks (sorted by signal size) ///
   Int_t n_other_peaks;
   vector<Float_t> other_peaks_area;
   vector<Float_t> other_peaks_area_fraction_top;
   vector<Float_t> other_peaks_area_fraction_leaky_pmts;
   vector<Float_t> other_peaks_height;
   vector<Float_t> other_peaks_top_hitpattern_spread;
   vector<Float_t> other_peaks_bottom_hitpattern_spread;
   vector<Float_t> other_peaks_center_time;
   vector<Float_t> other_peaks_hit_time_mean;
   vector<Float_t> other_peaks_hit_time_std;
   vector<Int_t> other_peaks_n_hits;
   vector<Float_t> other_peaks_hits_fraction_top;
   vector<Float_t> other_peaks_hits_fraction_leaky_pmts;
   vector<Int_t> other_peaks_n_contributing_channels;
   vector<Int_t> other_peaks_n_contributing_channels_top;
   vector<Int_t> other_peaks_n_saturated_channels;
   vector<Double_t> other_peaks_range_area_decile_10p;
   vector<Double_t> other_peaks_range_area_decile_50p;
   vector<Int_t> other_peaks_n_noise_pulses;
   vector<vector<Float_t> > other_peaks_x;
   vector<vector<Float_t> > other_peaks_y;
   vector<vector<Float_t> > other_peaks_goodness_of_fit;
   vector<vector<Float_t> > other_peaks_ndf;
   vector<vector<TString> > other_peaks_algorithm;

};

#endif
