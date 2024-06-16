#ifndef gPlastOnlineSpectra_H
#define gPlastOnlineSpectra_H

#include "FairTask.h"
#include "TFolder.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TgPlastConfiguration.h"
#include <vector>

class TClonesArray;
class EventHeader;
class TCanvas;
class TH1F;
class TH2F;
class TFile;
class TFolder;
class TDirectory;

class gPlastOnlineSpectra : public FairTask
{
    public:
        gPlastOnlineSpectra();
        gPlastOnlineSpectra(const TString& name, Int_t verbose = 1);
        
        virtual ~gPlastOnlineSpectra();

        virtual void SetParContainers();

        virtual InitStatus Init();

        virtual void Exec(Option_t* option);
        
        virtual void FinishEvent();

        virtual void FinishTask();

        virtual void Reset_Histo();

        virtual void Snapshot_Histo();

        // range setters

    
    private:
        TClonesArray* fHitgPlastTwinpeaks;

        TgPlastConfiguration const* gplast_conf;
        std::map<std::pair<int, int>, int> gplast_map;
        std::map<int, std::pair<TgPlastConfiguration::Position, int>> gplast_pos;

        EventHeader* header;
        Int_t fNEvents;
        int total_time_microsecs = 0;

        int nDetectors;
        int nTamexBoards;
        int nPos;
        int nSignalsPlastic;

        TString screenshot_path = "/u/despec/screenshots/";

        // Canvases
        TCanvas* c_gplast_slowToT;
        TCanvas* c_gplast_fastToT;
        TCanvas* c_gplast_hitpatterns;
        TCanvas* c_gplast_tamex_card_hitpattern;
        TCanvas* c_gplast_position_hitpattern;
        TCanvas* c_gplast_position_hitpattern_spill[2];
        TCanvas* c_gplast_coincidences_matrix;
        TCanvas* c_gplast_time_differences;
        TCanvas* c_gplast_time_differences_position;
        TCanvas* c_gplast_time_differences_xy;
        TCanvas* c_gplast_fast_v_slow;
        TCanvas* c_gplast_time_spectra_sci41;
        TCanvas* c_gplast_time_vs_position;
        TCanvas* c_gplast_time_vs_position_sci41;
        TCanvas* c_gplast_tof_sci41;
        TCanvas* c_gplast_tof_sci42;
        TCanvas* c_gplast_tof_sci43;
        TCanvas* c_gplast_total_slowTot;
        TCanvas* c_gplast_total_slowTot_sci41;
        //TCanvas* c_gplast_total_slowTot_position;
        TCanvas* c_gplast_multiplicity;
        TCanvas* c_gplast_channel_multiplicity;
        TCanvas* c_gplast_wr_time_diff;
        TCanvas* c_gplast_snapshot;

        //Folders and files
        TFolder* histograms;
        TDirectory* dir_gplast;
        TDirectory* dir_gplast_slowToT;
        TDirectory* dir_gplast_fastToT;
        TDirectory* dir_gplast_hitpattern;
        TDirectory* dir_gplast_fast_vs_slow;
        TDirectory* dir_gplast_slowToT_total;
        TDirectory* dir_gplast_coincidences;
        TDirectory* dir_gplast_time_spectra;
        TDirectory* dir_gplast_ToF_sci;

        TFile* file_gplast_snapshot;

        // Histograms
        std::vector<TH1*> h1_gplast_slowToT;
        std::vector<TH1*> h1_gplast_fastToT;
        TH1* h1_gplast_hitpatterns;
        std::vector<TH1*> h1_gplast_position_hitpattern;
        std::vector<TH1*> h1_gplast_position_hitpattern_spill[2];
        std::vector<TH1*> h1_gplast_tamex_card_hitpattern;

        std::vector<TH2*> h2_gplast_fastToT_vs_slowToT;
        std::vector<TH1*> h1_gplast_time_spectra_sci41;
        std::vector<TH2*> h2_gplast_time_vs_position;
        std::vector<TH2*> h2_gplast_time_sci41_vs_position;
        TH1* h1_gplast_wr_time_diff;

        TH1* h1_gplast_tof_sci41;
        TH1* h1_gplast_tof_sci41_veto;
        TH1* h1_gplast_tof_sci42;
        TH1* h1_gplast_tof_sci42_veto;
        TH1* h1_gplast_tof_sci43;

        TH1* h1_gplast_slowToT_total;
        TH1* h1_gplast_slowToT_total_veto;
        TH2* h2_gplast_slowToT_total_sci41;
        //std::vector<TH1F*> h1_gplast_slowToT_position;

        // Detector Multiplicity
        TH1* h1_gplast_multiplicity;
        TH2* h2_gplast_coincidence_matrix;
        TH1* h1_gplast_coinc_time_diff;
        std::vector<TH1*> h1_gplast_coinc_time_diff_position;
        TH1* h1_gplast_coinc_time_diff_vs_x;
        TH1* h1_gplast_coinc_time_diff_vs_y;

        int event_multiplicity;
        int wr_t;
        int wr_prev;

    public:
        ClassDef(gPlastOnlineSpectra, 1)
};

#endif
