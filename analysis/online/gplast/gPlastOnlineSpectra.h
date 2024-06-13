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

        TString screenshot_path = "/u/despec/screenshots/";

        // Canvases
        TCanvas* c_gplast_slowToT;
        TCanvas* c_gplast_fastToT;
        TCanvas* c_gplast_hitpatterns;
        TCanvas* c_gplast_tamex_card_hitpattern;
        TCanvas* c_gplast_position_hitpattern;
        TCanvas* c_gplast_fast_v_slow;
        TCanvas* c_gplast_time_spectra;
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

        TFile* file_gplast_snapshot;

        // Histograms
        std::vector<TH1F*> h1_gplast_slowToT;
        std::vector<TH1F*> h1_gplast_fastToT;
        TH1F* h1_gplast_hitpatterns;
        std::vector<TH1F*> h1_gplast_position_hitpattern;
        std::vector<TH1F*> h1_gplast_tamex_card_hitpattern;

        std::vector<TH2F*> h2_gplast_fastToT_vs_slowToT;
        std::vector<TH1F*> h1_gplast_time_spectra;
        TH1F* h1_gplast_wr_time_diff;


        // Detector Multiplicity
        TH1F* h1_gplast_multiplicity;

        int event_multiplicity;
        int wr_t;
        int wr_prev;

    public:
        ClassDef(gPlastOnlineSpectra, 1)
};

#endif
