#ifndef gPlastNearlineSpectra_H
#define gPlastNearlineSpectra_H

#include "FairTask.h"
#include "TFolder.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TgPlastConfiguration.h"
#include <vector>

class TClonesArray;
class EventHeader;
class TH1F;
class TH2F;
class TFile;
class TFolder;
class TDirectory;

class gPlastNearlineSpectra : public FairTask
{
    public:
        gPlastNearlineSpectra();
        gPlastNearlineSpectra(const TString& name, Int_t verbose = 1);

        virtual ~gPlastNearlineSpectra();

        virtual InitStatus Init();

        virtual void Exec(Option_t* option);
        
        virtual void FinishEvent();

        virtual void FinishTask();

    
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

        //Folders and files
        TDirectory* dir_gplast;
        TDirectory* dir_gplast_slowToT;
        TDirectory* dir_gplast_fastToT;
        TDirectory* dir_gplast_hitpattern;
        TDirectory* dir_gplast_fast_vs_slow;
        TDirectory* dir_gplast_time_spectra;

        // Histograms
        std::vector<TH1F*> h1_gplast_slowToT;
        std::vector<TH1F*> h1_gplast_fastToT;
        TH1F* h1_gplast_hitpatterns;
        std::vector<TH1F*> h1_gplast_tamex_card_hitpattern;
        std::vector<TH1F*> h1_gplast_position_hitpattern;

        std::vector<TH2F*> h2_gplast_fastToT_vs_slowToT;
        std::vector<TH1F*> h1_gplast_time_spectra;

        // Detector Multiplicity
        TH1F* h1_gplast_multiplicity;

        int event_multiplicity;

    public:
        ClassDef(gPlastNearlineSpectra, 1)
};

#endif
