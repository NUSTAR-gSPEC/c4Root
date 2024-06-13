// FairRoot
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "FairSink.h"

// c4
#include "gPlastNearlineSpectra.h"
#include "EventHeader.h"
#include "gPlastTwinpeaksCalData.h"
#include "TgPlastConfiguration.h"

#include "c4Logger.h"

#include "TClonesArray.h"
#include "TMath.h"
#include "TRandom.h"
#include "TFile.h"
#include "TROOT.h"
#include <chrono>

gPlastNearlineSpectra::gPlastNearlineSpectra() : gPlastNearlineSpectra("gPlastNearlineSpectra")
{
}

gPlastNearlineSpectra::gPlastNearlineSpectra(const TString& name, Int_t verbose)
    : FairTask(name, verbose)
    , fHitgPlastTwinpeaks(NULL)
    , fNEvents(0)
    , header(nullptr)
{
}

gPlastNearlineSpectra::~gPlastNearlineSpectra()
{
    c4LOG(info, "");
    if (fHitgPlastTwinpeaks)
        delete fHitgPlastTwinpeaks;
}

InitStatus gPlastNearlineSpectra::Init()
{
    // set batch mode // ok but why here specifically?
    gROOT->SetBatch(kTRUE);

    FairRootManager* mgr = FairRootManager::Instance();
    c4LOG_IF(fatal, NULL == mgr, "FairRootManager not found");

    FairRunAna* run = FairRunAna::Instance();

    header = (EventHeader*)mgr->GetObject("EventHeader.");
    c4LOG_IF(error, !header, "Branch EventHeader. not found");

    fHitgPlastTwinpeaks = (TClonesArray*)mgr->GetObject("gPlastTwinpeaksCalData");
    c4LOG_IF(fatal, !fHitgPlastTwinpeaks, "Branch gPlastTwinpeaksCalData not found!");

    TDirectory* tmp = gDirectory;
    FairRootManager::Instance()->GetOutFile()->cd();
    dir_gplast = gDirectory->mkdir("gPlast");
    gDirectory->cd("gPlast");

    dir_gplast_slowToT = gDirectory->mkdir("Slow ToT");
    dir_gplast_fastToT = gDirectory->mkdir("Fast ToT");
    dir_gplast_hitpattern = gDirectory->mkdir("Hit Pattern");
    dir_gplast_fast_vs_slow = gDirectory->mkdir("Fast Vs. Slow");
    dir_gplast_time_spectra = gDirectory->mkdir("Time Spectra");

    // gPlast Configuration
    gplast_conf = TgPlastConfiguration::GetInstance();
    gplast_map = gplast_conf->Mapping();
    gplast_pos = gplast_conf->Positioning();
    nDetectors = gplast_conf->NDetectors();
    nTamexBoards = gplast_conf->NTamexBoards();
    constexpr int nPos = TgPlastConfiguration::NB_POS-1; // discard EXTRA

    // Setting histogram sizes
    h1_gplast_slowToT.resize(nDetectors+1); // index from 1 
    h1_gplast_fastToT.resize(nDetectors+1);
    h1_gplast_tamex_card_hitpattern.resize(nTamexBoards);
    h1_gplast_position_hitpattern.resize(nPos);
    h2_gplast_fastToT_vs_slowToT.resize(nDetectors+1);
    h1_gplast_time_spectra.resize(nDetectors+1);

    // Slow ToT
    dir_gplast_slowToT->cd();
    for (int ihist = 1; ihist <= nDetectors; ihist++)
    {
        h1_gplast_slowToT[ihist] = new TH1F(Form("h1_gplast_slowToT_%d",ihist),Form("gPlastic Slow ToT %d",ihist),10000,0,3.5e3);
        h1_gplast_slowToT[ihist]->GetXaxis()->SetTitle("ToT (ns)");
    }

    // Fast ToT
    dir_gplast_fastToT->cd();
    for (int ihist = 1; ihist <= nDetectors; ihist++)
    {
        h1_gplast_fastToT[ihist] = new TH1F(Form("h1_gplast_fastToT_%d",ihist),Form("gPlastic Fast ToT %d",ihist),10000,0,3.5e3);
        h1_gplast_fastToT[ihist]->GetXaxis()->SetTitle("ToT (ns)"); 
    }

    // Hit Pattern
    dir_gplast_hitpattern->cd();

    h1_gplast_hitpatterns = new TH1F("h1_gplast_hitpattern","gPlast hit patterns",64,1,65);
    h1_gplast_hitpatterns->GetXaxis()->SetTitle("Detector ID");
    h1_gplast_hitpatterns->GetYaxis()->SetTitle("Counts");

    // tamex card hit pattern
    for (int ihist = 0; ihist < nTamexBoards; ihist++)
    {
        h1_gplast_tamex_card_hitpattern[ihist] = new TH1F(Form("h1_gplast_tamex_card_hitpattern_%d",ihist),Form("gPlastic Tamex Card Hit Pattern %d",ihist),16,1,17);
        h1_gplast_tamex_card_hitpattern[ihist]->GetXaxis()->SetTitle("Channel");
        h1_gplast_tamex_card_hitpattern[ihist]->GetYaxis()->SetTitle("Counts");
    }

    // position hit pattern
    for (int ihist = 0; ihist < nPos; ihist++){
        auto str = TgPlastConfiguration::PositionToString(static_cast<TgPlastConfiguration::Position>(ihist+1)); //EXTRA = 0
        h1_gplast_position_hitpattern[ihist] = new TH1F(Form("h1_gplast_hitpattern_%s",str.data()),Form("gPlast Hit Pattern %s",str.data()),16,1,17); // hard-coded : number of channels per position
        h1_gplast_position_hitpattern[ihist]->GetXaxis()->SetTitle("Channel");
        h1_gplast_position_hitpattern[ihist]->GetYaxis()->SetTitle("Counts");
    }

    // Time spectra
    dir_gplast_time_spectra->cd();
    for (int ihist = 1; ihist <= nDetectors; ihist++)
    {
        h1_gplast_time_spectra[ihist] = new TH1F(Form("h1_gplast_time_spectra_%d",ihist),Form("gPlast Time spectrum detector %d",ihist),10000,0,6e10);
        h1_gplast_time_spectra[ihist]->GetXaxis()->SetTitle("Time (ns)");
    }

    dir_gplast_fast_vs_slow->cd();
    for (int ihist = 1; ihist <= nDetectors; ihist++)
    {
        h2_gplast_fastToT_vs_slowToT[ihist] = new TH2F(Form("h1_gplast_fast_v_slow_%d",ihist),Form("gplast fast vs. slow detector %d",ihist),1000,0,3.5e3,1000,0,3.5e3);
        h2_gplast_fastToT_vs_slowToT[ihist]->GetXaxis()->SetTitle("Slow ToT (ns)");
        h2_gplast_fastToT_vs_slowToT[ihist]->GetYaxis()->SetTitle("Fast ToT (ns)");
        h2_gplast_fastToT_vs_slowToT[ihist]->SetOption("COLZ");
    }

    dir_gplast_hitpattern->cd();

    h1_gplast_multiplicity = new TH1F("Multiplicity","gPlast multiplicity",128,1,129);
    h1_gplast_multiplicity->GetXaxis()->SetTitle("Event Multilplicity");
    h1_gplast_multiplicity->GetYaxis()->SetTitle("Counts");
    // log y?

    dir_gplast->cd();
    gDirectory = tmp;

    return kSUCCESS;
}

void gPlastNearlineSpectra::Exec(Option_t* option)
{   
    auto start = std::chrono::high_resolution_clock::now();

    if (fHitgPlastTwinpeaks && fHitgPlastTwinpeaks->GetEntriesFast() > 0)
    {
        event_multiplicity = 0;
        Int_t nHits = fHitgPlastTwinpeaks->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {   
            gPlastTwinpeaksCalData* hit = (gPlastTwinpeaksCalData*)fHitgPlastTwinpeaks->At(ihit);
            if (!hit) continue;
            
            double slowToT = hit->Get_slow_ToT();
            double fastToT = hit->Get_fast_ToT();
            double fast_lead_time = hit->Get_fast_lead_time();
            
            int detector_id = hit->Get_detector_id();

            if (detector_id > nDetectors || detector_id < 0) continue;

            // Lead and Trail spectra -- empty for now
            h1_gplast_time_spectra[detector_id]->Fill(fast_lead_time);

            // Fast and Slow Tot spectra
            h1_gplast_slowToT[detector_id]->Fill(slowToT);
            h1_gplast_fastToT[detector_id]->Fill(fastToT);

            h2_gplast_fastToT_vs_slowToT[detector_id]->Fill(slowToT,fastToT);

            // Hit pattern spectra
            // the hit pattern spectrum is generated by filling the histogram with the detector ID of the hit
            auto pos = gplast_pos.at(detector_id);
            if (pos.first != TgPlastConfiguration::EXTRA){
                h1_gplast_hitpatterns->Fill(detector_id);
                h1_gplast_position_hitpattern[pos.first-1]->Fill(pos.second);
                event_multiplicity++;
            }

            for (const auto& entry : gplast_map)
            {
                if (entry.second == detector_id)
                {
                    for (int i = 0; i < nTamexBoards; i++)
                    {
                        if (entry.first.first == i)
                        {
                            h1_gplast_tamex_card_hitpattern[i]->Fill(entry.first.second);
                        }
                    }
                }

            }


        }
        h1_gplast_multiplicity->Fill(event_multiplicity);
    }
    fNEvents += 1;

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    total_time_microsecs += duration.count();
}

void gPlastNearlineSpectra::FinishEvent()
{
    if (fHitgPlastTwinpeaks)
    {
        fHitgPlastTwinpeaks->Clear();
    }
}

void gPlastNearlineSpectra::FinishTask()
{
    if(fNEvents == 0)
    { 
        c4LOG(warn, "No events found, not saving histograms!");
        return;
    }

    TDirectory* tmp = gDirectory;
    FairRootManager::Instance()->GetOutFile()->cd();
    dir_gplast->Write();
    gDirectory = tmp;

    c4LOG(info, "Average execution time: " << (double)total_time_microsecs/fNEvents << " microseconds.");

}

ClassImp(gPlastNearlineSpectra)
