// FairRoot
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"

// c4
#include "gPlastOnlineSpectra.h"
#include "EventHeader.h"
#include "gPlastTwinpeaksCalData.h"
#include "TgPlastConfiguration.h"

#include "c4Logger.h"

#include "TCanvas.h"
#include "TClonesArray.h"
#include "THttpServer.h"
#include "TMath.h"
#include "TRandom.h"
#include "TFile.h"
#include "TROOT.h"
#include <chrono>
#include <TKey.h>

gPlastOnlineSpectra::gPlastOnlineSpectra() : gPlastOnlineSpectra("gPlastOnlineSpectra")
{
}

gPlastOnlineSpectra::gPlastOnlineSpectra(const TString& name, Int_t verbose)
    : FairTask(name, verbose)
    , fHitgPlastTwinpeaks(NULL)
    , fNEvents(0)
    , header(nullptr)
{
}

gPlastOnlineSpectra::~gPlastOnlineSpectra()
{
    c4LOG(info, "");
    if (fHitgPlastTwinpeaks) delete fHitgPlastTwinpeaks;
}

void gPlastOnlineSpectra::SetParContainers()
{
    FairRuntimeDb *rtdb = FairRuntimeDb::instance();
    c4LOG_IF(fatal, NULL == rtdb, "FairRuntimeDb not found.");
}

InitStatus gPlastOnlineSpectra::Init()
{
    FairRootManager* mgr = FairRootManager::Instance();
    c4LOG_IF(fatal, NULL == mgr, "FairRootManager not found");

    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);

    header = (EventHeader*)mgr->GetObject("EventHeader.");
    c4LOG_IF(error, !header, "Branch EventHeader not found");

    fHitgPlastTwinpeaks = (TClonesArray*)mgr->GetObject("gPlastTwinpeaksCalData");
    c4LOG_IF(fatal, !fHitgPlastTwinpeaks, "Branch gPlastTwinpeaksCalData not found");

    histograms = (TFolder*)mgr->GetObject("Histograms");

    TDirectory::TContext ctx(nullptr);

    dir_gplast = new TDirectory("gPlast","gPlast", "", 0);
    // mgr->Register("gPlast", "gPlast Directory", dir_gplast, false); // allow other tasks to access directory.
    histograms->Add(dir_gplast);

    dir_gplast_slowToT = dir_gplast->mkdir("Slow ToT");
    dir_gplast_fastToT = dir_gplast->mkdir("Fast ToT");
    dir_gplast_hitpattern = dir_gplast->mkdir("Hit Pattern");
    dir_gplast_fast_vs_slow = dir_gplast->mkdir("Fast Vs. Slow");

    // gPlast Configuration
    gplast_conf = TgPlastConfiguration::GetInstance();
    gplast_map = gplast_conf->Mapping();
    gplast_pos = gplast_conf->Positioning();
    nDetectors = gplast_conf->NDetectors();
    nTamexBoards = gplast_conf->NTamexBoards();
    nPos = TgPlastConfiguration::NB_POS-1; // discard EXTRA
    
    // Setting histogram sizes
    h1_gplast_slowToT.resize(nDetectors+1); // index from 1 
    h1_gplast_fastToT.resize(nDetectors+1);
    h1_gplast_tamex_card_hitpattern.resize(nTamexBoards);
    h1_gplast_position_hitpattern.resize(nPos);
    h2_gplast_fastToT_vs_slowToT.resize(nDetectors+1);

    // Slow ToT
    dir_gplast_slowToT->cd(); 
    c_gplast_slowToT  = new TCanvas("c_gplast_slowToT","slow ToT gPlast spectra",1200,800);
    c_gplast_slowToT->Divide(5,(nDetectors%5==0) ? (nDetectors/5) : (nDetectors/5 + 1));
    for (int ihist = 1; ihist <= nDetectors; ihist++){
        c_gplast_slowToT->cd(ihist);
        h1_gplast_slowToT[ihist] = new TH1F(Form("h1_gplast_slowToT_%d",ihist),Form("gPlastic Slow ToT %d",ihist),10000,0,3.5e3);
        h1_gplast_slowToT[ihist]->GetXaxis()->SetTitle("ToT (ns)");
        h1_gplast_slowToT[ihist]->Draw();
    }
    c_gplast_slowToT->cd(0);
    dir_gplast_slowToT->Append(c_gplast_slowToT);

    // Fast ToT
    dir_gplast_fastToT->cd();
    c_gplast_fastToT  = new TCanvas("c_gplast_fastToT","Fast ToT gPlast spectra",1200,800);
    c_gplast_fastToT->Divide(5,(nDetectors%5==0) ? (nDetectors/5) : (nDetectors/5 + 1));
    for (int ihist = 1; ihist <= nDetectors; ihist++)
    {
        c_gplast_fastToT->cd(ihist);
        h1_gplast_fastToT[ihist] = new TH1F(Form("h1_gplast_fastToT_%d",ihist),Form("gPlast Fast ToT %d",ihist),10000,0,3.5e3);
        h1_gplast_fastToT[ihist]->GetXaxis()->SetTitle("ToT (ns)");
        h1_gplast_fastToT[ihist]->Draw();   
    }
    c_gplast_fastToT->cd(0);
    dir_gplast_fastToT->Append(c_gplast_fastToT);

    // Hit Pattern
    dir_gplast_hitpattern->cd();
    c_gplast_hitpatterns = new TCanvas("c_gplast_hitpatterns","gPlast Hit Pattern",1200,800);
    h1_gplast_hitpatterns = new TH1F("h1_gplast_hitpattern","Detector hit patterns",nDetectors,1,nDetectors+1);
    h1_gplast_hitpatterns->GetXaxis()->SetTitle("Detector ID");
    h1_gplast_hitpatterns->GetYaxis()->SetTitle("Counts");
    h1_gplast_hitpatterns->Draw();
    c_gplast_hitpatterns->cd();
    dir_gplast_hitpattern->Append(c_gplast_hitpatterns);

    // tamex card hit pattern
    c_gplast_tamex_card_hitpattern  = new TCanvas("c_gplast_tamex_card_hitpattern","gPlast Tamex Card Hit Pattern",1200,800);
    c_gplast_tamex_card_hitpattern->Divide(5,(nTamexBoards%5==0) ? (nTamexBoards/5) : (nTamexBoards/5 + 1));
    for (int ihist = 0; ihist < nTamexBoards; ihist++){
        c_gplast_tamex_card_hitpattern->cd(ihist+1);
        h1_gplast_tamex_card_hitpattern[ihist] = new TH1F(Form("h1_gplast_tamex_card_hitpattern_%d",ihist),Form("gPlast Hit Pattern Tamex Card %d",ihist),16,1,17); // hard-coded : number of channels per TAMEX
        h1_gplast_tamex_card_hitpattern[ihist]->GetXaxis()->SetTitle("Channel");
        h1_gplast_tamex_card_hitpattern[ihist]->GetYaxis()->SetTitle("Counts");
        h1_gplast_tamex_card_hitpattern[ihist]->Draw();
    }
    c_gplast_tamex_card_hitpattern->cd(0);
    dir_gplast_hitpattern->Append(c_gplast_tamex_card_hitpattern);

    // position hit pattern
    c_gplast_position_hitpattern  = new TCanvas("c_gplast_position_hitpattern","gPlast Hit Pattern vs Position",1200,800);
    c_gplast_position_hitpattern->Divide(5,(nPos%5==0) ? (nPos/5) : (nPos/5 + 1));
    for (int ihist = 0; ihist < nPos; ihist++){
        c_gplast_position_hitpattern->cd(ihist+1);
        auto str = TgPlastConfiguration::PositionToString(static_cast<TgPlastConfiguration::Position>(ihist+1)); //EXTRA = 0
        h1_gplast_position_hitpattern[ihist] = new TH1F(Form("h1_gplast_hitpattern_%s",str.data()),Form("gPlast Hit Pattern %s",str.data()),16,1,17); // hard-coded : number of channels per position
        h1_gplast_position_hitpattern[ihist]->GetXaxis()->SetTitle("Channel");
        h1_gplast_position_hitpattern[ihist]->GetYaxis()->SetTitle("Counts");
        h1_gplast_position_hitpattern[ihist]->Draw();
    }
    c_gplast_position_hitpattern->cd(0);
    dir_gplast_hitpattern->Append(c_gplast_position_hitpattern);

    dir_gplast_fast_vs_slow->cd();
    c_gplast_fast_v_slow  = new TCanvas("c_gplast_fast_v_slow","fast vs slow ToT gplast spectra",1200,800);
    c_gplast_fast_v_slow->Divide(5,(nDetectors%5==0) ? (nDetectors/5) : (nDetectors/5 + 1));
    for (int ihist = 1; ihist <= nDetectors; ihist++){
        c_gplast_fast_v_slow->cd(ihist);
        h2_gplast_fastToT_vs_slowToT[ihist] = new TH2F(Form("h1_gplast_fast_v_slow_%d",ihist),Form("gplast fast vs. slow detector %d",ihist),1000,0,3.5e3,1000,0,3.5e3);
        h2_gplast_fastToT_vs_slowToT[ihist]->GetXaxis()->SetTitle("Slow ToT (ns)");
        h2_gplast_fastToT_vs_slowToT[ihist]->GetYaxis()->SetTitle("Fast ToT (ns)");
        h2_gplast_fastToT_vs_slowToT[ihist]->SetOption("COLZ");
    }
    c_gplast_fast_v_slow->cd(0);
    dir_gplast_fast_vs_slow->Append(c_gplast_fast_v_slow);

    dir_gplast_hitpattern->cd();

    c_gplast_multiplicity = new TCanvas("c_gplast_multiplicity","gPlast multiplicity spectrum",1200,800);
    h1_gplast_multiplicity = new TH1F("h1_gplast_multiplicity","gPlast multiplicity (>=2)",nDetectors,1,nDetectors+1);
    h1_gplast_multiplicity->GetXaxis()->SetTitle("Channel Multiplicity");
    h1_gplast_multiplicity->GetYaxis()->SetTitle("Counts");
    h1_gplast_multiplicity->Draw();
    c_gplast_multiplicity->SetLogy();
    c_gplast_multiplicity->cd();
    dir_gplast_hitpattern->Append(c_gplast_multiplicity);

    c_gplast_wr_time_diff  = new TCanvas("c_gplast_wr_time_diff","gPlast WR time difference",1200,800);
    h1_gplast_wr_time_diff = new TH1F("h1_gplast_wr_time_diff","gPlast WR time difference",100,-1e2,80e3);
    h1_gplast_wr_time_diff->GetXaxis()->SetTitle("White Rabbit Event Time Difference (ns)");
    h1_gplast_wr_time_diff->Draw();
    c_gplast_wr_time_diff->cd();
    dir_gplast_hitpattern->Append(c_gplast_wr_time_diff);

    dir_gplast->cd();

    run->GetHttpServer()->RegisterCommand("Reset_gPlast_Histo", Form("/Objects/%s/->Reset_Histo()", GetName()));
    run->GetHttpServer()->RegisterCommand("Snapshot_gPlast_Histo", Form("/Objects/%s/->Snapshot_Histo()", GetName()));

    return kSUCCESS;
}

void gPlastOnlineSpectra::Reset_Histo()
{
    c4LOG(info, "Resetting gPlast histograms.");
    for (int ihist = 1; ihist<=nDetectors; ihist++) h1_gplast_slowToT[ihist]->Reset();
    for (int ihist = 1; ihist<=nDetectors; ihist++) h1_gplast_fastToT[ihist]->Reset();
    for (int ihist = 1; ihist<=nDetectors; ihist++) h2_gplast_fastToT_vs_slowToT[ihist]->Reset();
    h1_gplast_hitpatterns->Reset();
    h1_gplast_multiplicity->Reset();
    h1_gplast_wr_time_diff->Reset();
    for (int ihist = 0; ihist < nTamexBoards; ihist++) h1_gplast_tamex_card_hitpattern[ihist]->Reset();
    for (int ihist = 0; ihist < nPos; ihist++) h1_gplast_position_hitpattern[ihist]->Reset();
    
    c4LOG(info, "gPlast histograms reset.");
}

// make a date and time stamped folder with pngs of the histograms and .root file and save them
void gPlastOnlineSpectra::Snapshot_Histo()
{
    c4LOG(info, "Snapshotting gPlast histograms.");
    // date and time stamp

    time_t now = time(0);
    tm *ltm = localtime(&now);
    TDirectory *tmp = gDirectory;
    tmp = dir_gplast;
    // make a folder with the date and time
    const char* snapshot_dir = Form("gPlast_Snapshots_%d_%d_%d_%d_%d_%d", 1900 + ltm->tm_year, 1 + ltm->tm_mon, ltm->tm_mday, ltm->tm_hour, ltm->tm_min, ltm->tm_sec);
    gSystem->cd(screenshot_path);
    gSystem->mkdir(snapshot_dir);
    gSystem->cd(snapshot_dir);

    // save histograms
    c_gplast_snapshot = new TCanvas("c","c",1200,800);
    for (int ihist = 1; ihist<=nDetectors; ihist++) {
        h1_gplast_slowToT[ihist]->Draw();
        c_gplast_snapshot->SaveAs(Form("Slow_ToT_%d.png", ihist));
        c_gplast_snapshot->Clear();
        h1_gplast_fastToT[ihist]->Draw();
        c_gplast_snapshot->SaveAs(Form("Fast_ToT_%d.png", ihist));
        c_gplast_snapshot->Clear();
        h2_gplast_fastToT_vs_slowToT[ihist]->Draw("COLZ");
        c_gplast_snapshot->SaveAs(Form("Fast_vs._Slow_ToT_%d.png", ihist));
        c_gplast_snapshot->Clear();
    }
    
    for (int ihist = 0; ihist < nTamexBoards; ihist++){
        h1_gplast_tamex_card_hitpattern[ihist]->Draw();
        c_gplast_snapshot->SaveAs(Form("gPlast_Tamex_Card_HitPattern_%d.png", ihist));
        c_gplast_snapshot->Clear();
    }

    for (int ihist = 0; ihist < nPos; ihist++){
        h1_gplast_position_hitpattern[ihist]->Draw();
        c_gplast_snapshot->SaveAs(Form("gPlast_Position_HitPattern_%d.png", ihist));
        c_gplast_snapshot->Clear();
    }

    // save hit patterns
    c_gplast_hitpatterns->SaveAs("gPlast_HitPatterns.png");
    c_gplast_multiplicity->SaveAs("gPlast_Multiplicity.png");
    c_gplast_wr_time_diff->SaveAs("gPlast_WR_Time_Difference.png");

    gSystem->cd("..");
    c4LOG(info, "gPlast snapshot saved in:" << screenshot_path + snapshot_dir);
}

void gPlastOnlineSpectra::Exec(Option_t* option)
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
            
            // event_multiplicity_downstream ++;
            double slowToT = hit->Get_slow_ToT();
            double fastToT = hit->Get_fast_ToT();
            double fast_lead_time = hit->Get_fast_lead_time();
            
            int detector_id = hit->Get_detector_id();

            if (detector_id > nDetectors || detector_id < 0) continue;

            // Fast and Slow Tot spectra
            h1_gplast_slowToT[detector_id]->Fill(slowToT);
            h1_gplast_fastToT[detector_id]->Fill(fastToT);

            h2_gplast_fastToT_vs_slowToT[detector_id]->Fill(slowToT,fastToT);

            // White Rabbit time difference
            gPlastTwinpeaksCalData* wr_hit = (gPlastTwinpeaksCalData*)fHitgPlastTwinpeaks->At(0);
            wr_t = wr_hit->Get_wr_t();
            if (wr_t != wr_prev)
            {
                int dt = wr_t - wr_prev;
                h1_gplast_wr_time_diff->Fill(dt);
            } 
            wr_prev = wr_t;

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
        if(event_multiplicity >= 2) h1_gplast_multiplicity->Fill(event_multiplicity);
    }
    fNEvents += 1;

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    total_time_microsecs += duration.count();
}

void gPlastOnlineSpectra::FinishEvent()
{
    if (fHitgPlastTwinpeaks)
    {
        fHitgPlastTwinpeaks->Clear();
    }
}

void gPlastOnlineSpectra::FinishTask()
{
    if(fNEvents == 0)
    { 
        c4LOG(warn, "No events found, not saving histograms!");
        return;
    }
   
    c4LOG(info, "Average execution time: " << (double)total_time_microsecs/fNEvents << " microseconds.");
    
}

ClassImp(gPlastOnlineSpectra)
