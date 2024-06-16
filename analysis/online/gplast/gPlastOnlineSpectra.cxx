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

#include "AnalysisTools.h"
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
    dir_gplast_slowToT_total = dir_gplast->mkdir("Total Slow ToT");
    dir_gplast_coincidences = dir_gplast->mkdir("Coincidences");
    dir_gplast_time_spectra = dir_gplast->mkdir("Timing");
    dir_gplast_ToF_sci = dir_gplast->mkdir("ToF vs FRS");

    // gPlast Configuration
    gplast_conf = TgPlastConfiguration::GetInstance();
    gplast_map = gplast_conf->Mapping();
    gplast_pos = gplast_conf->Positioning();
    nDetectors = gplast_conf->NDetectors();
    nTamexBoards = gplast_conf->NTamexBoards();
    nPos = TgPlastConfiguration::NB_POS-1; // discard EXTRA
    nSignalsPlastic = nDetectors - gplast_conf->NExtraSignals();
    
    // Setting histogram sizes
    h1_gplast_slowToT.resize(nDetectors+1); // index from 1 
    h1_gplast_fastToT.resize(nDetectors+1);
    h1_gplast_tamex_card_hitpattern.resize(nTamexBoards);
    h1_gplast_tamex_card_hitpattern.resize(nTamexBoards);
    h1_gplast_position_hitpattern.resize(nPos);
    h1_gplast_position_hitpattern_spill[0].resize(nPos);
    h1_gplast_position_hitpattern_spill[1].resize(nPos);
    h2_gplast_fastToT_vs_slowToT.resize(nDetectors+1);
    h1_gplast_time_spectra_sci41.resize(nDetectors+1);
    h2_gplast_time_vs_position.resize(nPos);
    h1_gplast_coinc_time_diff_position.resize(nPos);
    h2_gplast_time_sci41_vs_position.resize(nPos);

    // ======== Slow ToT ======== //
    c_gplast_slowToT  = new TCanvas("c_gplast_slowToT","slow ToT gPlast spectra",1200,800);
    c_gplast_slowToT->Divide(5,(nDetectors%5==0) ? (nDetectors/5) : (nDetectors/5 + 1));
    for (int ihist = 1; ihist <= nDetectors; ihist++){
        c_gplast_slowToT->cd(ihist);
        h1_gplast_slowToT[ihist] = MakeTH1(dir_gplast_slowToT, "F", Form("h1_gplast_slowToT_%d",ihist),Form("gPlast Slow ToT %d",ihist), 1e4, 0, 3.5e3, "ToT [ns]", kSpring, kBlue+2);
        h1_gplast_slowToT[ihist]->Draw();
    }
    c_gplast_slowToT->cd(0);
    dir_gplast_slowToT->Append(c_gplast_slowToT);

    // ======== Fast ToT ======== //
    c_gplast_fastToT  = new TCanvas("c_gplast_fastToT","Fast ToT gPlast spectra",1200,800);
    c_gplast_fastToT->Divide(5,(nDetectors%5==0) ? (nDetectors/5) : (nDetectors/5 + 1));
    for (int ihist = 1; ihist <= nDetectors; ihist++)
    {
        c_gplast_fastToT->cd(ihist);
        h1_gplast_fastToT[ihist] = MakeTH1(dir_gplast_fastToT, "F", Form("h1_gplast_fastToT_%d",ihist),Form("gPlast Fast ToT %d",ihist),10000,0,3.5e3, "ToT [ns]", kSpring, kBlue+2);
        h1_gplast_fastToT[ihist]->Draw();
    }
    c_gplast_fastToT->cd(0);
    dir_gplast_fastToT->Append(c_gplast_fastToT);

    // ======== Fast vs Slow ToT ======== //
    c_gplast_fast_v_slow  = new TCanvas("c_gplast_fast_v_slow","fast vs slow ToT gplast spectra",1200,800);
    c_gplast_fast_v_slow->Divide(5,(nDetectors%5==0) ? (nDetectors/5) : (nDetectors/5 + 1));
    for (int ihist = 1; ihist <= nDetectors; ihist++){
        c_gplast_fast_v_slow->cd(ihist);
        h2_gplast_fastToT_vs_slowToT[ihist] = MakeTH2(dir_gplast_fast_vs_slow, "F", Form("h2_gplast_fast_v_slow_%d",ihist),Form("gplast fast vs. slow detector %d",ihist),1000,0,3.5e3,1000,0,3.5e3, "Slow ToT [ns]", "Fast ToT [ns]");
        h2_gplast_fastToT_vs_slowToT[ihist]->Draw();
    }
    c_gplast_fast_v_slow->cd(0);
    dir_gplast_fast_vs_slow->Append(c_gplast_fast_v_slow);

    // ======== Total Slow ToT ======== //
    c_gplast_total_slowTot  = new TCanvas("c_gplast_total_slowTot","gPlast total slow ToT",1200,800);
    c_gplast_total_slowTot->Divide(2);
    c_gplast_total_slowTot->cd(1);
    h1_gplast_slowToT_total = MakeTH1(dir_gplast_slowToT_total, "F", "h1_gplast_slowToT_total","gPlast total slow ToT - no veto",10000,0,3.5e3, "Total slow ToT [ns]", kSpring+1, kBlue+2); // hard-coded
    h1_gplast_slowToT_total->Draw();
    c_gplast_total_slowTot->cd(2);
    h1_gplast_slowToT_total_veto = MakeTH1(dir_gplast_slowToT_total, "F", "h1_gplast_slowToT_total_veto","gPlast total slow ToT - vetoed",10000,0,3.5e3, "Total slow ToT [ns]", kSpring+1, kBlue+2); // hard-coded
    h1_gplast_slowToT_total_veto->Draw();
    c_gplast_total_slowTot->cd();
    dir_gplast_slowToT_total->Append(c_gplast_total_slowTot);

    c_gplast_total_slowTot_sci41 = new TCanvas("c_gplast_total_slowTot_sci41","gPlast total slow ToT vs sci41",1200,800);
    h2_gplast_slowToT_total_sci41 = MakeTH2(dir_gplast_slowToT_total, "F", "h2_gplast_slowToT_total_sci41","gPlast total slow ToT vs sci41", 10000, 0., 3.5e3, 2000, 0, 200, "Total slow ToT [ns]", "Time difference with sci41 [ns]");
    h2_gplast_slowToT_total_sci41->Draw();
    c_gplast_total_slowTot_sci41->cd();
    dir_gplast_slowToT_total->Append(c_gplast_total_slowTot_sci41);

    // ======== Hit Pattern ======== //
    c_gplast_hitpatterns = new TCanvas("c_gplast_hitpatterns","gPlast Hit Pattern",1200,800);
    h1_gplast_hitpatterns = MakeTH1(dir_gplast_hitpattern, "F", "h1_gplast_hitpattern","Detector hit patterns",nDetectors,1,nDetectors+1, "Detector ID", kRed-3, kBlack);
    h1_gplast_hitpatterns->Draw();
    c_gplast_hitpatterns->cd();
    dir_gplast_hitpattern->Append(c_gplast_hitpatterns);

    // tamex card hit pattern
    c_gplast_tamex_card_hitpattern  = new TCanvas("c_gplast_tamex_card_hitpattern","gPlast Tamex Card Hit Pattern",1200,800);
    c_gplast_tamex_card_hitpattern->Divide(5,(nTamexBoards%5==0) ? (nTamexBoards/5) : (nTamexBoards/5 + 1));
    for (int ihist = 0; ihist < nTamexBoards; ihist++){
        c_gplast_tamex_card_hitpattern->cd(ihist+1);
        h1_gplast_tamex_card_hitpattern[ihist] = MakeTH1(dir_gplast_hitpattern, "I", Form("h1_gplast_tamex_card_hitpattern_%d",ihist),Form("gPlast Tamex Card Hit Pattern %d",ihist),16,1,17, "Channel", kRed-3, kBlack); // hard-coded : number of channels per TAMEX
        h1_gplast_tamex_card_hitpattern[ihist]->Draw();
    }
    c_gplast_tamex_card_hitpattern->cd(0);
    dir_gplast_hitpattern->Append(c_gplast_tamex_card_hitpattern);

    // position hit pattern
    c_gplast_position_hitpattern  = new TCanvas("c_gplast_position_hitpattern","gPlast Hit Pattern vs Position",1200,800);
    c_gplast_position_hitpattern->Divide(nPos);
    for (int ihist = 0; ihist < nPos; ihist++){
        c_gplast_position_hitpattern->cd(ihist+1);
        auto str_pos = TgPlastConfiguration::PositionToString(static_cast<TgPlastConfiguration::Position>(ihist+1)); //EXTRA = 0
        h1_gplast_position_hitpattern[ihist] = MakeTH1(dir_gplast_hitpattern, "I", Form("h1_gplast_hitpattern_%s",str_pos.data()),Form("gPlast Hit Pattern %s",str_pos.data()),16,1,17, "Position", kRed-3, kBlack); // hard-coded : number of channels per TAMEX
        h1_gplast_position_hitpattern[ihist]->Draw();
    }
    c_gplast_position_hitpattern->cd(0);
    dir_gplast_hitpattern->Append(c_gplast_position_hitpattern);

    // hit pattern on/off spill
    for (int s=0;s<2;s++){
        const char* tag[2] = {"off","on"};
        c_gplast_position_hitpattern_spill[s]  = new TCanvas(Form("c_gplast_position_hitpattern_spill_%s",tag[s]),Form("gPlast Hit Pattern vs Position Spill %s",tag[s]),1200,800);
        c_gplast_position_hitpattern_spill[s]->Divide(nPos);
        for (int ihist = 0; ihist < nPos; ihist++){
            c_gplast_position_hitpattern_spill[s]->cd(ihist+1);
            auto str_pos = TgPlastConfiguration::PositionToString(static_cast<TgPlastConfiguration::Position>(ihist+1)); //EXTRA = 0
            h1_gplast_position_hitpattern_spill[s][ihist] = MakeTH1(dir_gplast_hitpattern, "I", Form("h1_gplast_hitpattern_%s_spill_%s",str_pos.data(),tag[s]),Form("gPlast Hit Pattern %s (spill %s)",str_pos.data(),tag[s]),16,1,17, "Position", kRed-3, kBlack); // hard-coded : number of channels per TAMEX
            h1_gplast_position_hitpattern_spill[s][ihist]->Draw();
        }
        c_gplast_position_hitpattern_spill[s]->cd(0);
        dir_gplast_hitpattern->Append(c_gplast_position_hitpattern_spill[s]);
    }

    // WR
    c_gplast_wr_time_diff  = new TCanvas("c_gplast_wr_time_diff","gPlast WR time difference",1200,800);
    h1_gplast_wr_time_diff = MakeTH1(dir_gplast_hitpattern, "F", "h1_gplast_wr_time_diff","gPlast WR time difference",1e3,-1e2,5e5, "White Rabbit Event Time Difference [ns]", kViolet, kBlue+2);
    h1_gplast_wr_time_diff->Draw();
    c_gplast_wr_time_diff->cd();
    dir_gplast_hitpattern->Append(c_gplast_wr_time_diff);

    // ======== Coincidences ======== //
    // multiplicity
    c_gplast_multiplicity = new TCanvas("c_gplast_multiplicity","gPlast multiplicity spectrum",1200,800);
    h1_gplast_multiplicity = MakeTH1(dir_gplast_coincidences, "F", "h1_gplast_multiplicity","gPlast multiplicity",nDetectors,1,nDetectors+1, "Channel Multiplicity", kRed-3, kBlack);
    h1_gplast_multiplicity->Draw();
    c_gplast_multiplicity->SetLogy();
    c_gplast_multiplicity->cd();
    dir_gplast_coincidences->Append(c_gplast_multiplicity);

    // coincidences
    c_gplast_coincidences_matrix = new TCanvas("c_gplast_coincidences_matrix","gPlast matrix of coincidences",1200,800);
    h2_gplast_coincidence_matrix = MakeTH2(dir_gplast_coincidences, "F", "h2_gplast_coincidence_matrix","gPlast coincidences",nSignalsPlastic,1,nSignalsPlastic+1,nSignalsPlastic,1,nSignalsPlastic+1, "Position (T-B-L-R)", "Position (T-B-L-R)");
    h2_gplast_coincidence_matrix->Draw();
    c_gplast_coincidences_matrix->cd();
    dir_gplast_coincidences->Append(c_gplast_coincidences_matrix);

    // dt distribution (all hits)
    c_gplast_time_differences = new TCanvas("c_gplast_time_differences","gPlast time differences (coincidences)",1200,800);
    h1_gplast_coinc_time_diff = MakeTH1(dir_gplast_coincidences, "F", "h1_gplast_coinc_time_diff","gPlast time difference between hits (mult >= 2)",1e3,-500,500, "Hit Time Difference [ns]", kMagenta, kBlue+2);
    h1_gplast_coinc_time_diff->Draw();
    c_gplast_time_differences->cd();
    dir_gplast_coincidences->Append(c_gplast_time_differences);

    // dt distribution vs position
    c_gplast_time_differences_position = new TCanvas("c_gplast_time_differences_position","gPlast time differences (coincidences) vs Position",1200,800);
    c_gplast_time_differences_position->Divide(nPos);
    for (int ihist = 0; ihist < nPos; ihist++){
        c_gplast_time_differences_position->cd(ihist+1);
        auto str_pos = TgPlastConfiguration::PositionToString(static_cast<TgPlastConfiguration::Position>(ihist+1)); //EXTRA = 0
        h1_gplast_coinc_time_diff_position[ihist] = MakeTH1(dir_gplast_coincidences, "F", Form("h1_gplast_coinc_time_diff_%s",str_pos.data()),Form("gPlast time difference between hits (mult >= 2) - %s",str_pos.data()),1e3,-500,500, "Hit Time Difference [ns]", kMagenta, kBlue+2);
        h1_gplast_coinc_time_diff_position[ihist]->Draw();
    }
    c_gplast_time_differences_position->cd();
    dir_gplast_coincidences->Append(c_gplast_time_differences_position);

    // dt distribution left/right (x) and top/bottom (y)
    c_gplast_time_differences_xy = new TCanvas("c_gplast_time_differences_xy","gPlast time differences (coincidences) vs X/Y",1200,800);
    c_gplast_time_differences_xy->Divide(2);
    c_gplast_time_differences_xy->cd(1);
    h1_gplast_coinc_time_diff_vs_x = MakeTH1(dir_gplast_coincidences, "F", "h1_gplast_coinc_time_diff_vs_x","gPlast time difference between hits (mult >= 2) - Left vs Right",1e3,-500,500, "Hit Time Difference Left-Right [ns]", kMagenta, kBlue+2);
    h1_gplast_coinc_time_diff_vs_x->Draw();
    c_gplast_time_differences_xy->cd(2);
    h1_gplast_coinc_time_diff_vs_y = MakeTH1(dir_gplast_coincidences, "F", "h1_gplast_coinc_time_diff_vs_y","gPlast time difference between hits (mult >= 2) - Top vs Bottom",1e3,-500,500, "Hit Time Difference Top-Bottom [ns]", kMagenta, kBlue+2);
    h1_gplast_coinc_time_diff_vs_y->Draw();
    c_gplast_time_differences_xy->cd();
    dir_gplast_coincidences->Append(c_gplast_time_differences_xy);

    // ======== Time spectra ======== //
    c_gplast_time_spectra_sci41  = new TCanvas("c_bplast_time_spectra","Fast bPlast time spectra vs sci41",1200,800);
    c_gplast_time_spectra_sci41->Divide(5,(nDetectors%5==0) ? (nDetectors/5) : (nDetectors/5 + 1));
    for (int ihist = 1; ihist <= nDetectors; ihist++)
    {
        c_gplast_time_spectra_sci41->cd(ihist);
        h1_gplast_time_spectra_sci41[ihist] = MakeTH1(dir_gplast_time_spectra, "F", Form("h1_gplast_time_spectra_vs_sci41_%d",ihist),Form("gPlast time difference (sci41) %d",ihist),10000,-1000,1000, "Time difference [ns]", kMagenta, kBlue+2);
        // hard-coded -> binning should be given in macro !
        h1_gplast_time_spectra_sci41[ihist]->Draw();
    }
    c_gplast_time_spectra_sci41->cd(0);
    dir_gplast_time_spectra->Append(c_gplast_time_spectra_sci41);

    c_gplast_time_vs_position = new TCanvas("c_gplast_time_vs_position","gPlast Time Differences vs Position",1200,800);
    c_gplast_time_vs_position->Divide(nPos);
    for (int ihist = 0; ihist < nPos; ihist++){
        c_gplast_time_vs_position->cd(ihist+1);
        auto str_pos = TgPlastConfiguration::PositionToString(static_cast<TgPlastConfiguration::Position>(ihist+1)); //EXTRA = 0
        h2_gplast_time_vs_position[ihist] = MakeTH2(dir_gplast_time_spectra, "F", Form("h2_gplast_time_vs_position_%s",str_pos.data()),Form("gPlast Time vs Position %s",str_pos.data()),6000,-100,500,16,1,17,"Time difference with first hit [ns]","Position"); // hard-coded
        h2_gplast_time_vs_position[ihist]->Draw();
    }
    c_gplast_time_vs_position->cd(0);
    dir_gplast_time_spectra->Append(c_gplast_time_vs_position);

    c_gplast_time_vs_position_sci41 = new TCanvas("c_gplast_time_vs_position_sci41","gPlast Time Differences (sci41) vs Position",1200,800);
    c_gplast_time_vs_position_sci41->Divide(nPos);
    for (int ihist = 0; ihist < nPos; ihist++){
        c_gplast_time_vs_position_sci41->cd(ihist+1);
        auto str_pos = TgPlastConfiguration::PositionToString(static_cast<TgPlastConfiguration::Position>(ihist+1)); //EXTRA = 0
        h2_gplast_time_sci41_vs_position[ihist] = MakeTH2(dir_gplast_time_spectra, "F", Form("h2_gplast_time_sci41_vs_position_%s",str_pos.data()),Form("gPlast Time (sci41) vs Position %s",str_pos.data()),10000,-100,900,16,1,17,"Time difference with sci41 [ns]","Position"); // hard-coded
        h2_gplast_time_sci41_vs_position[ihist]->Draw();
    }
    c_gplast_time_vs_position_sci41->cd(0);
    dir_gplast_time_spectra->Append(c_gplast_time_vs_position_sci41);

    // ======== Time of Flight vs SCI4x ======== //
    c_gplast_tof_sci41 = new TCanvas("c_gplast_tof_sci41","gPlast ToF vs sci41",1200,800);
    c_gplast_tof_sci41->Divide(2);
    c_gplast_tof_sci41->cd(1);
    h1_gplast_tof_sci41 = MakeTH1(dir_gplast_ToF_sci, "F", "h1_gplast_tof_sci41","gPlast ToF vs sci41 - no veto",1000,0,500, "Time of Flight with respect to sci41 [ns]", kMagenta+1, kBlue+2); //hard-coded
    h1_gplast_tof_sci41->GetXaxis()->SetTitle(" (ns)");
    h1_gplast_tof_sci41->Draw();
    c_gplast_tof_sci41->cd(2);
    h1_gplast_tof_sci41_veto = MakeTH1(dir_gplast_ToF_sci, "F", "h1_gplast_tof_sci41_veto","gPlast ToF vs sci41 - vetoed",1000,0,500, "Time of Flight with respect to sci41 [ns]", kMagenta+1, kBlue+2); //hard-coded
    h1_gplast_tof_sci41_veto->Draw();
    c_gplast_tof_sci41->cd();
    dir_gplast_ToF_sci->Append(c_gplast_tof_sci41);

    c_gplast_tof_sci42 = new TCanvas("c_gplast_tof_sci42","gPlast ToF vs sci42",1200,800);
    c_gplast_tof_sci42->Divide(2);
    c_gplast_tof_sci42->cd(1);
    h1_gplast_tof_sci42 = MakeTH1(dir_gplast_ToF_sci, "F", "h1_gplast_tof_sci42","gPlast ToF vs sci42 - no veto",1000,0,500, "Time of Flight with respect to sci42 [ns]", kMagenta+1, kBlue+2); //hard-coded
    h1_gplast_tof_sci42->Draw();
    c_gplast_tof_sci42->cd(2);
    h1_gplast_tof_sci42_veto = MakeTH1(dir_gplast_ToF_sci, "F", "h1_gplast_tof_sci42_veto","gPlast ToF vs sci42 - vetoed",1000,0,500, "Time of Flight with respect to sci42 [ns]", kMagenta+1, kBlue+2); //hard-coded
    h1_gplast_tof_sci42_veto->Draw();
    c_gplast_tof_sci42->cd();
    dir_gplast_ToF_sci->Append(c_gplast_tof_sci42);

    c_gplast_tof_sci43 = new TCanvas("c_gplast_tof_sci43","gPlast ToF vs sci43",1200,800);
    h1_gplast_tof_sci43 = MakeTH1(dir_gplast_ToF_sci, "F", "h1_gplast_tof_sci43","gPlast ToF vs sci43",1000,-500,0, "Time of Flight with respect to sci43 [ns]", kMagenta+1, kBlue+2); //hard-coded
    h1_gplast_tof_sci43->Draw();
    c_gplast_tof_sci43->cd();
    dir_gplast_ToF_sci->Append(c_gplast_tof_sci43);

    dir_gplast->cd();

    run->GetHttpServer()->RegisterCommand("Reset_gPlast_Histo", Form("/Objects/%s/->Reset_Histo()", GetName()));
    run->GetHttpServer()->RegisterCommand("Snapshot_gPlast_Histo", Form("/Objects/%s/->Snapshot_Histo()", GetName()));

    return kSUCCESS;
}

void gPlastOnlineSpectra::Reset_Histo()
{
    c4LOG(info, "Resetting gPlast histograms.");
    std::vector<TDirectory*> dirs = {dir_gplast_slowToT, dir_gplast_fastToT, dir_gplast_fast_vs_slow, dir_gplast_hitpattern, dir_gplast_slowToT_total, dir_gplast_coincidences, dir_gplast_time_spectra, dir_gplast_ToF_sci};
    for (const auto* dir:dirs){
        for(TObject* obj: *dir->GetList()){
            if (obj->InheritsFrom("TH1")){
                TH1* hist = dynamic_cast<TH1*>(obj);
                if (hist)
                    hist->Reset();
            }
        }
    }
    c4LOG(info, "gPlast histograms reset.");
}

// make a date and time stamped folder with pngs of the histograms and .root file and save them
void gPlastOnlineSpectra::Snapshot_Histo()
{
    c4LOG(info, "Snapshotting gPlast histograms.");
    c4LOG(warn,"THIS IS NOT UP TO DATE !");
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

    bool has_sci41_signal = false;
    bool has_sci42_signal = false;
    bool has_sci43_signal = false;
    double time_sci41_signal = 0.;
    double time_sci42_signal = 0.;
    double time_sci43_signal = 0.;
    double total_slowTot = 0.;
    double time_first_hit = -1;
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

            if (detector_id == gplast_conf->SC41L() || detector_id == gplast_conf->SC41R()){
                if (!has_sci41_signal || time_sci41_signal > fast_lead_time)
                    time_sci41_signal = fast_lead_time;
                has_sci41_signal = true;
            }
            else if (detector_id == gplast_conf->SC42L() || detector_id == gplast_conf->SC42R()){
                if (!has_sci42_signal || time_sci42_signal > fast_lead_time)
                    time_sci42_signal = fast_lead_time;
                has_sci42_signal = true;
            }
            else if (detector_id == gplast_conf->SC43L() || detector_id == gplast_conf->SC43R()){
                if (!has_sci43_signal || time_sci43_signal > fast_lead_time)
                    time_sci43_signal = fast_lead_time;
                has_sci43_signal = true;
            }

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
                total_slowTot += slowToT;
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

            // get the timestamp of the first hit from this event
            if (time_first_hit < 0) time_first_hit = fast_lead_time;
            else{
                if (time_first_hit > fast_lead_time){
                    time_first_hit = fast_lead_time;
                }
            }
        }

        if(event_multiplicity >= 1) {
            h1_gplast_multiplicity->Fill(event_multiplicity);

            if (has_sci43_signal){
                // veto
                h1_gplast_slowToT_total_veto->Fill(total_slowTot);
                if (has_sci41_signal)
                    h1_gplast_tof_sci41_veto->Fill(time_first_hit-time_sci41_signal);
                if (has_sci42_signal)
                    h1_gplast_tof_sci42_veto->Fill(time_first_hit-time_sci42_signal);
                h1_gplast_tof_sci43->Fill(time_first_hit-time_sci43_signal);
            }
            else{
                h1_gplast_slowToT_total->Fill(total_slowTot);
                if (has_sci41_signal){
                    h1_gplast_tof_sci41->Fill(time_first_hit-time_sci41_signal);
                    h2_gplast_slowToT_total_sci41->Fill(total_slowTot,time_first_hit-time_sci41_signal);
                }
                if (has_sci42_signal)
                    h1_gplast_tof_sci42->Fill(time_first_hit-time_sci42_signal);
            }

            for (Int_t ihit = 0; ihit < nHits; ihit++)
            {
                gPlastTwinpeaksCalData* hit = (gPlastTwinpeaksCalData*)fHitgPlastTwinpeaks->At(ihit);
                if (!hit) continue;

                int detector_id = hit->Get_detector_id();
                if (detector_id > nDetectors || detector_id < 0) continue;

                auto pos = gplast_pos.at(detector_id);
                if (pos.first != TgPlastConfiguration::EXTRA){
                    double t_hit = hit->Get_fast_lead_time();
                    if (has_sci41_signal){
                        // spill on
                        // - dt vs sci (1d/2d)
                        // - hit pattern
                        // - tof between gplast and sci41/42/43
                        double delta_t = t_hit - time_sci41_signal;
                        h1_gplast_time_spectra_sci41[detector_id]->Fill(delta_t);
                        h1_gplast_position_hitpattern_spill[1][pos.first-1]->Fill(pos.second);
                        h2_gplast_time_sci41_vs_position[pos.first-1]->Fill(delta_t,pos.second);
                    }
                    else{
                        // spill off
                        // - hit pattern
                        h1_gplast_position_hitpattern_spill[0][pos.first-1]->Fill(pos.second);
                    }

                    if (event_multiplicity > 1){
                        // coincidences (between gPlast channels)
                        //? should this be moved to nearline ?
                        for (Int_t ihit2 = ihit+1; ihit2 < nHits; ihit2++)
                        {
                            gPlastTwinpeaksCalData* hit2 = (gPlastTwinpeaksCalData*)fHitgPlastTwinpeaks->At(ihit2);
                            if (!hit2) continue;

                            int detector2_id = hit2->Get_detector_id();
                            auto pos2 = gplast_pos.at(detector2_id);
                            if (pos2.first != TgPlastConfiguration::EXTRA){
                                // fill (i,j) and (j,i)
                                h2_gplast_coincidence_matrix->Fill((pos.first-1)*16+pos.second,(pos2.first-1)*16+pos2.second);
                                h2_gplast_coincidence_matrix->Fill((pos2.first-1)*16+pos2.second,(pos.first-1)*16+pos.second);
                            }
                            h1_gplast_coinc_time_diff->Fill(hit->Get_fast_lead_time()-hit2->Get_fast_lead_time());
                            if (pos.first == pos2.first){
                                h1_gplast_coinc_time_diff_position[pos.first-1]->Fill(hit->Get_fast_lead_time()-hit2->Get_fast_lead_time());
                            }
                            else{
                                if (pos.first == TgPlastConfiguration::TOP && pos2.first == TgPlastConfiguration::BOTTOM){
                                    h1_gplast_coinc_time_diff_vs_y->Fill(hit->Get_fast_lead_time()-hit2->Get_fast_lead_time());
                                }
                                else if (pos2.first == TgPlastConfiguration::TOP && pos.first == TgPlastConfiguration::BOTTOM){
                                    h1_gplast_coinc_time_diff_vs_y->Fill(hit2->Get_fast_lead_time()-hit->Get_fast_lead_time());
                                }
                                else if (pos.first == TgPlastConfiguration::LEFT && pos2.first == TgPlastConfiguration::RIGHT){
                                    h1_gplast_coinc_time_diff_vs_x->Fill(hit->Get_fast_lead_time()-hit2->Get_fast_lead_time());
                                }
                                else if (pos2.first == TgPlastConfiguration::LEFT && pos.first == TgPlastConfiguration::RIGHT){
                                    h1_gplast_coinc_time_diff_vs_x->Fill(hit2->Get_fast_lead_time()-hit->Get_fast_lead_time());
                                }
                            }
                        }
                        h2_gplast_time_vs_position[pos.first-1]->Fill(t_hit-time_first_hit,pos.second); // is this good ?
                    }

                }
            }
        }
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
