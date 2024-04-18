#include "FairRootManager.h"
#include "FairRunOnline.h"
#include "FairTask.h"

#include "FrsAidaCorrelationsOnline.h"
#include "FrsHitData.h"
#include "AidaHitData.h"
#include "TAidaConfiguration.h"
#include "c4Logger.h"

#include "TDirectory.h"
#include "THttpServer.h"

FrsAidaCorrelationsOnline::FrsAidaCorrelationsOnline(std::vector<FrsGate*> fg)
    :   FairTask()
    ,   fNEvents(0)
    ,   header(nullptr)
    ,   fFrsHitArray(new TClonesArray("FrsHitData"))
    ,   fAidaImplants(new std::vector<AidaHit>)
{
    FrsGates = fg;
    aida_config = TAidaConfiguration::GetInstance();
    frs_config = TFrsConfiguration::GetInstance();
    correl_config = TCorrelationsConfiguration::GetInstance();
    Correl = correl_config->CorrelationsMap();
}

FrsAidaCorrelationsOnline::FrsAidaCorrelationsOnline(const TString& name, Int_t verbose)
    :   FairTask(name, verbose)
    ,   fNEvents(0)
    ,   header(nullptr)
    ,   fFrsHitArray(new TClonesArray("FrsHitData"))
    ,   fAidaImplants(new std::vector<AidaHit>)
{
    aida_config = TAidaConfiguration::GetInstance();
    frs_config = TFrsConfiguration::GetInstance();
    correl_config = TCorrelationsConfiguration::GetInstance();
    Correl = correl_config->CorrelationsMap();
}

FrsAidaCorrelationsOnline::~FrsAidaCorrelationsOnline()
{
    c4LOG(info, "Deleting FrsAidaCorrelationsOnline task.");
    if (fFrsHitArray) delete fFrsHitArray;
}

InitStatus FrsAidaCorrelationsOnline::Init()
{
    FairRootManager* mgr = FairRootManager::Instance();
    c4LOG_IF(fatal, NULL == mgr, "FairRootManager not found");

    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);

    header = (EventHeader*)mgr->GetObject("EventHeader.");
    c4LOG_IF(error, !header, "EventHeader. not found!");

    fFrsHitArray = (TClonesArray*)mgr->GetObject("FrsHitData");
    c4LOG_IF(fatal, !fFrsHitArray, "FrsHitData branch not found!");

    fAidaImplants = mgr->InitObjectAs<decltype(fAidaImplants)>("AidaImplantHits");
    c4LOG_IF(fatal, !fAidaImplants, "Branch AidaImplantHits not found!");

    histograms = (TFolder*)mgr->GetObject("Histograms");
    
    TDirectory::TContext ctx(nullptr);

    // look for FRS directory, create it if not found
    dir_frs = (TDirectory*)mgr->GetObject("FRS");
    if (dir_frs == nullptr) 
    {
        dir_frs = new TDirectory("FRS Online", "FRS Online", "", 0);
        mgr->Register("FRS", "FRS Online Directory", dir_frs, false); // allow other tasks to find this
        histograms->Add(dir_frs);
    }

    dir_frs_aida_corr = dir_frs->mkdir("FRS-AIDA Correlations");

    // init tgraphs/histograms here
    dir_frs_aida_corr->cd();
    int bins = FrsGates.size() * aida_config->DSSDs();
    h1_stopped_implanted_passed_gate = new TH1I("h1_stopped_implanted_passed_gate", "Stopped ions by DSSD, passing FRS PID Gates", bins, 0, bins);

    for (int i = 0; i < bins; i++)
    {
        const char* labelDSSD;
        const char* labelPID;
        int DSSD = (i % aida_config->DSSDs() + 1);

        if (i % aida_config->DSSDs() == 0)
        {
            std::string stringPID = "PID: " + FrsGates[i/aida_config->DSSDs()]->GetName();
            labelPID = stringPID.c_str();
            h1_stopped_implanted_passed_gate->GetXaxis()->SetBinLabel(i+1, Form("#splitline{DSSD %i}{%s}",DSSD,labelPID));
        }
        else h1_stopped_implanted_passed_gate->GetXaxis()->SetBinLabel(i+1, Form("DSSD %i", DSSD));

    }

    stopped_implant_passed_gate_count = new int*[FrsGates.size()];
    for (int i = 0; i < FrsGates.size(); i++) stopped_implant_passed_gate_count[i] = new int[aida_config->DSSDs()];

    return kSUCCESS;

}


void FrsAidaCorrelationsOnline::Exec(Option_t* option)
{ 
    int multFRS = 0;
    multFRS = fFrsHitArray->GetEntriesFast();
    if (multFRS == 0) return;

    if (fFrsHitArray && fFrsHitArray->GetEntriesFast() > 0)
    {
        Int_t nHits = fFrsHitArray->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            FrsHitData* FrsHit = (FrsHitData*)fFrsHitArray->At(ihit);
            if (!FrsHit) continue;

            for (auto & i : *fAidaImplants)
            {
                AidaHit hit = i;
                
                if (hit.Time - FrsHit->Get_wr_t() > Correl["FRS-AIDA WR Gate"][0] && hit.Time - FrsHit->Get_wr_t() < Correl["FRS-AIDA WR Gate"][1])
                {
                    if (hit.Stopped)
                    {
                        if (!FrsGates.empty())
                        {
                            for (int gate = 0; gate < FrsGates.size(); gate++)
                            {
                                if (FrsGates[gate]->PassedGate(FrsHit->Get_ID_z(), FrsHit->Get_ID_z2(), FrsHit->Get_ID_x2(), FrsHit->Get_ID_x4(), FrsHit->Get_ID_AoQ(), FrsHit->Get_ID_dEdeg()))
                                {
                                    h1_stopped_implanted_passed_gate->Fill(gate * aida_config->DSSDs() + hit.DSSD - 1);

                                    stopped_implant_passed_gate_count[gate][hit.DSSD - 1]++;
                                }
                            }


                            // calculate ratio
                        }
                        
                    } // Stopped

                } // FRS-AIDA WR Gate
            } // Aida Implants
        
        } // nHits

    } // GetEntries() > 0

    fNEvents++;
}

void FrsAidaCorrelationsOnline::FinishEvent()
{

}

void FrsAidaCorrelationsOnline::FinishTask()
{

}


ClassImp(FrsAidaCorrelationsOnline)