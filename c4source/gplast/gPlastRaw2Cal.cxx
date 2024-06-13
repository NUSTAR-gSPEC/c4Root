// FairRoot
#include "FairTask.h"
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"

// c4
#include "gPlastReader.h"
#include "gPlastTwinpeaksData.h"
#include "gPlastTwinpeaksCalData.h"
#include "TimeMachineData.h"

#include "c4Logger.h"

#include "TClonesArray.h"

#include "gPlastRaw2Cal.h"
#include <string>
#include <chrono>

/*
empty constructor required for FairRoot.
*/
gPlastRaw2Cal::gPlastRaw2Cal()
: FairTask(), 
fNEvents(0),
header(nullptr),
fOnline(kFALSE),
funcal_data(new TClonesArray("gPlastTwinpeaksData")),
fcal_data(new TClonesArray("gPlastTwinpeaksCalData")),
ftime_machine_array(new TClonesArray("TimeMachineData"))
{
    gplast_config = TgPlastConfiguration::GetInstance();
}

/*
Named constructor with verbosity level.
*/
gPlastRaw2Cal::gPlastRaw2Cal(const TString& name, Int_t verbose) 
    : FairTask(name, verbose),
    fNEvents(0),
    header(nullptr),
    fOnline(kFALSE),
    funcal_data(new TClonesArray("gPlastTwinpeaksData")),
    fcal_data(new TClonesArray("gPlastTwinpeaksCalData")),
    ftime_machine_array(new TClonesArray("TimeMachineData"))
{
    gplast_config = TgPlastConfiguration::GetInstance();
}

/*  
Clearing old constructed objects.
*/
gPlastRaw2Cal::~gPlastRaw2Cal(){
    c4LOG(info, "Deleting gPlastRaw2Cal task");
    if (funcal_data) delete funcal_data;
    if (fcal_data) delete fcal_data;
    if (ftime_machine_array) delete ftime_machine_array;
}

/*
This is called AFTER the detector mapping. This picks out the two timemachine channels and writes them to the TimeMachine structure.
*/
void gPlastRaw2Cal::SetTimeMachineChannels(int ftime_machine_undelayed_detector_id, int ftime_machine_delayed_detector_id)
{
    time_machine_delayed_detector_id = ftime_machine_delayed_detector_id;
    time_machine_undelayed_detector_id = ftime_machine_undelayed_detector_id;
}



void gPlastRaw2Cal::SetParContainers()
{
    FairRuntimeDb *rtdb = FairRuntimeDb::instance();
    c4LOG_IF(fatal, NULL == rtdb, "FairRuntimeDb not found.");
}

/*
Initialiser called by the FairRoot manager. Gets the required FairRootManager objects to read and register the data to be written to the tree.
*/
InitStatus gPlastRaw2Cal::Init()
{
    FairRootManager* mgr = FairRootManager::Instance();
    c4LOG_IF(fatal, NULL == mgr, "FairRootManager not found");

    header = (EventHeader*)mgr->GetObject("EventHeader.");
    c4LOG_IF(error, !header, "Branch EventHeader. not found");

    funcal_data = (TClonesArray*)mgr->GetObject("gPlastTwinpeaksData");
    c4LOG_IF(fatal, !funcal_data, "gPlast branch of gPlastTwinpeaksData not found.");
    
    //needs to have the name of the detector subsystem here:
    FairRootManager::Instance()->Register("gPlastTwinpeaksCalData", "gPlast Cal Data", fcal_data, !fOnline);
    FairRootManager::Instance()->Register("gPlastTimeMachineData", "gPlast Time Machine Data", ftime_machine_array, !fOnline);

    fcal_data->Clear();
    funcal_data->Clear();

    DetectorMap_loaded = gplast_config->MappingLoaded();
    if (!DetectorMap_loaded){
        c4LOG(error, "gPlast mapping not loaded !");
    }
    else{
        detector_mapping = gplast_config->Mapping();
        int ntam_map = gplast_config->NTamexBoards();
        int ntam_data = gPlastReader::GetNumberOfTAMEXBoards();
        c4LOG_IF(error, ntam_map != ntam_data, Form("Unconsistent number of TAMEX cards in mapping file (%d) and in the data (%d)",ntam_map,ntam_data));
    }

    return kSUCCESS;
}

/*
Reads a file containing the detector calibrations. To be called before Init. Assumed structure of the file is:
    - aribtrary lines of comments starting with #
    - each entry is a line with four number: (fatima detector id) (a1/slope) (a0/offset)

Raises a fatal error if the detector numbers are not unique.
*/
Bool_t gPlastRaw2Cal::SetDetectorCalFile(TString filename){
    c4LOG(info, "Reading Calibration coefficients.");

    std::ifstream cal_map_file (filename);

    int rdetector_id; // temp read variables
    double a0,a1;
    
    // loading and reading detector calibration file. Assumes the first line in the file is num-modules used
    while(!cal_map_file.eof()){
        if(cal_map_file.peek()=='#') cal_map_file.ignore(256,'\n');
        else{
            cal_map_file >> rdetector_id >> a1 >> a0;
            std::pair<double,double> cals = {a0,a1};
            calibration_coeffs.insert(std::pair<int,std::pair<double,double>>{rdetector_id,cals});
            cal_map_file.ignore(256,'\n');
            
            auto it = calibration_coeffs.find(rdetector_id);
            if (it != calibration_coeffs.end()) c4LOG(fatal,Form("Calibration coefficients not unique. Multiple entries of (detector id = %i)",rdetector_id));
        }
    }
    DetectorCal_loaded = 1;
    cal_map_file.close();  
    return 0; 
};

/*
Writes the detector map to console.
*/
void gPlastRaw2Cal::PrintDetectorMap()
{
    if (DetectorMap_loaded){
        auto position_map = gplast_config->Positioning();
        for (const auto &entry : detector_mapping){
            std::cout << "TAM_id: " << entry.first.first
                      << " TAM_ch: " << entry.first.second;
            auto pos = position_map.at(entry.second);
            std::cout << " DETECTOR_id: " << entry.second
                      << " POSITION: " << TgPlastConfiguration::PositionToString(pos.first)
                      << "_" << pos.second << "\n";
        }
        std::cout << std::endl;
    }
    else{
        c4LOG(warn, "Detector map is not load. Cannot print.");
    }
}

/*
Write the detector calibration to console.
*/
void gPlastRaw2Cal::PrintDetectorCal(){
    if (DetectorCal_loaded){
        for (const auto& entry : calibration_coeffs){
            std::cout << "DETECTORID: " << entry.first;
            std::cout << " a0: " << entry.second.first << " a1: " << entry.second.second << "\n";
        }
    }
    else{
        c4LOG(warn, "Cal map is not load. Cannot print.");
    }
}        

/*
The event loop executable. This is where the events are analyzed. Only used implicitly by FairRoot during Run().

Further the hits are matched slow + fast and assigned from the internal Twinpeaks channnel number to the detector number if DetectorMap is loaded.
Assumes that fast hits always preceedes slow hits. 

Writes the times in ns!
*/
void gPlastRaw2Cal::Exec(Option_t* option)
{
    
    auto start = std::chrono::high_resolution_clock::now();

    if (funcal_data && funcal_data->GetEntriesFast() > 1){ // only get events with two hits.or more
        Int_t event_multiplicity = funcal_data->GetEntriesFast();
        // event mulitiplicity loop
        for (Int_t ihit = 0; ihit < event_multiplicity; ihit++){

            gPlastTwinpeaksData* first_hit_in_fast_channel = (gPlastTwinpeaksData*)funcal_data->At(ihit);

            // under the assumption fast-slow always follows:
            //assume that only matched lead-trail hits are written.
            if (first_hit_in_fast_channel->Get_ch_ID()%2==0) {continue;} //get the first odd numbered channel

            int hits_in_fast_channel = 1;
            int hits_in_slow_channel = 0;

            int look_ahead_counter = 1;
            bool all_hits_in_fast_slow_found = false;
            while (!all_hits_in_fast_slow_found)
            {
                if (ihit+look_ahead_counter >= event_multiplicity) break;
                gPlastTwinpeaksData* this_hit = (gPlastTwinpeaksData*)funcal_data->At(ihit+look_ahead_counter);

                //c4LOG(info,this_hit->Get_ch_ID());
                if (this_hit->Get_ch_ID() == first_hit_in_fast_channel->Get_ch_ID() && this_hit->Get_board_id() == first_hit_in_fast_channel->Get_board_id()) hits_in_fast_channel++;
                else if (this_hit->Get_ch_ID() == first_hit_in_fast_channel->Get_ch_ID()+1 && this_hit->Get_board_id() == first_hit_in_fast_channel->Get_board_id()) hits_in_slow_channel++;
                else all_hits_in_fast_slow_found = true;
                look_ahead_counter++;
            }


            //c4LOG(info,Form("fast hits = %i, slow hits = %i, look ahead counter = %i",hits_in_fast_channel,hits_in_slow_channel,look_ahead_counter));
            if (hits_in_fast_channel != hits_in_slow_channel) {
                //break condition - cant recover.
                ihit = hits_in_fast_channel + hits_in_slow_channel - 1 + ihit;
                continue;
            }


            for (int hitnr = 0; hitnr<hits_in_fast_channel; hitnr++)
            {

                funcal_hit = (gPlastTwinpeaksData*)funcal_data->At(ihit+hitnr);
                funcal_hit_next = (gPlastTwinpeaksData*)funcal_data->At(ihit+hitnr+hits_in_fast_channel);
            
                if (funcal_hit_next->Get_ch_ID() != funcal_hit->Get_ch_ID()+1)
                { // this assumption seems empirically true - no events are filled when reverse order is put.
                    fNunmatched++; continue;
                }

                if (funcal_hit_next->Get_board_id() != funcal_hit->Get_board_id()){
                    continue;
                }


                //from here the funcalhitpartner is the slow branch and funcal_hit the fast:
                detector_id = -1;
                if (DetectorMap_loaded)
                {
                    std::pair<int, int> unmapped_det { funcal_hit->Get_board_id(), (funcal_hit->Get_ch_ID()+1)/2 };

                    if (auto result_find = detector_mapping.find(unmapped_det); result_find != detector_mapping.end())
                    {
                        detector_id = result_find->second;
                        if (detector_id == -1) { fNunmatched++; continue; } // if only one event is left
                    }
                    // else c4LOG(warn, "Detector mapping is not complete!"); // comment it for now to avoid flooding the output
                }

                if (funcal_hit_next->Get_trail_epoch_counter() == 0) continue; // missing trail in either

                // I am slightly worried about round-off errors by this method (but as far as i can see the maximum epoch counter values is not so large that the digits are suppressed but it is something to keep in mind). However constructing the times like this makes it very easy to use.
                fast_lead_time = static_cast<double>(funcal_hit->Get_lead_epoch_counter()) * 10.24e3
                            + static_cast<double>(funcal_hit->Get_lead_coarse_T()) * 5.0
                            - static_cast<double>(funcal_hit->Get_lead_fine_T());

                fast_trail_time = static_cast<double>(funcal_hit->Get_trail_epoch_counter()) * 10.24e3
                                + static_cast<double>(funcal_hit->Get_trail_coarse_T()) * 5.0
                                - static_cast<double>(funcal_hit->Get_trail_fine_T());

                slow_lead_time = static_cast<double>(funcal_hit_next->Get_lead_epoch_counter()) * 10.24e3
                                + static_cast<double>(funcal_hit_next->Get_lead_coarse_T()) * 5.0
                                - static_cast<double>(funcal_hit_next->Get_lead_fine_T());

                slow_trail_time = static_cast<double>(funcal_hit_next->Get_trail_epoch_counter()) * 10.24e3
                                + static_cast<double>(funcal_hit_next->Get_trail_coarse_T()) * 5.0
                                - static_cast<double>(funcal_hit_next->Get_trail_fine_T());

                fast_ToT =  fast_trail_time - fast_lead_time;
                slow_ToT =  slow_trail_time - slow_lead_time;

            //if (detector_id == 0 || detector_id == 1) c4LOG(info,Form("id = %i, fast lead = %f, fast trail = %f, fast ToT = %f",detector_id,fast_lead_time,fast_trail_time,fast_ToT));

            // TODO: add channel calibration if needed.
            /*if (gplast_config->MappingLoaded())
            {
                if (gplast_config->CalibrationCoefficientsLoaded())
                {
                    // etc etc
                }
            }*/
            
            if (((detector_id == gplast_config->TM_Delayed()) || (detector_id == gplast_config->TM_Undelayed())) && gplast_config->TM_Delayed() != 0 && gplast_config->TM_Undelayed() != 0)
            { // currently only gets the TM if it also matches it slow-fast...
                new ((*ftime_machine_array)[ftime_machine_array->GetEntriesFast()]) TimeMachineData((detector_id == gplast_config->TM_Undelayed()) ? (fast_lead_time) : (0), (detector_id == gplast_config->TM_Undelayed()) ? (0) : (fast_lead_time), funcal_hit->Get_wr_subsystem_id(), funcal_hit->Get_wr_t() );
            }

            // else ?
            new ((*fcal_data)[fcal_data->GetEntriesFast()]) gPlastTwinpeaksCalData(
                funcal_hit->Get_board_id(),
                (int)((funcal_hit->Get_ch_ID()+1)/2),
                detector_id,
                slow_lead_time,
                slow_trail_time,
                fast_lead_time,
                fast_trail_time,
                fast_ToT,
                slow_ToT,
                funcal_hit->Get_wr_subsystem_id(),
                funcal_hit->Get_wr_t());
            
            
            fNEvents++;
            //ihit++; //increment it by one extra.
            }
        }
        fExecs++; // count once every time we do something with gplast
    } // if gplast hit exists

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    total_time_microsecs += duration.count();


}


/*
THIS FUNCTION IS EXTREMELY IMPORTANT!!!!

Clears the TClonesArray used in the function. If they are not cleared after each event they will eat all your RAM.
*/
void gPlastRaw2Cal::FinishEvent()
{
    // reset output array
    funcal_data->Clear();
    fcal_data->Clear();
    ftime_machine_array->Clear();
};

/*
Some stats are written when finishing.
*/
void gPlastRaw2Cal::FinishTask()
{
    c4LOG(info, Form("Wrote %i events.",fNEvents));
    c4LOG(info, Form("%i events are unmatched (not written).",fNunmatched));
    c4LOG(info, "Average execution time: " << (double)total_time_microsecs/fExecs << " microseconds.");
}


ClassImp(gPlastRaw2Cal)
