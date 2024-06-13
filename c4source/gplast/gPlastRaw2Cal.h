#ifndef gPlastRaw2Cal_H
#define gPlastRaw2Cal_H

#include "FairTask.h"
#include "TgPlastConfiguration.h"


class TClonesArray;
class EventHeader;
class gPlastTwinpeaksData;
class gPlastTwinpeaksCalData;
class TimeMachineData;

class gPlastRaw2Cal : public FairTask
{
    public:
        gPlastRaw2Cal();

        gPlastRaw2Cal(const TString& name, Int_t verbose);

        ~gPlastRaw2Cal();

        void PrintDetectorMap();
        void PrintDetectorCal();
        
        Bool_t SetDetectorCalFile(TString);

        void SetTimeMachineChannels(int ftime_machine_delayed_detector_id, int ftime_machine_undelayed_detector_id);

        void Exec(Option_t* option);

        void FinishEvent();
        void FinishTask();

        void SetOnline(Bool_t set_online){fOnline = set_online;}

        virtual void SetParContainers();

        virtual InitStatus Init();

    private:

        TgPlastConfiguration const* gplast_config;

        Bool_t fOnline;

        TClonesArray* fcal_data;
        TClonesArray* funcal_data;
        TClonesArray* ftime_machine_array;

        gPlastTwinpeaksData* funcal_hit;
        gPlastTwinpeaksData* funcal_hit_next; // this is the slow or fast hit corresponding to the fast or slow hit :)
        gPlastTwinpeaksCalData* fcal_hit;

        uint16_t detector_id;
        
        double slow_lead_time;
        double slow_trail_time;

        double fast_lead_time;
        double fast_trail_time;

        double fast_ToT;
        double slow_ToT;

        int fNunmatched = 0;

        EventHeader * header;
        Int_t fNEvents = 0;
        Int_t fExecs = 0;
        int total_time_microsecs = 0;

        //internal status flags for detector map and calibration map:
        Bool_t DetectorMap_loaded = 0;
        Bool_t DetectorCal_loaded = 0;

        // time machine variables:
        int tm_delayed;
        int tm_undelayed;
        int time_machine_delayed_detector_id;
        int time_machine_undelayed_detector_id;
        //maps:
        std::map<std::pair<int,int>, int> detector_mapping; // [board_id][channel_id] -> [detector_id] (we do not care for position here)
        std::map<int,std::pair<double,double>> calibration_coeffs; // key: [detector id] -> [a0][a1]

    public:
        ClassDef(gPlastRaw2Cal, 1);
};

#endif
