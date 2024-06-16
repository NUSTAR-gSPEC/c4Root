#ifndef gPlastTwinpeaksCalData_H
#define gPlastTwinpeaksCalData_H

#include "TObject.h"

class gPlastTwinpeaksCalData : public TObject
{
    public:
        // Default Constructor
        gPlastTwinpeaksCalData();

        gPlastTwinpeaksCalData(
        uint16_t board_id,
        int ch_ID, // uint16_t
        uint16_t detector_id,
        
        double slow_lead_time,
        double slow_trail_time,

        double fast_lead_time,
        double fast_trail_time,

        double fast_ToT,
        double slow_ToT,
        uint16_t wr_subsystem_id,
        uint64_t wr_t);

        // Destructor
        virtual ~gPlastTwinpeaksCalData() {}

        // Getters

        inline const uint16_t Get_board_id() const {return fboard_id; }
        inline const uint16_t Get_ch_ID() const {return fch_ID; }
        inline const uint16_t Get_detector_id() const {return fdetector_id; }
        inline const double Get_slow_lead_time() const {return fslow_lead_time; }
        inline const double Get_slow_trail_time() const {return fslow_trail_time; }
        inline const double Get_fast_lead_time() const {return ffast_lead_time; }
        inline const double Get_fast_trail_time() const {return ffast_trail_time; }
        inline const double Get_fast_ToT() const {return ffast_ToT; }
        inline const double Get_slow_ToT() const {return fslow_ToT; }
        inline const uint16_t Get_wr_subsystem_id() const {return fwr_subsystem_id; }
        inline const uint64_t Get_wr_t() const {return fwr_t; }


        // Setters
        void Set_board_id(uint16_t v){fboard_id = v; }
        void Set_ch_ID(uint16_t v){ fch_ID = v; }
        void Set_detector_id(uint16_t v){ fdetector_id = v; }
        void Set_slow_lead_time(double v){ fslow_lead_time = v; }
        void Set_slow_trail_time(double v){ fslow_trail_time = v; }
        void Set_fast_lead_time(double v){ ffast_lead_time = v; }
        void Set_fast_trail_time(double v){ ffast_trail_time = v; }
        void Set_fast_ToT(double v){ ffast_ToT = v; }
        void Set_slow_ToT(double v){ fslow_ToT = v; }
        void Set_wr_subsystem_id(uint16_t v){ fwr_subsystem_id = v; }
        void Set_wr_t(uint64_t v){ fwr_t = v; }

    protected:
        // Data items
        uint16_t fboard_id;
        uint16_t fch_ID;
        uint16_t fdetector_id;
        
        double fslow_lead_time;
        double fslow_trail_time;

        double ffast_lead_time;
        double ffast_trail_time;

        double ffast_ToT;
        double fslow_ToT;
        uint16_t fwr_subsystem_id;
        uint64_t fwr_t;    

    public:
        ClassDef(gPlastTwinpeaksCalData, 1);
};

#endif