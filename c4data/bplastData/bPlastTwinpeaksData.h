#ifndef bPlastTwinpeaksData_H
#define bPlastTwinpeaksData_H

#include "TObject.h"

class bPlastTwinpeaksData : public TObject
{
    public:
        // Default Constructor
        bPlastTwinpeaksData();

        bPlastTwinpeaksData(
        uint16_t board_id,
        uint32_t ch_ID,
        uint32_t lead_epoch_counter,
        uint32_t lead_coarse_T,
        double lead_fine_T,

        uint32_t trail_epoch_counter,
        uint32_t trail_coarse_T,
        double trail_fine_T,
        uint16_t wr_subsystem_id,
        uint64_t wr_t);

        // Destructor
        virtual ~bPlastTwinpeaksData() {}

        // Getters

        inline const uint16_t Get_board_id() const {return fboard_id; }
        inline const uint32_t Get_ch_ID() const {return fch_ID; }

        inline const uint32_t Get_lead_epoch_counter() const {return flead_epoch_counter; }
        inline const uint32_t Get_lead_coarse_T() const {return flead_coarse_T; }
        inline const double Get_lead_fine_T() const {return flead_fine_T; }
    
        inline const uint32_t Get_trail_epoch_counter() const {return ftrail_epoch_counter; }
        inline const uint32_t Get_trail_coarse_T() const {return ftrail_coarse_T; }
        inline const double Get_trail_fine_T() const {return ftrail_fine_T; }


        inline const uint16_t Get_wr_subsystem_id() const { return fwr_subsystem_id; }
        inline const uint64_t Get_wr_t() const { return fwr_t; }



        // Setters
        void Set_board_id(uint16_t v){fboard_id = v;}
        void Set_ch_ID(uint32_t v){fch_ID = v;}
        void Set_lead_epoch_counter(uint32_t v){flead_epoch_counter = v;}
        void Set_lead_coarse_T(uint32_t v){flead_coarse_T = v;}
        void Set_lead_fine_T(double v){flead_fine_T = v;}
        void Set_trail_epoch_counter(uint32_t v){ftrail_epoch_counter = v;}
        void Set_trail_coarse_T(uint32_t v){ftrail_coarse_T = v;}
        void Set_trail_fine_T(double v){ftrail_fine_T = v;}
        void Set_wr_subsystem_id(uint32_t v) { fwr_subsystem_id = v; }
        void Set_wr_t(uint64_t v) { fwr_t = v; }

    protected:
        // Data items

        uint16_t fboard_id;
        uint32_t fch_ID;

        uint32_t flead_epoch_counter;
        uint32_t flead_coarse_T;
        double flead_fine_T;

        //if these guys are zero, no trail was found:
        uint32_t ftrail_epoch_counter;
        uint32_t ftrail_coarse_T;
        double ftrail_fine_T;

        //whiterabbit
        uint32_t fwr_subsystem_id;
        uint64_t fwr_t;
    

    public:
        ClassDef(bPlastTwinpeaksData, 1);
};

#endif
