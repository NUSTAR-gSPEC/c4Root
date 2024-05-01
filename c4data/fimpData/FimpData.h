#ifndef FimpData_H
#define FimpData_H

#include "TObject.h"

class FimpItem : public TObject
{
    public:
        FimpItem();

        virtual ~FimpItem() {}

        void SetAll(uint64_t wr, 
                    uint16_t id,
                    uint64_t trig_t, 
                    uint16_t chan, 
                    uint16_t lead_ct, 
                    uint16_t lead_ft, 
                    uint16_t trail_ct, 
                    uint16_t trail_ft);
        void Reset();

        uint64_t Get_wr_t() const;
        uint16_t Get_wr_subsystem_id() const;
        uint64_t Get_trig_time_long() const;
        uint16_t Get_channel() const;
        uint16_t Get_lead_coarse_time() const;
        uint16_t Get_lead_fine_time() const;
        uint16_t Get_trail_coarse_time() const;
        uint16_t Get_trail_fine_time() const;

        ClassDef(FimpItem, 2);
    
    private:
        // members
        uint64_t wr_t;
        uint16_t wr_id;
        uint64_t trig_time_long; 
        uint16_t channel;
        uint16_t lead_coarse_time;
        uint16_t lead_fine_time;
        uint16_t trail_coarse_time;
        uint16_t trail_fine_time;
        bool leadOrTrail;
};

inline uint64_t FimpItem::Get_wr_t() const
{
    return wr_t;
}

inline uint16_t FimpItem::Get_wr_subsystem_id() const
{
    return wr_id;
}

inline uint64_t FimpItem::Get_trig_time_long() const
{
    return trig_time_long;
}

inline uint16_t FimpItem::Get_channel() const
{
    return channel;
}

inline uint16_t FimpItem::Get_lead_coarse_time() const
{
    return lead_coarse_time;
}

inline uint16_t FimpItem::Get_lead_fine_time() const
{
    return lead_fine_time;
}

inline uint16_t FimpItem::Get_trail_coarse_time() const
{
    return trail_coarse_time;
}

inline uint16_t FimpItem::Get_trail_fine_time() const
{
    return trail_fine_time;
}


#endif