#ifndef FrsHitData_H
#define FrsHitData_H

#include "../../config/setup.h"
#include <vector>
#include "TObject.h"

class FrsHitData : public TObject
{
    public:
        FrsHitData();

        ~FrsHitData();

        inline const uint64_t Get_wr_t() const { return fwr_t; }
        inline const uint16_t Get_tpat() const { return ftpat; }
        
        inline const Float_t Get_ID_x2() const { return fID_x2; }
        inline const Float_t Get_ID_y2() const { return fID_x2; }
        inline const Float_t Get_ID_a2() const { return fID_a2; }
        inline const Float_t Get_ID_b2() const { return fID_b2; }
        inline const Float_t Get_ID_x4() const { return fID_x4; }
        inline const Float_t Get_ID_y4() const { return fID_y4; }
        inline const Float_t Get_ID_a4() const { return fID_a4; }
        inline const Float_t Get_ID_b4() const { return fID_b4; }
        inline const Float_t Get_ID_AoQ() const { return fID_AoQ; }
        inline const Float_t Get_ID_AoQ_corr() const { return fID_AoQ_corr; }
        inline const Float_t Get_ID_z() const { return fID_z; }
        inline const Float_t Get_ID_z2() const { return fID_z2; }
        inline const Float_t Get_ID_beta() const { return fID_beta; }
        inline const Float_t Get_ID_dEdegoQ() const { return fID_dEdegoQ; }
        inline const Float_t Get_ID_dEdeg() const { return fID_dEdeg; }

        inline const Float_t Get_music_dE(int channel) const { return fmusic_dE[channel]; }
        inline const Float_t Get_sci_e(int channel) const { return fsci_e[channel]; }
        inline const Float_t Get_sci_tof2() const { return fsci_tof2; }
        inline const Float_t Get_sci_l(int channel) const { return fsci_l[channel]; }
        inline const Float_t Get_sci_r(int channel) const { return fsci_r[channel]; }
        inline const Float_t Get_sci_tx(int channel) const { return fsci_tx[channel]; }

        inline const Float_t Get_ID_z_mhtdc(int mhit) const 
        {
            if (mhit < fID_z_mhtdc.size())
            {
                return fID_z_mhtdc[mhit];
            }
            else
            {
                return 0.0;
            }
        }
        inline const Float_t Get_ID_z2_mhtdc(int mhit) const 
        {
            if (mhit < fID_z2_mhtdc.size())
            {
                return fID_z2_mhtdc[mhit];
            }
            else
            {
                return 0.0;
            }
        }
        inline const Float_t Get_ID_AoQ_mhtdc(int mhit) const 
        { 
            if (mhit < fID_AoQ_mhtdc.size())
            {
                return fID_AoQ_mhtdc[mhit];
            }
            else
            {
                return 0.0;
            } 
        }
        inline const Float_t Get_ID_AoQ_corr_mhtdc(int mhit) const 
        {
            if (mhit < fID_AoQ_corr_mhtdc.size())
            {
                return fID_AoQ_corr_mhtdc[mhit];
            }
            else
            {
                return 0.0;
            } 
        }
        inline const Float_t Get_ID_dEdeg_mhtdc(int mhit) const 
        {
            if (mhit < fID_dEdeg_mhtdc.size())
            {
                return fID_dEdeg_mhtdc[mhit];
            }
            else
            {
                return 0.0;
            } 
        }
        inline const Float_t Get_ID_dEdegoQ_mhtdc(int mhit) const 
        { 
            if (mhit < fID_dEdegoQ_mhtdc.size())
            {
                return fID_dEdegoQ_mhtdc[mhit];
            }
            else
            {
                return 0.0;
            } 
        }

        inline uint32_t Get_time_in_ms() { return ftime_in_ms; }
        inline uint32_t Get_ibin_for_s() { return fibin_for_s; }
        inline uint32_t Get_ibin_for_100ms() { return fibin_for_100ms; }
        inline uint32_t Get_ibin_for_spill() { return fibin_for_spill; }
        inline uint32_t* Get_increase_sc_temp_main() { return fincrease_sc_temp_main; }
        inline uint32_t* Get_increase_sc_temp_user() { return fincrease_sc_temp_user; }
        inline uint32_t Get_increase_sc_temp2() { return fincrease_sc_temp2; }
        inline uint32_t Get_increase_sc_temp3() { return fincrease_sc_temp3; }
        inline uint32_t Get_extraction_time_ms() { return fextraction_time_ms; }
        inline uint32_t Get_ibin_clean_for_s() { return fibin_clean_for_s; }
        inline uint32_t Get_ibin_clean_for_100ms() { return fibin_clean_for_100ms; }
        inline uint32_t Get_ibin_clean_for_spill() { return fibin_clean_for_spill; }


        // ----- Setter -------------------------------
        inline void Set_wr_t(ULong64_t ts) { fwr_t = ts; }
        inline void Set_tpat(uint16_t tpat) { ftpat = tpat; }
        
        inline void Set_time_in_ms(uint32_t time) { ftime_in_ms = time; }
        inline void Set_ibin_for_s(uint32_t ibin) { fibin_for_s = ibin; }
        inline void Set_ibin_for_100ms(uint32_t ibin) { fibin_for_100ms = ibin; }
        inline void Set_ibin_for_spill(uint32_t ibin) { fibin_for_spill = ibin; }
        inline void Set_increase_sc_temp_user(int index, uint32_t increase) { fincrease_sc_temp_user[index] = increase; }
        inline void Set_increase_sc_temp_main(int index, uint32_t increase) { fincrease_sc_temp_main[index] = increase; }
        inline void Set_increase_sc_temp2(uint32_t increase) { fincrease_sc_temp2 = increase; }
        inline void Set_increase_sc_temp3(uint32_t increase) { fincrease_sc_temp3 = increase; }
        inline void Set_extraction_time_ms(uint32_t time) { fextraction_time_ms = time; }
        inline void Set_ibin_clean_for_s(uint32_t ibin) { fibin_clean_for_s = ibin; }
        inline void Set_ibin_clean_for_100ms(uint32_t ibin) { fibin_clean_for_100ms = ibin; }
        inline void Set_ibin_clean_for_spill(uint32_t ibin) { fibin_clean_for_spill = ibin; }

        inline void Set_music_dE(int index, Float_t dE) { fmusic_dE[index] = dE; }
        inline void Set_music_dE_cor(int index, Float_t dE_cor) { fmusic_dE_cor[index] = dE_cor; }
        inline void Set_sci_l(int index, Float_t sci_l) { fsci_l[index] = sci_l; }
        inline void Set_sci_r(int index, Float_t sci_r) { fsci_r[index] = sci_r; }
        inline void Set_sci_tx(int index, Float_t sci_tx) { fsci_tx[index] = sci_tx; }
        inline void Set_sci_e(int index, Float_t sci_e) { fsci_e[index] = sci_e; }
        inline void Set_sci_tof2(Float_t sci_tof2) { fsci_tof2 = sci_tof2; }
        inline void Set_sci_tof(int index, Float_t sci_tof) { fsci_tof[index] = sci_tof; }
        inline void Set_sci_tof_calib(int index, Float_t sci_tof_calib) { fsci_tof_calib[index] = sci_tof_calib; }

        inline void Set_ID_beta(Float_t beta) { fID_beta = beta; }
        inline void Set_ID_AoQ(Float_t aoq) { fID_AoQ = aoq; }
        inline void Set_ID_AoQ_corr(Float_t aoq_corr) { fID_AoQ_corr = aoq_corr; }
        inline void Set_ID_z(Float_t z) { fID_z = z; }
        inline void Set_ID_z2(Float_t z2) { fID_z2 = z2; }
        inline void Set_ID_x2(Float_t x2) { fID_x2 = x2; }
        inline void Set_ID_y2(Float_t y2) { fID_y2 = y2; }
        inline void Set_ID_a2(Float_t a2) { fID_a2 = a2; }
        inline void Set_ID_b2(Float_t b2) { fID_b2 = b2; }
        inline void Set_ID_x4(Float_t x4) { fID_x4 = x4; }
        inline void Set_ID_y4(Float_t y4) { fID_y4 = y4; }
        inline void Set_ID_a4(Float_t a4) { fID_a4 = a4; }
        inline void Set_ID_b4(Float_t b4) { fID_b4 = b4; }
        inline void Set_ID_dEdegoQ(Float_t dEdegoQ) { fID_dEdegoQ = dEdegoQ; }
        inline void Set_ID_dEdeg(Float_t dEdeg) { fID_dEdeg = dEdeg; }
        inline void Set_ID_rho(int index, Float_t rho) { fID_rho[index] = rho; }
        inline void Set_ID_brho(int index, Float_t brho) { fID_brho[index] = brho; }

        inline void Set_ID_beta_mhtdc(std::vector<Float_t> beta) { fID_beta_mhtdc = beta; }
        inline void Set_ID_AoQ_mhtdc(std::vector<Float_t> aoq) { fID_AoQ_mhtdc = aoq; }
        inline void Set_ID_AoQ_corr_mhtdc(std::vector<Float_t> aoq_corr) { fID_AoQ_corr_mhtdc = aoq_corr; }
        inline void Set_ID_z_mhtdc(std::vector<Float_t> z) { fID_z_mhtdc = z; }
        inline void Set_ID_z2_mhtdc(std::vector<Float_t> z2) { fID_z2_mhtdc = z2; }
        inline void Set_ID_dEdegoQ_mhtdc(std::vector<Float_t> dEdegoQ) { fID_dEdegoQ_mhtdc = dEdegoQ; }
        inline void Set_ID_dEdeg_mhtdc(std::vector<Float_t> dEdeg) { fID_dEdeg_mhtdc = dEdeg; }


    protected:

        uint64_t fwr_t;
        uint16_t ftpat;
        Float_t fID_x2;
        Float_t fID_y2;
        Float_t fID_x4;
        Float_t fID_y4;
        Float_t fID_a2;
        Float_t fID_b2;
        Float_t fID_a4;
        Float_t fID_b4;
        Float_t fID_AoQ;
        Float_t fID_AoQ_corr;
        Float_t fID_z;
        Float_t fID_z2;
        Float_t fID_beta;
        Float_t fID_dEdegoQ;
        Float_t fID_dEdeg;
        Float_t fID_rho[2];
        Float_t fID_brho[2];
       
        Float_t fmusic_dE[2]; // CEJ: 2? 3 for more musics?
        Float_t fmusic_dE_cor[2];
        Float_t fsci_e[6]; // CEJ: by chance we need [6]. 5->2 but 10->5 (so 6 elements)
        Float_t fsci_l[6];
        Float_t fsci_r[6];
        Float_t fsci_tx[6];
        Float_t fsci_tof2;
        Float_t fsci_tof[6];
        Float_t fsci_tof_calib[6];

        // MHTDC
        std::vector<Float_t> fID_beta_mhtdc;
        std::vector<Float_t> fID_AoQ_mhtdc;
        std::vector<Float_t> fID_AoQ_corr_mhtdc;
        std::vector<Float_t> fID_z_mhtdc;
        std::vector<Float_t> fID_z2_mhtdc;
        std::vector<Float_t> fID_dEdeg_mhtdc;
        std::vector<Float_t> fID_dEdegoQ_mhtdc;

        uint32_t ftime_in_ms;
        uint32_t fibin_for_s;
        uint32_t fibin_for_100ms;
        uint32_t fibin_for_spill;
        uint32_t fincrease_sc_temp_main[32];
        uint32_t fincrease_sc_temp_user[32];
        uint32_t fincrease_sc_temp2;
        uint32_t fincrease_sc_temp3;
        uint32_t fextraction_time_ms;
        uint32_t fibin_clean_for_s;
        uint32_t fibin_clean_for_100ms;
        uint32_t fibin_clean_for_spill;
    
    public:
        ClassDef(FrsHitData, 1)

};


class FrsHitItem : public TObject
{
    public:
        FrsHitItem();

        void SetAll(uint64_t wr_t,
                    uint16_t tpat,
                    Float_t x2,
                    Float_t y2,
                    Float_t x4,
                    Float_t y4,
                    Float_t a2,
                    Float_t b2,
                    Float_t a4,
                    Float_t b4,
                    Float_t AoQ,
                    Float_t AoQ_corr,
                    Float_t z,
                    Float_t z2,
                    Float_t beta,
                    Float_t dEdeg,
                    Float_t dEdegoQ,
                    Float_t* rho,
                    Float_t* brho,
                    Float_t* music_dE,
                    Float_t* music_dE_cor,
                    Float_t* sci_e,
                    Float_t* sci_l,
                    Float_t* sci_r,
                    Float_t* sci_x,
                    Float_t* sci_tof,
                    Float_t* sci_tof_calib,
                    Float_t sci_tof2,
                    uint32_t time_in_ms,
                    uint32_t ibin_for_s,
                    uint32_t ibin_for_100ms,
                    uint32_t ibin_for_spill,
                    uint32_t* increase_sc_temp_main,
                    uint32_t* increase_sc_temp_user,
                    uint32_t increase_sc_temp2,
                    uint32_t increase_sc_temp3,
                    uint32_t extraction_time_ms,
                    uint32_t ibin_clean_for_s,
                    uint32_t ibin_clean_for_100ms,
                    uint32_t ibin_clean_for_spill);
        void Reset();

        uint64_t Get_wr_t() const;
        uint16_t Get_tpat() const;
        Float_t Get_ID_x2() const;
        Float_t Get_ID_y2() const;
        Float_t Get_ID_x4() const;
        Float_t Get_ID_y4() const;
        Float_t Get_ID_a2() const;
        Float_t Get_ID_b2() const;
        Float_t Get_ID_a4() const;
        Float_t Get_ID_b4() const;
        Float_t Get_ID_AoQ() const;
        Float_t Get_ID_AoQ_corr() const;
        Float_t Get_ID_z() const;
        Float_t Get_ID_z2() const;
        Float_t Get_ID_beta() const;
        Float_t Get_ID_dEdeg() const;
        Float_t Get_ID_dEdegoQ() const;
        Float_t* Get_ID_rho() const;
        Float_t* Get_ID_brho() const;
        Float_t* Get_music_dE() const;
        Float_t* Get_music_dE_cor() const;
        Float_t* Get_sci_e() const;
        Float_t* Get_sci_l() const;
        Float_t* Get_sci_r() const;
        Float_t* Get_sci_x() const;
        Float_t* Get_sci_tof() const;
        Float_t* Get_sci_tof_calib() const;
        Float_t Get_sci_tof2() const;
        uint32_t Get_time_in_ms() const;
        uint32_t Get_ibin_for_s() const;
        uint32_t Get_ibin_for_100ms() const;
        uint32_t Get_ibin_for_spill() const;
        uint32_t* Get_increase_sc_temp_main() const;
        uint32_t* Get_increase_sc_temp_user() const;
        uint32_t Get_increase_sc_temp2() const;
        uint32_t Get_increase_sc_temp3() const;
        uint32_t Get_extraction_time_ms() const;
        uint32_t Get_ibin_clean_for_s() const;
        uint32_t Get_ibin_clean_for_100ms() const;
        uint32_t Get_ibin_clean_for_spill() const;

        ClassDefNV(FrsHitItem, 2);

    private:

        uint64_t fwr_t;
        uint16_t ftpat;

        // ids
        Float_t fID_x2;
        Float_t fID_y2;
        Float_t fID_x4;
        Float_t fID_y4;
        Float_t fID_a2;
        Float_t fID_b2;
        Float_t fID_a4;
        Float_t fID_b4;
        Float_t fID_AoQ;
        Float_t fID_AoQ_corr;
        Float_t fID_z;
        Float_t fID_z2;
        Float_t fID_beta;
        Float_t fID_dEdeg;
        Float_t fID_dEdegoQ;
        Float_t* fID_rho; // [2]
        Float_t* fID_brho; // [2]

        // scitillators / energy loss
        Float_t* fmusic_dE; // [2] // CEJ: 2? 3 for more musics?
        Float_t* fmusic_dE_cor; // [2]
        Float_t* fsci_e; // [6] // CEJ: by chance we need [6]. 5->2 but 10->5 (so 6 elements)
        Float_t* fsci_l; // [6]
        Float_t* fsci_r; // [6]
        Float_t* fsci_x; // [6]
        Float_t fsci_tof2; // should we look at passing all the tof calcs?
        Float_t* fsci_tof; // [6]
        Float_t* fsci_tof_calib; // [6]

        // scalers
        uint32_t ftime_in_ms;
        uint32_t fibin_for_s;
        uint32_t fibin_for_100ms;
        uint32_t fibin_for_spill;
        uint32_t* fincrease_sc_temp_main; // [32]
        uint32_t* fincrease_sc_temp_user; // [32]
        uint32_t fincrease_sc_temp2;
        uint32_t fincrease_sc_temp3;
        uint32_t fextraction_time_ms;
        uint32_t fibin_clean_for_s;
        uint32_t fibin_clean_for_100ms;
        uint32_t fibin_clean_for_spill;

};

class FrsMultiHitItem : public TObject
{
    public: 
        FrsMultiHitItem();

        void SetAll(Float_t beta, 
                    Float_t AoQ, 
                    Float_t AoQ_corr, 
                    Float_t z,
                    Float_t z2,
                    Float_t dEdeg,
                    Float_t dEdegoQ);
        void Reset();

        Float_t Get_ID_beta_mhtdc() const;
        Float_t Get_ID_AoQ_mhtdc() const;
        Float_t Get_ID_AoQ_corr_mhtdc() const;
        Float_t Get_ID_z_mhtdc() const;
        Float_t Get_ID_z2_mhtdc() const;
        Float_t Get_ID_dEdeg_mhtdc() const;
        Float_t Get_ID_dEdegoQ_mhtdc() const;

        ClassDefNV(FrsMultiHitItem, 2);

    private:

        Float_t fID_beta_mhtdc;
        Float_t fID_AoQ_mhtdc;
        Float_t fID_AoQ_corr_mhtdc;
        Float_t fID_z_mhtdc;
        Float_t fID_z2_mhtdc;
        Float_t fID_dEdeg_mhtdc;
        Float_t fID_dEdegoQ_mhtdc;

};


inline uint64_t FrsHitItem::Get_wr_t() const
{
    return fwr_t;
}

inline uint16_t FrsHitItem::Get_tpat() const
{
    return ftpat;
}

inline Float_t FrsHitItem::Get_ID_x2() const
{
    return fID_x2;
}

inline Float_t FrsHitItem::Get_ID_y2() const
{
    return fID_y2;
}

inline Float_t FrsHitItem::Get_ID_x4() const
{
    return fID_x4;
}

inline Float_t FrsHitItem::Get_ID_y4() const
{
    return fID_y4;
}

inline Float_t FrsHitItem::Get_ID_a2() const
{
    return fID_a2;
}

inline Float_t FrsHitItem::Get_ID_b2() const
{
    return fID_b2;
}

inline Float_t FrsHitItem::Get_ID_a4() const
{
    return fID_a4;
}

inline Float_t FrsHitItem::Get_ID_b4() const
{
    return fID_b4;
}

inline Float_t FrsHitItem::Get_ID_AoQ() const
{
    return fID_AoQ;
}

inline Float_t FrsHitItem::Get_ID_AoQ_corr() const
{
    return fID_AoQ_corr;
}

inline Float_t FrsHitItem::Get_ID_z() const
{
    return fID_z;
}

inline Float_t FrsHitItem::Get_ID_z2() const
{
    return fID_z2;
}

inline Float_t FrsHitItem::Get_ID_beta() const
{
    return fID_beta;
}

inline Float_t FrsHitItem::Get_ID_dEdeg() const
{
    return fID_dEdeg;
}

inline Float_t FrsHitItem::Get_ID_dEdegoQ() const
{
    return fID_dEdegoQ;
}

inline Float_t* FrsHitItem::Get_ID_rho() const
{
    return fID_rho;
}

inline Float_t* FrsHitItem::Get_ID_brho() const
{
    return fID_brho;
}

inline Float_t* FrsHitItem::Get_music_dE() const
{
    return fmusic_dE;
}

inline Float_t* FrsHitItem::Get_music_dE_cor() const
{
    return fmusic_dE_cor;
}

inline Float_t* FrsHitItem::Get_sci_e() const
{
    return fsci_e;
}

inline Float_t* FrsHitItem::Get_sci_l() const
{
    return fsci_l;
}

inline Float_t* FrsHitItem::Get_sci_r() const
{
    return fsci_r;
}

inline Float_t* FrsHitItem::Get_sci_x() const
{
    return fsci_x;
}

inline Float_t* FrsHitItem::Get_sci_tof() const
{
    return fsci_tof;
}

inline Float_t* FrsHitItem::Get_sci_tof_calib() const
{
    return fsci_tof_calib;
}

inline Float_t FrsHitItem::Get_sci_tof2() const
{
    return fsci_tof2;
}

inline uint32_t FrsHitItem::Get_time_in_ms() const
{
    return ftime_in_ms;
}

inline uint32_t FrsHitItem::Get_ibin_for_s() const
{
    return fibin_for_s;
}

inline uint32_t FrsHitItem::Get_ibin_for_100ms() const
{
    return fibin_for_100ms;
}

inline uint32_t FrsHitItem::Get_ibin_for_spill() const
{
    return fibin_for_spill;
}

inline uint32_t* FrsHitItem::Get_increase_sc_temp_main() const
{
    return fincrease_sc_temp_main;
}

inline uint32_t* FrsHitItem::Get_increase_sc_temp_user() const
{
    return fincrease_sc_temp_user;
}

inline uint32_t FrsHitItem::Get_increase_sc_temp2() const
{
    return fincrease_sc_temp2;
}

inline uint32_t FrsHitItem::Get_increase_sc_temp3() const
{
    return fincrease_sc_temp3;
}

inline uint32_t FrsHitItem::Get_extraction_time_ms() const
{
    return fextraction_time_ms;
}
inline uint32_t FrsHitItem::Get_ibin_clean_for_s() const
{
    return fibin_clean_for_s;
}

inline uint32_t FrsHitItem::Get_ibin_clean_for_100ms() const
{
    return fibin_clean_for_100ms;
}

inline uint32_t FrsHitItem::Get_ibin_clean_for_spill() const
{
    return fibin_clean_for_spill;
}

inline Float_t FrsMultiHitItem::Get_ID_beta_mhtdc() const
{
    return fID_beta_mhtdc;
}

inline Float_t FrsMultiHitItem::Get_ID_AoQ_mhtdc() const
{
    return fID_AoQ_mhtdc;
}

inline Float_t FrsMultiHitItem::Get_ID_AoQ_corr_mhtdc() const
{
    return fID_AoQ_corr_mhtdc;
}

inline Float_t FrsMultiHitItem::Get_ID_z_mhtdc() const
{
    return fID_z_mhtdc;
}

inline Float_t FrsMultiHitItem::Get_ID_z2_mhtdc() const
{
    return fID_z2_mhtdc;
}

inline Float_t FrsMultiHitItem::Get_ID_dEdeg_mhtdc() const
{
    return fID_dEdeg_mhtdc;
}

inline Float_t FrsMultiHitItem::Get_ID_dEdegoQ_mhtdc() const
{
    return fID_dEdegoQ_mhtdc;
}


#endif
