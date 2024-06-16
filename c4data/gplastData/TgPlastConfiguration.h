#ifndef TgPlastConfiguration_H
#define TgPlastConfiguration_H

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <set>

class TgPlastConfiguration
{
public:
    enum Position
    {
        EXTRA = 0,
        TOP = 1,
        BOTTOM = 2,
        LEFT = 3,
        RIGHT = 4,
        NB_POS = 5
    };

public:
    static TgPlastConfiguration const *GetInstance();
    static void Create();
    static void SetDetectorMapFile(std::string fp) { filepath = fp; }

    std::map<std::pair<int, int>, int> Mapping() const;
    std::map<int,std::pair<Position, int>> Positioning() const;
    int NDetectors() const;
    int NTopDetectors() const;
    std::set<int> TopDetectors() const;
    int NBottomDetectors() const;
    std::set<int> BottomDetectors() const;
    int NLeftDetectors() const;
    std::set<int> LeftDetectors() const;
    int NRightDetectors() const;
    std::set<int> RightDetectors() const;
    int NTamexBoards() const;
    bool MappingLoaded() const;
    int TM_Undelayed() const;
    int TM_Delayed() const;
    int SC41L() const;
    int SC41R() const;
    int SC42L() const;
    int SC42R() const;
    int SC43L() const;
    int SC43R() const;
    int FRS_ACCEPT() const;
    int NExtraSignals() const;
    std::set<int> ExtraSignals() const;

    static std::string PositionToString(Position);

private:
    static std::string filepath;
    TgPlastConfiguration();
    void ReadConfiguration();

    static TgPlastConfiguration *instance;

    std::map<std::pair<int, int>, int> detector_mapping; // TAMEX --> det id
    std::map<int, std::pair<Position,int>> position_mapping; // det id --> position
    std::map<Position, std::set<int>> detectors;
    // std::set<int> top_detectors;
    // std::set<int> bottom_detectors;
    // std::set<int> left_detectors;
    // std::set<int> right_detectors;
    // std::set<int> extra_signals;

    int total_detectors;
    int num_tamex_boards;
    std::map<Position, int> num_detectors;
    // int num_top_detectors;
    // int num_bottom_detectors;
    // int num_left_detectors;
    // int num_right_detectors;

    int tm_undelayed;
    int tm_delayed;
    int sc41l_d;
    int sc41r_d;
    int sc42l_d;
    int sc42r_d;
    int sc43l_d;
    int sc43r_d;
    int frs_accept;

    bool DetectorMap_loaded = 0;

    static std::map<Position,std::string> PosToString;
};

inline TgPlastConfiguration const *TgPlastConfiguration::GetInstance()
{
    if (!instance)
    {
        TgPlastConfiguration::Create();
    }
    return instance;
}

inline void TgPlastConfiguration::Create()
{
    delete instance;
    instance = new TgPlastConfiguration();
}

inline std::map<std::pair<int, int>, int> TgPlastConfiguration::Mapping() const
{
    return detector_mapping;
}

inline std::map<int, std::pair<TgPlastConfiguration::Position, int>> TgPlastConfiguration::Positioning() const
{
    return position_mapping;
}

inline bool TgPlastConfiguration::MappingLoaded() const
{
    return DetectorMap_loaded;
}

inline int TgPlastConfiguration::NDetectors() const
{
    return total_detectors;
}

inline int TgPlastConfiguration::NTamexBoards() const
{
    return num_tamex_boards;
}

inline int TgPlastConfiguration::NTopDetectors() const
{
    return num_detectors.at(TOP);
}

inline std::set<int> TgPlastConfiguration::TopDetectors() const
{
    return detectors.at(TOP);
}

inline int TgPlastConfiguration::NBottomDetectors() const
{
    return num_detectors.at(BOTTOM);
}

inline std::set<int> TgPlastConfiguration::BottomDetectors() const
{
    return detectors.at(BOTTOM);
}

inline int TgPlastConfiguration::NLeftDetectors() const
{
    return num_detectors.at(LEFT);
}

inline std::set<int> TgPlastConfiguration::LeftDetectors() const
{
    return detectors.at(LEFT);
}

inline int TgPlastConfiguration::NRightDetectors() const
{
    return num_detectors.at(RIGHT);
}

inline std::set<int> TgPlastConfiguration::RightDetectors() const
{
    return detectors.at(RIGHT);
}

inline int TgPlastConfiguration::TM_Undelayed() const
{
    return tm_undelayed;
}

inline int TgPlastConfiguration::TM_Delayed() const
{
    return tm_delayed;
}

inline int TgPlastConfiguration::SC41L() const
{
    return sc41l_d;
}

inline int TgPlastConfiguration::SC41R() const
{
    return sc41r_d;
}

inline int TgPlastConfiguration::SC42L() const
{
    return sc42l_d;
}

inline int TgPlastConfiguration::SC42R() const
{
    return sc42r_d;
}

inline int TgPlastConfiguration::SC43L() const
{
    return sc43l_d;
}

inline int TgPlastConfiguration::SC43R() const
{
    return sc43r_d;
}

inline int TgPlastConfiguration::FRS_ACCEPT() const
{
    return frs_accept;
}

inline int TgPlastConfiguration::NExtraSignals() const
{
    return num_detectors.at(EXTRA);
}

inline std::set<int> TgPlastConfiguration::ExtraSignals() const
{
    return detectors.at(EXTRA);
}

inline std::string TgPlastConfiguration::PositionToString(Position p)
{
    return PosToString.at(p);
}

#endif
