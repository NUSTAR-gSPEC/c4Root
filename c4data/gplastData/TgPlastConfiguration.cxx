#include "TgPlastConfiguration.h"

#include "c4Logger.h"

#include <iostream>
#include <sstream>
#include <string>
#include <set>

#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"

TgPlastConfiguration* TgPlastConfiguration::instance = nullptr;
std::string TgPlastConfiguration::filepath = "gplast_alloc.txt";

std::map<TgPlastConfiguration::Position,std::string> TgPlastConfiguration::PosToString = {
	{TgPlastConfiguration::EXTRA, "Extra"},
	{TgPlastConfiguration::TOP, "Top"},
	{TgPlastConfiguration::BOTTOM, "Bottom"},
	{TgPlastConfiguration::LEFT, "Left"},
	{TgPlastConfiguration::RIGHT, "Right"},
	{TgPlastConfiguration::NB_POS, "<ERROR>"},
};


TgPlastConfiguration::TgPlastConfiguration()
    :   total_detectors(0),   num_tamex_boards(0)
{
    num_detectors[TOP] = 0;
    num_detectors[BOTTOM] = 0;
    num_detectors[LEFT] = 0;
    num_detectors[RIGHT] = 0;
    num_detectors[EXTRA] = 0;
    sc41l_d = sc41r_d = sc42l_d = sc42r_d = sc43l_d = sc43r_d = -1;
    frs_accept = -1;
    tm_delayed = tm_undelayed = -1;
    ReadConfiguration();
}

void TgPlastConfiguration::ReadConfiguration()
{
	c4LOG(info,"Reading gPlast configuration file in " + filepath);

    std::ifstream detector_map_file(filepath);
    std::string line;
    std::set<int> tamex_boards;
    std::set<int> detector_ids;
	std::set<std::pair<Position,int>> all_positions;

	const int NC = 4;
	
    if (detector_map_file.fail()) c4LOG(fatal, "Could not open gPlast Twinpeaks allocation");

    int nline = 0;
    while (std::getline(detector_map_file, line)){
		nline++;
		TString str(line);
		str.Strip(TString::kLeading,' '); // remove leading spaces
		if (str.Length()==0 || str[0] == '#') continue;

		TObjArray *arr = str.Tokenize("\t "); // accept whitespaces and tab
		Int_t nb_cols = arr->GetEntries();
		if (nb_cols < NC){
			c4LOG(fatal, Form("Bad mapping format in file %s (line %d)", filepath.data(), nline));  
		}
		else{
			int tamex_board = -1, tamex_channel = -1, detector_id = -1;
			Position pos = EXTRA;
			int pos_id = -1;
			for (int k=0;k<std::min(nb_cols,NC+1);k++){
				TString signal = ((TObjString*) arr->At(k))->GetString();
				switch (k) {
				case 0:
					if(signal.IsDigit()) tamex_board = signal.Atoi();
					else c4LOG(error, Form("Bad format in the mapping file (line %d) : col 1 (tamex board) should be an integer",nline));
					// consistency between data and mapping will be checked later on (gPlastRaw2Cal)
					break;
				case 1:
					if(signal.IsDigit()) tamex_channel = signal.Atoi();
					else c4LOG(error, Form("Bad format in the mapping file (line %d) : col 2 (tamex channel) should be an integer",nline));
					c4LOG_IF(error, tamex_channel > 16, Form("Line %d : tamex card does not have more than 16 channels",nline)); // hard coded --> should be a global constexpr defined elsewhere ?
					break;
				case 2:
					if(signal.IsDigit()) detector_id = signal.Atoi();
					else c4LOG(error, Form("Bad format in the mapping file (line %d) : col 3 (detector id) should be an integer",nline));
					break;
				case 3:
					signal.ToUpper();
					switch (signal[0]){
					case 'T': pos = TOP; break;
					case 'B': pos = BOTTOM; break;
					case 'L': pos = LEFT; break;
					case 'R': pos = RIGHT; break;
					default: 
						// some additional signal (pos == EXTRA)
						if (signal == "TIMEMACHINEU") tm_undelayed = detector_id;
						else if (signal == "TIMEMACHINED") tm_delayed = detector_id;
						else if (signal == "SC41L_D") sc41l_d = detector_id;
						else if (signal == "SC41R_D") sc41r_d = detector_id;
						else if (signal == "SC42L_D") sc42l_d = detector_id;
						else if (signal == "SC42R_D") sc42r_d = detector_id;
						else if (signal == "SC43L_D") sc43l_d = detector_id;
						else if (signal == "SC43R_D") sc43r_d = detector_id;
						else if (signal == "FRS_ACCEPT") frs_accept = detector_id;
						else c4LOG(warn, Form("Unknown auxiliary channel at line %d (detector %d) will be ignored",nline, detector_id));
					}
					if (pos != EXTRA){
						signal.Remove(0,1);
						if(signal.IsDigit()) pos_id = signal.Atoi();
						else c4LOG(error, Form("Bad format in the mapping file (line %d) : col 4 (position) should contain an integer",nline));
					}
					break;
				default:
					// extra column, might be a comment
					// if not, print a warning
					c4LOG_IF(warn, signal[0] != '#' , Form("Extra column at line %d (detector %d) is ignored",nline, detector_id));
				}
			} // end loop on columns
			if (detector_id > -1){
				// do mapping + extra checks for indexes appearing twice
				detectors[pos].insert(detector_id);
				std::pair<int, int> tamex_mc = {tamex_board, tamex_channel};
				auto map_ret = detector_mapping.insert(std::make_pair(tamex_mc, detector_id));
				c4LOG_IF(warn, !map_ret.second, Form("TAMEX channel in line %d has already been registered by detector %d",nline, map_ret.first->second));

				// positioning
				std::pair<Position, int> pos_det = {pos, pos_id};
				auto det_ret = position_mapping.insert(std::make_pair(detector_id, pos_det));
				c4LOG_IF(warn, !det_ret.second, Form("Detector %d (line %d) has already been registered",detector_id,nline));
				
				if (pos != EXTRA){
					auto pos_ret = all_positions.insert(pos_det);
					c4LOG_IF(warn, !pos_ret.second, Form("Detector position in line %d has already been registered",nline));
				}
			}
			if (tamex_board > -1){
				tamex_boards.insert(tamex_board);
			}
		}
    } // end loop on lines

    num_tamex_boards = tamex_boards.size();

	for (int k=0;k<NB_POS;k++){
		Position p = static_cast<Position>(k);
		int ndet = detectors[p].size();
		num_detectors[p] = ndet;
		total_detectors += ndet;
	}

	c4LOG_IF(error, total_detectors != (int) position_mapping.size(),"Something is not correct with the mapping");
	// this could happen when the same detector id is set to different positions --> see warnings that have probably been displayed beforehand

	// detector id should not be larger than the number of detectors 
	for (const auto& det: detector_mapping){
		c4LOG_IF(error, det.second > total_detectors, Form("Found detector with id=%d, but there are only %d channels in the mapping", det.second, total_detectors));
	}

	// // Online histograms assume idet in range [1,nbdet]
	// int nb_true_signals = total_detectors - num_detectors[EXTRA];
	// for (const auto& det: position_mapping){
	// 	c4LOG_IF(error, det.first > nb_true_signals || det.first < 1, Form("gPlast signals should be in the first %d channels for online histograms (found id=%d)", nb_true_signals, det.first));
	// }

    DetectorMap_loaded = 1;
    detector_map_file.close();
    return;
}
