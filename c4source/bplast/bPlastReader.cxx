// FairRoot
#include "FairLogger.h"
#include "FairRootManager.h"

#include "TH1.h"
#include "TFile.h"

#include <iomanip>

// c4
#include "bPlastTwinpeaksData.h"
#include "bPlastReader.h"
#include "c4Logger.h"

#include "TClonesArray.h"
#include "ext_data_struct_info.hh"

extern "C"
{
    #include "ext_data_client.h"
    #include "ext_h101_bplast.h"
}

/*
Constructor. Needs the input data structure obtained from ucesb with the following args:
 --ntuple=UNPACK:fatima,STRUCT_HH,ext_h101_fatima.h (add NOEVENTHEAD also). 

The resulting file must be in the same folder as this file.

*/

bPlastReader::bPlastReader(EXT_STR_h101_bplast_onion* data, size_t offset)

    : c4Reader("bPlastReader")
    , fNEvent(0)
    , fData(data)
    , fOffset(offset)
    , fOnline(kFALSE)
    , fArray(new TClonesArray("bPlastTwinpeaksData"))
{
}

/*
Deletes the arrays allocated.
And prints some statistics for the run.
*/
bPlastReader::~bPlastReader() { 

    if (fPrintStatistics == true) 
    {
        PrintStatistics();
    }

    for (int i = 0; i < NBoards; i++) {
        for (int j = 0; j < NChannels; j++) {
            if (fine_time_calibration_coeffs[i][j] != nullptr) {
                delete[] fine_time_calibration_coeffs[i][j];
                fine_time_calibration_coeffs[i][j] = nullptr;
            }
        }

        if (fine_time_calibration_coeffs[i] != nullptr) {
            delete[] fine_time_calibration_coeffs[i];
            fine_time_calibration_coeffs[i] = nullptr;
        }

        if (fine_time_hits[i] != nullptr) {
            delete[] fine_time_hits[i];
            fine_time_hits[i] = nullptr;
        }
    }

    if (fine_time_hits != nullptr) {
        delete[] fine_time_hits;
        fine_time_hits = nullptr;
    }

    if (fine_time_calibration_coeffs != nullptr) {
        delete[] fine_time_calibration_coeffs;
        fine_time_calibration_coeffs = nullptr;
    }

    if (fArray != nullptr) delete fArray;

    c4LOG(info, "Destroyed bPlastReader properly.");

}

/*
Required set up before reading. 

If the fine time calibration histograms are read from disk, it loads these.
If not new histograms are allocated and ready for filling and calibration.

*/
Bool_t bPlastReader::Init(ext_data_struct_info* a_struct_info)
{
    Int_t ok;
    c4LOG(info, "");

    EXT_STR_h101_bplast_ITEMS_INFO(ok, *a_struct_info, fOffset, EXT_STR_h101_bplast, 0);

    if (!ok)
    {
        c4LOG(error, "Failed to set up structure information");
        return kFALSE;
    }

    //allocated 3D array for the fine time calibration look-up table: fine_time_calibration_coeffs[board nr][channel nr][tdc value nr]
    fine_time_calibration_coeffs = new double**[NBoards];
    for (int i = 0; i < NBoards; i++) {
        fine_time_calibration_coeffs[i] = new double*[NChannels];    

    }

    // initialising the last hits to zero
    for (int i = 0; i < NBoards; i++) 
    {

        for (int j = 0; j < NChannels; j++) {
            fine_time_calibration_coeffs[i][j] = new double[Nbins_fine_time];
            for (int k = 0; k<Nbins_fine_time; k++) fine_time_calibration_coeffs[i][j][k] = 0.0;
            
            
            fNtrails_read[i][j] = 0;
            fNleads_read[i][j] = 0;
            fNmatched[i][j] = 0;
            fNevents_lacking_epoch[i][j] = 0;
            fNevents_TAMEX_fail[i][j] = 0;
            fNevents_second_lead_seen[i][j] = 0;
            fNevents_trail_seen_no_lead[i][j] = 0;
        }
    }


    // fine time calibration if specified by the user. otherwise make the fine time histograms.

    if (fine_time_calibration_read_from_file)
    {
        ReadFineTimeHistosFromFile();
        DoFineTimeCalibration();
        fine_time_calibration_set = true;
        c4LOG(info,"Fine Time calibration set from file.");
    }
    else
    {
        fine_time_hits = new TH1I**[NBoards];
        for (int i = 0; i < NBoards; i++) {
            fine_time_hits[i] = new TH1I*[NChannels];
            for (int j = 0; j<NChannels; j++){
                fine_time_hits[i][j] = new TH1I(Form("fine_time_hits_%i_%i", i, j),Form("fine_time_hits_%i_%i", i, j),Nbins_fine_time,0,Nbins_fine_time);
            }
        }
        fine_time_calibration_set = false;
    }

    // Register output array in a tree
    FairRootManager::Instance()->Register("bPlastTwinpeaksData", "bPlastTwinpeaksDataFolder", fArray, !fOnline);
    fArray->Clear();

    memset(fData, 0, sizeof *fData);

    c4LOG(info,"bPlastReader init setup.");

    return kTRUE;
}

/*
This does the fine time calibrations of the fine times and builds the required look-up table in fine_time_calibrations[board][channel][tdc channel] = time (ns)

This can be called explicitly if desired - but will be done automatically by the rest of the code. Throws a warning if there are channels without hits. These are mapped tdc_channel = tdc_channel/512*5 ns
*/
void bPlastReader::DoFineTimeCalibration(){
    c4LOG(info, "Doing fine time calibrations.");
    std::vector<std::pair<int>> warning_channels;
    int warning_counter = 0;
    for (int i = 0; i < NBoards; i++) {
        for (int j = 0; j < NChannels; j++) {
            int running_sum = 0;
            int total_counts = fine_time_hits[i][j]->GetEntries();
            if (total_counts == 0) 
            {
                std::pair<int, int> pair = std::make_pair(j, i); // channel, board
                warning_channels.emplace_back(pair); // dump to a log file in future
                warning_counter++;
                // c4LOG(warning,Form("Channel %i on board %i does not have any fine time hits in the interval.",j,i));
            
            }
            for (int k = 0; k < Nbins_fine_time; k++) {
                running_sum += fine_time_hits[i][j]->GetBinContent(k+1); //bin 0 is the underflow bin, hence we start at [1,Nbins_fine_time].
                //no counts?
                if (total_counts == 0) {
                    fine_time_calibration_coeffs[i][j][k] = k*TAMEX_fine_time_clock/(double)Nbins_fine_time;
                    continue;
                }

                fine_time_calibration_coeffs[i][j][k] = TAMEX_fine_time_clock*(double)running_sum/(double)total_counts;
            }
        }
    }

    if (warning_counter > 0) c4LOG(warning, Form("%i channels do not have any fine time hits in the interval.", warning_counter));
    fine_time_calibration_set = true;
}

/*
Uses the conversion table to look-up fine times in ns.
*/
double bPlastReader::GetFineTime(int tdc_fine_time_channel, int board_id, int channel_id){
    return fine_time_calibration_coeffs[board_id][channel_id][tdc_fine_time_channel];
}

/*
Fine time histograms are stored as ROOT TH1I histograms as they are efficient and compresses when written.
This function saves the fine time hits directly to file. On restart this file can then be read and the look-up table reconstructed.
*/
void bPlastReader::WriteFineTimeHistosToFile(){

    if (!fine_time_calibration_set) {
        c4LOG(info,"Fine time calibrations not set, cannot write to file.");
        return;
    }
    c4LOG(info,"Fine time calibrations starting to write to file.");

    TFile * outputfile = TFile::Open(Form("%s",fine_time_histo_outfile.Data()),"recreate");
    
    std::cout << outputfile << std::endl;
    
    if (!outputfile->IsOpen()) {c4LOG(warning, "File to write histos not opened, skipping writing. "); return;}
    // i think the easiest is to write the histograms of the fine times themselves - please let me know if you disagree (JEL)
    for (int i = 0; i < NBoards; i++) {
        for (int j = 0; j < NChannels; j++) {
            if (fine_time_hits[i][j] != nullptr) 
            {
                outputfile->WriteObject(fine_time_hits[i][j], Form("fine_time_hits_%i_%i", i, j),"RECREATE");
            }
        }
    }
    c4LOG(info,Form("Written fine time calibrations (i.e. raw fine time histograms) to  %s",fine_time_histo_outfile.Data()));

    outputfile->Close();

}

/*
Read the fine time histograms which are written by the function above.

Assumes the names of the histograms are fine_time_hist_module_channel and that they have 1024 bins (since the TAMEX fine time is written with 2^10 bits).
*/
void bPlastReader::ReadFineTimeHistosFromFile() {

    TFile* inputfile = TFile::Open(fine_time_histo_infile, "READ");
    if (!inputfile || inputfile->IsZombie()) {
        c4LOG(fatal, "File to read histos not opened.");
    }

    c4LOG(info,"Reading the histograms used in fine time calibration from file.");
    fine_time_hits = new TH1I**[NBoards];
    for (int i = 0; i < NBoards; i++) {
        fine_time_hits[i] = new TH1I*[NChannels];
        for (int j = 0; j < NChannels; j++) fine_time_hits[i][j] = nullptr;
    }


    if (!fine_time_hits) {c4LOG(fatal,"fine_time_hits not declared.");}
    for (int i = 0; i < NBoards; i++) {
        if (!fine_time_hits[i]) {c4LOG(fatal,"fine_time_hits[i] not declared.");}
        for (int j = 0; j < NChannels; j++) {
            if (fine_time_hits[i][j] != nullptr) {c4LOG(fatal,"fine_time_hits[i][j] not declared.");}
            TH1I* a = nullptr;
            inputfile->GetObject(Form("fine_time_hits_%i_%i", i, j), a);
            if (a) {
                c4LOG(debug1, Form("Accessing i = %i, j = %i",i,j));
                fine_time_hits[i][j] = (TH1I*)a->Clone();
                c4LOG_IF(fatal,fine_time_hits==nullptr,"Failed reading the file for fine time calibration histograms");
                delete a;
                a = nullptr;
                fine_time_hits[i][j]->SetDirectory(0); // this was the trick to access the histos after closing the inputfile.
            }
        }
    }

    inputfile->Close();
    c4LOG(info, Form("Read fine time calibrations (i.e. raw fine time histograms) from %s", fine_time_histo_infile.Data()));
}

/*
Event read loop. Only called by the FairRootManager.

If no fine time calibrations are set, the fine time histograms are filled until a set number of events have been read. Then the FineTimeCalibration is automatically done and the writing of hits start.
If the fine time calibration is set from the start, hits are written immediately.

Some assumptions:
    - Each time a new channel is hit, it must be preceeded by an epoch word.
    - Each time the epoch word changes it is written again.
    - Time hits are time ordered as they are written by the TAMEX module.

*/
Bool_t bPlastReader::Read() //do fine time here:
{
    c4LOG(debug1, "Event Data");

    if (!fData) return kTRUE;

    if ((fNEvent==fine_time_calibration_after)  & (!fine_time_calibration_set))
    {
        c4LOG(info, "Doing fine time calibration now.");
        DoFineTimeCalibration();
        if (fine_time_calibration_save) WriteFineTimeHistosToFile();
    }

    //whiterabbit timestamp:
    wr_t = (((uint64_t)fData->bplast_ts_t[3]) << 48) + (((uint64_t)fData->bplast_ts_t[2]) << 32) + (((uint64_t)fData->bplast_ts_t[1]) << 16) + (uint64_t)(fData->bplast_ts_t[0]);
    

    for (int it_board_number = 0; it_board_number < NBoards; it_board_number++)
    { //per board:
        
        //for (int i = 0; i<32; i++) last_epoch[i] = 0;

        if (fData->bplast_tamex[it_board_number].event_size == 0) continue; // empty event skip
        
        last_word_read_was_epoch = false;
        last_channel_read = 0;

        uint32_t previous_epoch_word = 0; // last seen epoch word

        bplast_last_lead_hit_struct last_tdc_hit; // initialized and reset - this keeps the information of the last lead that was seen
        last_tdc_hit.hit=false;
        last_tdc_hit.lead_epoch_counter = 0;
        last_tdc_hit.lead_coarse_T = 0;
        last_tdc_hit.lead_fine_T = 0;

        //c4LOG(info,"\n\n\n\n New event:");
        //c4LOG(info,Form("Board_id = %i",it_board_number));
        //c4LOG(info,Form("event size: %i",fData->bplast_tamex[it_board_number].event_size));
        //c4LOG(info,Form("time_epoch =  %i, time_coarse = %i, time_fine = %i, time_edge = %i, time_channel = %i",fData->bplast_tamex[it_board_number].time_epoch,fData->bplast_tamex[it_board_number].time_coarse,fData->bplast_tamex[it_board_number].time_fine,fData->bplast_tamex[it_board_number].time_edge,fData->bplast_tamex[it_board_number].time_channel));
        //c4LOG(info,"Here comes data:");

       

        for (int it_hits = 0; it_hits < fData->bplast_tamex[it_board_number].event_size/4 - 3 ; it_hits++)
        {

            //if (it_hits + 1 == fData->bplast_tamex[it_board_number].event_size/4 - 3) { c4LOG(fatal, "Data ends on a epoch..."); }
            //if (fData->bplast_tamex[it_board_number].time_channelv[it_hits] != 1) continue;

            //c4LOG(info, Form("epoch = %i",fData->bplast_tamex[it_board_number].time_epochv[it_hits]));
            //c4LOG(info, Form("coarse = %i",fData->bplast_tamex[it_board_number].time_coarsev[it_hits]));
            //c4LOG(info, Form("fine = %i",fData->bplast_tamex[it_board_number].time_finev[it_hits]));
            //c4LOG(info, Form("channel = %i",fData->bplast_tamex[it_board_number].time_channelv[it_hits]));
            //c4LOG(info, Form("edge = %i \n",fData->bplast_tamex[it_board_number].time_edgev[it_hits]));

            //continue;

            // now operating under the assumption 16.01.2024:
            /*
            the TAMEX readout happens channel after channel with increasing channel number. Multiple hits within one channel are time ordered.
            e.g.: ch1-hit1, ch1-hit2, ch1-hit3, ch2-hit1, ch2-hit2, ch3-hit1, ....., chN-hitN.
            from M Reese 08.12.23
            */


            // check if arrays are overflowing. These could be placed outside this inner loop, but now it sums up the number of lost events in case they are not equal. This indicates an error in UCESB/MBS i believe if this fails.
            if (fData->bplast_tamex[it_board_number].time_epoch <= it_hits) {fNevents_skipped++; continue;}
            if (fData->bplast_tamex[it_board_number].time_fine <= it_hits) {fNevents_skipped++; continue;}
            if (fData->bplast_tamex[it_board_number].time_coarse <= it_hits) {fNevents_skipped++; continue;}
            //the length of the arrays should be equal. 
            if (!(fData->bplast_tamex[it_board_number].time_coarse == fData->bplast_tamex[it_board_number].time_fine && fData->bplast_tamex[it_board_number].time_coarse == fData->bplast_tamex[it_board_number].time_epoch)) {fNevents_skipped++; continue;}

            //any time you see an epoch - this epoch will apply to the next time data words until a new epoch word is written. If the following data in the buffer is another channel then it must be preceeded with another epoch. 
            if (fData->bplast_tamex[it_board_number].time_epochv[it_hits] != 0)
            {
                previous_epoch_word = fData->bplast_tamex[it_board_number].time_epochv[it_hits] & 0xFFFFFFF;
                //if (it_board_number == 1) c4LOG(info,Form("Found epoch for ch = %i, e = %i",next_channel,fData->bplast_tamex[it_board_number].time_epochv[it_hits] & 0xFFFFFFF));
                fNepochwordsread++;
                last_word_read_was_epoch = true;
                continue;
            }

            //from this point we should have seen an epoch for channel id.

            uint32_t channelid = fData->bplast_tamex[it_board_number].time_channelv[it_hits] & 0x7F; // 1-32
            if (channelid == 0) {continue;} // skip channel 0 for now. This is the trigger information.
            
            //if (it_board_number == 1) c4LOG(info,Form("ch = %i, coarse = %i, edge = %i", channelid, fData->bplast_tamex[it_board_number].time_coarsev[it_hits], fData->bplast_tamex[it_board_number].time_edgev[it_hits]));

            if (channelid != 0 && channelid != last_channel_read && !last_word_read_was_epoch) { fNevents_lacking_epoch[it_board_number][channelid-1]++; continue; } // c4LOG(warning, "Event lacking epoch."); continue;}
        

            if (fData->bplast_tamex[it_board_number].time_finev[it_hits] == 0x3FF) {fNevents_TAMEX_fail[it_board_number][channelid-1]++; continue;} // this happens if TAMEX loses the fine time - skip it

            if (channelid != 0 && channelid != last_channel_read && !last_word_read_was_epoch){fNevents_lacking_epoch[it_board_number][channelid-1]++; c4LOG(warning, Form("Event lacking epoch in channel %d",channelid));} // if the channel has changed but no epoch word was seen in between, channel 0 is always the first one so dont check if that s the case.

            if (!(channelid >= last_channel_read)) {c4LOG(fatal, Form("Data format is inconcistent with assumption: Channels are not read out in increasing order. This channel = %i, last channel = %i",channelid,last_channel_read));}

    
            last_word_read_was_epoch = false;
            last_channel_read = channelid;

            //Fill fine times and skip.
            if (!fine_time_calibration_set) 
            {
                fine_time_hits[it_board_number][channelid-1]->Fill(fData->bplast_tamex[it_board_number].time_finev[it_hits]);
                continue;
            }

            bool is_leading = fData->bplast_tamex[it_board_number].time_edgev[it_hits] & 0x1;
            uint32_t coarse_T = fData->bplast_tamex[it_board_number].time_coarsev[it_hits] & 0x7FF;
            double fine_T = GetFineTime(fData->bplast_tamex[it_board_number].time_finev[it_hits],it_board_number,channelid-1);

            //if(it_board_number == 1) c4LOG(info,fData->bplast_tamex[it_board_number].time_edgev[it_hits]);

            if (is_leading)
            { // rise signal:
                //if (it_board_number == 1) c4LOG(info,Form("Found rise: ch = %i, le = %i, lc = %i, lf = %f", channelid, previous_epoch_word,coarse_T,fine_T));
                
                //count number of double leads
                if (last_tdc_hit.hit) {fNevents_second_lead_seen[it_board_number][channelid-1]++;}
                
                last_tdc_hit.hit = true;
                last_tdc_hit.lead_epoch_counter = previous_epoch_word;
                last_tdc_hit.lead_coarse_T = coarse_T;
                last_tdc_hit.lead_fine_T = fine_T;

                
                fNleads_read[it_board_number][channelid-1]++;
                continue;
            }
            else if (!is_leading && last_tdc_hit.hit)
            { 
                //trail and rise are matched
                //if (it_board_number == 1) c4LOG(info,Form("Writing: ch = %i, le = %i lc = %i, lf = %f, te = %i tc = %i, tf = %f ",channelid,last_hits[it_board_number][channelid-1].lead_epoch_counter, last_hits[it_board_number][channelid-1].lead_coarse_T, last_hits[it_board_number][channelid-1].lead_fine_T,last_epoch[channelid-1],coarse_T,fine_T));

                new ((*fArray)[fArray->GetEntriesFast()]) bPlastTwinpeaksData(
                    it_board_number,
                    channelid,
                    //last_hits[it_board_number][channelid-1].lead_epoch_counter,
                    last_tdc_hit.lead_epoch_counter,
                    //last_hits[it_board_number][channelid-1].lead_coarse_T,
                    last_tdc_hit.lead_coarse_T,
                    //last_hits[it_board_number][channelid-1].lead_fine_T,
                    last_tdc_hit.lead_fine_T,

                    previous_epoch_word,
                    coarse_T,
                    fine_T,
                    fData->bplast_ts_subsystem_id,
                    wr_t);
            
                    //reset:

                    last_tdc_hit.hit=false;
                    last_tdc_hit.lead_epoch_counter = 0;
                    last_tdc_hit.lead_coarse_T = 0;
                    last_tdc_hit.lead_fine_T = 0;
                
                fNmatched[it_board_number][channelid-1]++;
                fNtrails_read[it_board_number][channelid-1]++;

                continue;
            }
            else
            {
                // do nothing, trail found with no rise.
                fNevents_trail_seen_no_lead[it_board_number][channelid-1]++;
                fNtrails_read[it_board_number][channelid-1]++;
            }
        }
    } // boards
    
    fNEvent += 1;
    return kTRUE;
}


#define ANSI_COLOR_RESET "\033[0m"
#define ANSI_COLOR_BLUE "\033[44m"
#define ANSI_COLOR_RED     "\033[41m"
#define ANSI_COLOR_GREEN   "\033[42m"


/*
Playing with colors :D
*/
void bPlastReader::PrintStatistics(){
    std::ostringstream oss;
    // Print column labels
    oss << "\n";
    oss << ANSI_COLOR_GREEN << "bPLAST TAMEX4 STATISTICS - boards from 0-N. TDC channels are numbered from 1-32 fast(odd)-slow(even)" << ANSI_COLOR_RESET;
    oss << "\n";
    oss << "Number of epoch words read in total: " << fNepochwordsread << "\n";
    oss << "Number of events skipped due to c4/mbs inconsitency of array size: " << fNevents_skipped << "\n";
    oss << ANSI_COLOR_RED << "Number of trails read per channel." << ANSI_COLOR_RESET;
    oss << std::setw(12)  << ANSI_COLOR_BLUE << "\n      Channel "; // Empty space for the top-left cell
    for (int j = 0; j < NChannels; ++j) {
        oss << std::setw(12) << j + 1; // Assuming channels are 1-indexed
    }
    oss << ANSI_COLOR_RESET <<  '\n';

    // Print array values with row labels
    for (int i = 0; i < NBoards; ++i) {
        oss << std::setw(12) << ANSI_COLOR_BLUE << "Board " << i << ANSI_COLOR_RESET; // Assuming boards are 1-indexed
        for (int j = 0; j < NChannels; ++j) {
            oss << std::setw(12) << fNtrails_read[i][j];
        }
        oss << '\n';
    }
    oss << "\n\n";

    oss << ANSI_COLOR_RED << "Number of leads read per channel." << ANSI_COLOR_RESET;
    oss << std::setw(12)  << ANSI_COLOR_BLUE << "\n      Channel "; // Empty space for the top-left cell
    for (int j = 0; j < NChannels; ++j) {
        oss << std::setw(12) << j + 1; // Assuming channels are 1-indexed
    }
    oss << ANSI_COLOR_RESET <<  '\n';

    // Print array values with row labels
    for (int i = 0; i < NBoards; ++i) {
        oss << std::setw(12) << ANSI_COLOR_BLUE << "Board " << i << ANSI_COLOR_RESET; // Assuming boards are 1-indexed
        for (int j = 0; j < NChannels; ++j) {
            oss << std::setw(12) << fNleads_read[i][j];
        }
        oss << '\n';
    }
    oss << "\n";

    oss << ANSI_COLOR_RED << "Number of matched lead-trail read per channel." << ANSI_COLOR_RESET;
    oss << std::setw(12)  << ANSI_COLOR_BLUE << "\n      Channel "; // Empty space for the top-left cell
    for (int j = 0; j < NChannels; ++j) {
        oss << std::setw(12) << j + 1; // Assuming channels are 1-indexed
    }
    oss << ANSI_COLOR_RESET <<  '\n';

    // Print array values with row labels
    for (int i = 0; i < NBoards; ++i) {
        oss << std::setw(12) << ANSI_COLOR_BLUE << "Board " << i << ANSI_COLOR_RESET; // Assuming boards are 1-indexed
        for (int j = 0; j < NChannels; ++j) {
            oss << std::setw(12) << fNmatched[i][j];
        }
        oss << '\n';
    }
    oss << "\n";

    oss << ANSI_COLOR_RED << "Number of events lacking epoch." << ANSI_COLOR_RESET;
    oss << std::setw(12)  << ANSI_COLOR_BLUE << "\n      Channel "; // Empty space for the top-left cell
    for (int j = 0; j < NChannels; ++j) {
        oss << std::setw(12) << j + 1; // Assuming channels are 1-indexed
    }
    oss << ANSI_COLOR_RESET <<  '\n';

    // Print array values with row labels
    for (int i = 0; i < NBoards; ++i) {
        oss << std::setw(12) << ANSI_COLOR_BLUE << "Board " << i << ANSI_COLOR_RESET; // Assuming boards are 1-indexed
        for (int j = 0; j < NChannels; ++j) {
            oss << std::setw(12) << fNevents_lacking_epoch[i][j];
        }
        oss << '\n';
    }
    oss << "\n";


    oss << ANSI_COLOR_RED << "Number of TAMEX fails where the error 0x3FF is set." << ANSI_COLOR_RESET;
    oss << std::setw(12)  << ANSI_COLOR_BLUE << "\n      Channel "; // Empty space for the top-left cell
    for (int j = 0; j < NChannels; ++j) {
        oss << std::setw(12) << j + 1; // Assuming channels are 1-indexed
    }
    oss << ANSI_COLOR_RESET <<  '\n';

    // Print array values with row labels
    for (int i = 0; i < NBoards; ++i) {
        oss << std::setw(12) << ANSI_COLOR_BLUE << "Board " << i << ANSI_COLOR_RESET; // Assuming boards are 1-indexed
        for (int j = 0; j < NChannels; ++j) {
            oss << std::setw(12) << fNevents_TAMEX_fail[i][j];
        }
        oss << '\n';
    }
    oss << "\n";


    oss << ANSI_COLOR_RED << "Number of times a second lead is seen (i.e. a lead-lead in the channel) keeping only the last lead.." << ANSI_COLOR_RESET;
    oss << std::setw(12)  << ANSI_COLOR_BLUE << "\n      Channel "; // Empty space for the top-left cell
    for (int j = 0; j < NChannels; ++j) {
        oss << std::setw(12) << j + 1; // Assuming channels are 1-indexed
    }
    oss << ANSI_COLOR_RESET <<  '\n';

    // Print array values with row labels
    for (int i = 0; i < NBoards; ++i) {
        oss << std::setw(12) << ANSI_COLOR_BLUE << "Board " << i << ANSI_COLOR_RESET; // Assuming boards are 1-indexed
        for (int j = 0; j < NChannels; ++j) {
            oss << std::setw(12) << fNevents_second_lead_seen[i][j];
        }
        oss << '\n';
    }
    oss << "\n";

    oss << ANSI_COLOR_RED << "Nnumber of times a trail is seen without a preceeding lead - skipping this event." << ANSI_COLOR_RESET;
    oss << std::setw(12)  << ANSI_COLOR_BLUE << "\n      Channel "; // Empty space for the top-left cell
    for (int j = 0; j < NChannels; ++j) {
        oss << std::setw(12) << j + 1; // Assuming channels are 1-indexed
    }
    oss << ANSI_COLOR_RESET <<  '\n';

    // Print array values with row labels
    for (int i = 0; i < NBoards; ++i) {
        oss << std::setw(12) << ANSI_COLOR_BLUE << "Board " << i << ANSI_COLOR_RESET; // Assuming boards are 1-indexed
        for (int j = 0; j < NChannels; ++j) {
            oss << std::setw(12) << fNevents_trail_seen_no_lead[i][j];
        }
        oss << '\n';
    }
    oss << "\n";

    
    // Get the formatted string
    TString formattedString = oss.str();
    
    c4LOG(info, formattedString);
}

/*
Memory management. Do not touch.
*/
void bPlastReader::Reset()
{
    // reset output array
    fArray->Clear();
       
}

ClassImp(bPlastReader)
