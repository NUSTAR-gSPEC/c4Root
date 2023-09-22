#ifndef AidaReader_h
#define AidaReader_h

#include "c4Reader.h"

#include <Rtypes.h>

extern "C"
{
    #include "ext_h101_aida.h" // name??
}

class AidaReader : public c4Reader
{
    public:
        AidaReader(EXT_STR_h101_aida_onion*, size_t);

        virtual ~AidaReader();

        virtual Bool_t Init(ext_data_struct_info*) override;

        virtual Bool_t Read() override;

        virtual void Reset() override;

        void SetOnline(Bool_t option) { fOnline = option; }
    
    private:
        unsigned int fNEvent;

        EXT_STR_h101_aida_onion* fData;

        size_t fOffset;

        Bool_t fOnline;

        //TClonesArray* fArray;
        // Data to register here
    
    public:
        ClassDefOverride(AidaReader, 0);
};

#endif
