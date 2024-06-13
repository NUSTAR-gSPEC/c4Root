// -*- C++ -*-

#include "../../common/whiterabbit.spec"
#include "../../common/gsi_tamex4.spec"

SUBEVENT(gplast_subev)
{
    select optional
    {
         ts = TIMESTAMP_WHITERABBIT_EXTENDED(id=0x1600);       
    };
    
    select optional
    {
        trigger_window = TAMEX4_HEADER();
    };
 
    select several 
    {
        padding = TAMEX4_PADDING();
    }
    select several
    {
        tamex[0] = TAMEX4_SFP(sfp=0,card=0);
	tamex[1] = TAMEX4_SFP(sfp=0,card=1);
	tamex[2] = TAMEX4_SFP(sfp=0,card=2);
        tamex[3] = TAMEX4_SFP(sfp=0,card=3);
    }  
}


EVENT
{   
    gplast = gplast_subev(type = 10, subtype = 1, procid = 80, control = 20);

    //ignore_unknown_subevent;
};

