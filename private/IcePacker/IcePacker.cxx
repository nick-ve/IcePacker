/**
 * class: IcePacker
 *
 * Version $Id: $
 *
 *
 * (c) 2009 IceCube Collaboration
 * @file IcePacker.cxx
 * @date $Date: $
 * @author duvoort
 **/
//Modified 3 Dec 2014: Add NcDevice for the common variables
//Add a signal slot for HiveSplit pulses
//Add characteristics to several tracks
//L. Brayeur IIHE-VUB

#include "../IcePackBase/IPBHeaders.h"
#include "IcePacker/IcePacker.h"
#include "dataclasses/calibration/I3Calibration.h"
#include "dataclasses/geometry/I3Geometry.h"
#include "dataclasses/I3Position.h"
#include "icetray/I3Units.h"
#include "dataclasses/physics/I3Trigger.h"
#include "dataclasses/physics/I3Particle.h"
#include "dataclasses/physics/I3MCTreeUtils.h"
#include "dataclasses/physics/I3MCTreePhysicsLibrary.hh"
//@@@ #include "dataclasses/I3TreeUtils.h"
#include "dataclasses/I3Tree.h"
#include "dataclasses/physics/I3TriggerHierarchy.h"
#include "dataclasses/TriggerKey.h"
#include "dataclasses/I3Double.h"
#include "dataclasses/I3Map.h"
#include "dataclasses/physics/I3EventHeader.h"
#include "dataclasses/physics/I3DOMLaunch.h"
#include "dataclasses/physics/I3Waveform.h"
#include "dataclasses/physics/I3RecoPulse.h"
#include "paraboloid/I3ParaboloidFitParams.h"
#include "gulliver/I3LogLikelihoodFitParams.h"
//@@@ #include "linefit/I3LineFitParams.h"
//@@@#include "linefit/I3LineFit.h"
#include "recclasses/I3LineFitParams.h"
//LB: Add necessary headers
//@@@ #include "cramer-rao/CramerRaoParams.h"
//@@@#include "cramer-rao/CramerRao.h"
#include "recclasses/CramerRaoParams.h"
//@@@@#include "finiteReco/I3FiniteCuts.h"
#include "recclasses/I3FiniteCuts.h"
/***
#include "CommonVariables/direct_hits/I3DirectHitsValues.h"
#include "CommonVariables/hit_multiplicity/I3HitMultiplicityValues.h"
#include "CommonVariables/hit_statistics/I3HitStatisticsValues.h"
#include "CommonVariables/track_characteristics/I3TrackCharacteristicsValues.h"
***/
/***
#include "CommonVariables/direct_hits/calculator.h"
#include "CommonVariables/hit_multiplicity/calculator.h"
#include "CommonVariables/hit_statistics/calculator.h"
#include "CommonVariables/track_characteristics/calculator.h"
***/
#include "recclasses/I3DirectHitsValues.h"
#include "recclasses/I3HitMultiplicityValues.h"
#include "recclasses/I3HitStatisticsValues.h"
#include "recclasses/I3TrackCharacteristicsValues.h"
//LB:------------------------

// Next statement to be used for older IceTray versions
// The I3FilterResult.h has been moved into "dataclasses/physics" afterwards
// #include "jebclasses/I3FilterResult.h"
#include "dataclasses/physics/I3FilterResult.h"
#include <boost/foreach.hpp>

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <TTree.h>
#include <TFile.h>

#include <TSystem.h>
#include <TROOT.h>

I3_MODULE(IcePacker);

using namespace std;

IcePacker::IcePacker(const I3Context& ctx)
: I3Module(ctx)
{
    // adding an outbox
    rootfilename_ = "IcePackdata.root";
    treename_ = "T";
    domlMapInIceName_ = "InIceRawData";
    domlMapIceTopName_ = "IceTopRawData";
    wvfrmMapInIceATWDName_ = "InIceCalibratedATWD";
    wvfrmMapInIceFADCName_ = "InIceCalibratedFADC";
    wvfrmMapIceTopATWDName_ = "IceTopCalibratedATWD";
    wvfrmMapIceTopFADCName_ = "IceTopCalibratedFADC";
    InIcePulsesName_ = "InIcePulseSeries";
    IceTopPulsesName_ = "IceTopPulseSeries";
    // No need to initialize InIcePulsesNames_
    // No need to initialize IceTopPulsesNames_
    // No need to initialize fSubeventSels
    // No need to initialize fFilterSels
    fFilterPrescale=false;
    fSelNamePattern=false;
    mctreeName_ = "I3MCTree";
    mcweightinfo_ = "I3MCWeightDict";
    trackextractOn_ = true;
    eventheaderName_ = "I3EventHeader";
    triggerName_ = "I3TriggerHierarchy";
    filterName_ = "FilterMask";
    maxtreesize = 200000000;
    multimctrack_=false;
    
    AddParameter("RootFileName", "IceEvent output file name", rootfilename_);
    AddParameter("TreeName", "root tree name", treename_);
    AddParameter("TrackExtract","Extracts the track information",trackextractOn_);
    AddParameter("TriggerName","Name of trigger hierarchy",triggerName_);
    AddParameter("FilterName","Name of filter data container",filterName_);
    AddParameter("EventHeaderName","Name of I3EventHeader",eventheaderName_);
    AddParameter("DOMmapInIceName","Name of InIce DOM map to be used",domlMapInIceName_);
    AddParameter("DOMmapIceTopName","Name of IceTop DOM map to be used",domlMapIceTopName_);
    AddParameter("WaveformmapInIceATWDName","Name of InIce ATWD Waveform map to be used",wvfrmMapInIceATWDName_);
    AddParameter("WaveformmapInIceFADCName","Name of InIce FADC Waveform map to be used",wvfrmMapInIceFADCName_);
    AddParameter("WaveformmapIceTopATWDName","Name of IceTop ATWD Waveform map to be used",wvfrmMapIceTopATWDName_);
    AddParameter("WaveformmapIceTopFADCName","Name of IceTop FADC Waveform map to be used",wvfrmMapIceTopFADCName_);
    AddParameter("InIcePulseSeriesName","Name of InIce RecoPulseSeries to be used",InIcePulsesName_);
    AddParameter("IceTopPulseSeriesName","Name of IceTop RecoPulseSeries to be used",IceTopPulsesName_);
    AddParameter("InIcePulseSeriesNames","Name(s) of InIce RecoPulseSeries to be used in a vector of strings",InIcePulsesNames_);
    AddParameter("IceTopPulseSeriesNames","Name(s) of IceTop RecoPulseSeries to be used in a vector of strings",IceTopPulsesNames_);
    AddParameter("SubeventSelections","Name(s) of the subevent stream(s) to be extracted in a vector of strings",fSubeventSels);
    AddParameter("FilterSelections","Name(s) of the filter(s) to be extracted in a vector of strings",fFilterSels);
    AddParameter("FilterPrescale","Set to true to require prescale flag to be set for selected filters",fFilterPrescale);
    AddParameter("SelNamePattern","Set to true to only match the name pattern for the selection names",fSelNamePattern);
    AddParameter("MCTreeName","Name of MCTree to get the MC particle from",mctreeName_);
    AddParameter("MCWeightInfo","Info of the MC production from the MCWeightDictionary",mcweightinfo_);
    AddParameter("MultiMCtrack","Set to true for conversion of all MC track information",multimctrack_);
    AddParameter("setmaxtreesize","Setting a maximum treesize",maxtreesize);
    
    AddOutBox("OutBox");
}

void IcePacker::Configure()
{
    cout<<"Configuring IcePacker"<<endl;
    GetParameter("RootFileName", rootfilename_);
    GetParameter("TreeName", treename_);
    GetParameter("TrackExtract",trackextractOn_);
    GetParameter("TriggerName",triggerName_);
    GetParameter("FilterName",filterName_);
    GetParameter("EventHeaderName",eventheaderName_);
    GetParameter("DOMmapInIceName",domlMapInIceName_);
    GetParameter("DOMmapIceTopName",domlMapIceTopName_);
    GetParameter("WaveformmapInIceATWDName",wvfrmMapInIceATWDName_);
    GetParameter("WaveformmapInIceFADCName",wvfrmMapInIceFADCName_);
    GetParameter("WaveformmapIceTopATWDName",wvfrmMapIceTopATWDName_);
    GetParameter("WaveformmapIceTopFADCName",wvfrmMapIceTopFADCName_);
    GetParameter("InIcePulseSeriesName",InIcePulsesName_);
    GetParameter("IceTopPulseSeriesName",IceTopPulsesName_);
    GetParameter("InIcePulseSeriesNames",InIcePulsesNames_);
    GetParameter("IceTopPulseSeriesNames",IceTopPulsesNames_);
    GetParameter("SubeventSelections",fSubeventSels);
    GetParameter("FilterSelections",fFilterSels);
    GetParameter("FilterPrescale",fFilterPrescale);
    GetParameter("SelNamePattern",fSelNamePattern);
    GetParameter("MCTreeName",mctreeName_);
    GetParameter("MCWeightInfo",mcweightinfo_);
    GetParameter("MultiMCtrack",multimctrack_);
    GetParameter("setmaxtreesize",maxtreesize);
    
    //Creating outputfile
    tfile_ = new TFile(rootfilename_.c_str(), "RECREATE");
    ttree_ = new TTree(treename_.c_str(), "T", 10000000);
    ttree_->SetMaxTreeSize(maxtreesize);
    
    evt= new IceEvent();
    evt->SetOwner();
    
    // Branch in the tree for the event structure
    if (ttree_) ttree_->Branch("IceEvent","IceEvent",&evt,32000,0);
    
    count_ = 0;
    log_info("Finished Configuring IcePacker");
}

void IcePacker::Reconfigure()
{
    log_fatal("IcePacker: Cannot reconfigure this module");
}

void IcePacker::Physics(I3FramePtr frame)
{
    count_++;
    log_info("IcePacker Physics: Event number: %i",count_);
    
    // Retrieving general event info from the header
    log_debug("ABOUT TO READ THE HEADER");
    I3EventHeaderConstPtr header = frame->Get<I3EventHeaderConstPtr>(eventheaderName_);
    Int_t eventtimemjd = 0;
    Int_t eventtimemjds = 0;
    Int_t eventtimemjdns = 0;
    Int_t eventid = 0;
    Int_t runid = 0;
    TString subeventstream="none";
    Int_t subeventid=0;
    if (header)
    {
        I3Time starttimeevt = header->GetStartTime();
        eventtimemjd = (Int_t)starttimeevt.GetModJulianDay();
        eventtimemjds = (Int_t)starttimeevt.GetModJulianSec();
        eventtimemjdns = (Int_t)starttimeevt.GetModJulianNanoSec();
        eventid = (Int_t)header->GetEventID();
        runid = (Int_t)header->GetRunID();
        subeventstream=header->GetSubEventStream();
        subeventid=(Int_t)header->GetSubEventID();
    }
    
    if(eventid==0) eventid = count_; //Sometimes a MC file has bad event indexing
    
    // Skip subevents that are not in the selection list
    Int_t iskip=0;
    TString subeventname;
    Int_t nsubs=fSubeventSels.size();
    if (nsubs) iskip=1;
    for (Int_t jsub=0; jsub<nsubs; jsub++)
    {
        subeventname=fSubeventSels[jsub];
        if (subeventstream==subeventname || subeventname=="all") iskip=0;
        if (fSelNamePattern && subeventstream.Contains(subeventname.Data())) iskip=0;
    }
    
    if (iskip) // Event was not one of the selected subevent streams
    {
        log_debug("IcePacker: Finished event number : %i",count_);
        PushFrame(frame,"OutBox");
        return;
    }
    
    // The event filter mask info
    // Skip events that do not satisfy one of the selected filters
    TString filtername="none";
    TString filtersel;
    Int_t nfilt=fFilterSels.size();
    if (nfilt) iskip=1;
    NcDevice filter;
    filter.SetNameTitle("Filter","IceCube event filter mask");
    if (frame->Has(filterName_))
    {
        I3FilterResultMapConstPtr filterResultMap=frame->Get<I3FilterResultMapConstPtr>(filterName_);
        NcSignal sx;
        Float_t condition=0;
        Float_t prescale=0;
        BOOST_FOREACH(I3FilterResultMap::const_reference filterResult,*filterResultMap)
        {
            filtername=filterResult.first;
            condition=0;
            if(filterResult.second.conditionPassed) condition=1;
            prescale=0;
            if(filterResult.second.prescalePassed) prescale=1;
            sx.Reset();
            sx.SetName(filtername);
            sx.AddNamedSlot("condition");
            sx.SetSignal(condition,"condition");
            sx.AddNamedSlot("prescale");
            sx.SetSignal(prescale,"prescale");
            
            filter.AddHit(sx);
            
            // Check whether this event satisfies one of the selected filters
            for (Int_t jfilt=0; jfilt<nfilt; jfilt++)
            {
                filtersel=fFilterSels[jfilt];
                if (condition<0.5) continue;
                if (fFilterPrescale && prescale<0.5) continue;
                if (filtername==filtersel || filtersel=="all") iskip=0;
                if (fSelNamePattern && filtername.Contains(filtersel.Data())) iskip=0;
            }
        }
    }
    
    if (iskip) // Event did not satisfy one of the selected filters
    {
        log_debug("IcePacker: Finished event number : %i",count_);
        PushFrame(frame,"OutBox");
        return;
    }
    
    // Initialise some standard devices
    daq.SetNameTitle("DAQ","DAQ system used");
    mcinfo.SetNameTitle("MCInfo","The MC production info");
    trig.SetNameTitle("Trigger","IceCube event trigger(s)");
    
    NcDevice splitter;
    splitter.SetNameTitle("Splitter","IceCube subevent name and ID");
    
    evt->Reset();
    
    evt->SetRunNumber(runid);
    evt->SetEventNumber(eventid);
    evt->SetMJD(eventtimemjd,eventtimemjds,eventtimemjdns);

    // Correct for faulty leap second insertion.
    // For the run periods indicated below, the IceCube clock was 1 sec. behind.
    // So, here we add 1 sec. to the IceCube time for those run periods.
    if (runid>=120398 && runid<=126377) evt->Add(0,1,0,0);
    
    // Entering the DAQ information
    daq.Reset(1);
    daq.AddNamedSlot("JEB");
    daq.SetSignal(1,"JEB");
    evt->AddDevice(daq);
    
    // Entering the subevent stream information
    splitter.AddNamedSlot(subeventstream);
    splitter.SetSignal(subeventid,subeventstream);
    evt->AddDevice(splitter);
    
    // Entering the filter information
    evt->AddDevice(filter);
    
    // Simulation data contains weighting information and production info
    // Extracting this from the MCWeightDict
    I3Double parWeight(0.);
    mcinfo.Reset(1);
    I3MapStringDoubleConstPtr mcweightmap = frame ->Get<I3MapStringDoubleConstPtr>(mcweightinfo_);
    if (mcweightmap)
    {
        Double_t mcweightitemvalue;
        
        log_debug("It found the MCInfo map");
        for (I3MapStringDouble::const_iterator mapiter = mcweightmap->begin(); mapiter != mcweightmap->end(); mapiter++)
        {
            TString mcinfoname = mapiter->first;
            mcinfo.AddNamedSlot(mcinfoname);
            mcweightitemvalue = mapiter->second;
            mcinfo.SetSignal(mcweightitemvalue,mcinfoname.Data());
            if (mcinfoname.Index("OneWeight")>=0) evt->SetWeight(mcweightitemvalue);
        }
        
    }
    
    //LB: Add weight for Atm nugen and Corsika
    else if (frame->Has("CorsikaWeightMap")){
        I3MapStringDoubleConstPtr corsweightmap = frame ->Get<I3MapStringDoubleConstPtr>("CorsikaWeightMap");
        if(corsweightmap){
            
            Double_t corsweightitemvalue;
            
            for (I3MapStringDouble::const_iterator mapiter = corsweightmap->begin(); mapiter != corsweightmap->end(); mapiter++)
            {
                TString corsinfoname = mapiter->first;
                mcinfo.AddNamedSlot(corsinfoname);
                corsweightitemvalue = mapiter->second;
                mcinfo.SetSignal(corsweightitemvalue,corsinfoname.Data());
            }
        }
        
    }
    if(frame->Has("AtmWeight")){
        parWeight = frame->Get<I3Double>((string)"AtmWeight");
        mcinfo.AddNamedSlot("AtmWeight");
        mcinfo.SetSignal(parWeight.value,"AtmWeight");
    }
    if(frame->Has("GaisserH3aWeight")){
        parWeight = frame->Get<I3Double>((string)"GaisserH3aWeight");
        mcinfo.AddNamedSlot("GaisserH3aWeight");
        mcinfo.SetSignal(parWeight.value,"GaisserH3aWeight");
    }
    if(frame->Has("HoerandelWeight")){
        parWeight = frame->Get<I3Double>((string)"HoerandelWeight");
        mcinfo.AddNamedSlot("HoerandelWeight");
        mcinfo.SetSignal(parWeight.value,"HoerandelWeight");
    }
    
    if(mcinfo.GetNslots()!=0) evt->AddDevice(mcinfo);
    
    
    
    //LB: Add a device with the Common variables
    NcDevice comVar;
    comVar.SetNameTitle("CommonVariables","Results of the Common Variables for the event");
    Int_t addcom=0; //NvE: To flag presence of Common variables
    if (frame->Has("HitStatisticsValues_HiveSplit"))
    {
     addcom=1;
     I3HitStatisticsValues statisticsValues=frame->Get<I3HitStatisticsValues>("HitStatisticsValues_HiveSplit");
     //NvE: Suffix _HS added to identify HiveSplit produced values and COG X,Y,Z coordinates added
     comVar.AddNamedSlot("COGX_HS");
     comVar.SetSignal(statisticsValues.GetCOG().GetX(),"COGX_HS");
     comVar.AddNamedSlot("COGY_HS");
     comVar.SetSignal(statisticsValues.GetCOG().GetY(),"COGY_HS");
     comVar.AddNamedSlot("COGZ_HS");
     comVar.SetSignal(statisticsValues.GetCOG().GetZ(),"COGZ_HS");
     comVar.AddNamedSlot("COGZSigma_HS");
     comVar.SetSignal(statisticsValues.GetCOGZSigma(),"COGZSigma_HS");
     comVar.AddNamedSlot("MinPulseTime_HS");
     comVar.SetSignal(statisticsValues.GetMinPulseTime(),"MinPulseTime_HS");
     comVar.AddNamedSlot("MaxPulseTime_HS");
     comVar.SetSignal(statisticsValues.GetMaxPulseTime(),"MaxPulseTime_HS");
     comVar.AddNamedSlot("QMaxDoms_HS");
     comVar.SetSignal(statisticsValues.GetQMaxDoms(),"QMaxDoms_HS");
     comVar.AddNamedSlot("QTotPulses_HS");
     comVar.SetSignal(statisticsValues.GetQTotPulses(),"QTotPulses_HS");
     comVar.AddNamedSlot("ZMin_HS");
     comVar.SetSignal(statisticsValues.GetZMin(),"ZMin_HS");
     comVar.AddNamedSlot("ZMax_HS");
     comVar.SetSignal(statisticsValues.GetZMax(),"ZMax_HS");
     comVar.AddNamedSlot("ZMean_HS");
     comVar.SetSignal(statisticsValues.GetZMean(),"ZMean_HS");
     comVar.AddNamedSlot("ZSigma_HS");
     comVar.SetSignal(statisticsValues.GetZSigma(),"ZSigma_HS");
     comVar.AddNamedSlot("ZTravel_HS");
     comVar.SetSignal(statisticsValues.GetZTravel(),"ZTravel_HS");
        
//NvE: Commented out, since it is invoked later        evt->AddDevice(comVar);
        
    }
    //LB:--------------------------

    //NvE: Extension of the Common variables with the CVMultiplicity from the pure L2 data
    if (frame->Has("CVMultiplicity"))
    {
     addcom=1;
     I3HitMultiplicityValues multiplicityValues=frame->Get<I3HitMultiplicityValues>("CVMultiplicity");
     comVar.AddNamedSlot("NHitStrings");
     comVar.SetSignal(multiplicityValues.GetNHitStrings(),"NHitStrings");
     comVar.AddNamedSlot("NHitDoms");
     comVar.SetSignal(multiplicityValues.GetNHitDoms(),"NHitDoms");
     comVar.AddNamedSlot("NHitDomsOnePulse");
     comVar.SetSignal(multiplicityValues.GetNHitDomsOnePulse(),"NHitDomsOnePulse");
     comVar.AddNamedSlot("NPulses");
     comVar.SetSignal(multiplicityValues.GetNPulses(),"NPulses");
    }

    //NvE: Extension of the Common variables with the CVStatistics from the pure L2 data
    if (frame->Has("CVStatistics"))
    {
     addcom=1;
     I3HitStatisticsValues statisticsValues=frame->Get<I3HitStatisticsValues>("CVStatistics");
     comVar.AddNamedSlot("COGX");
     comVar.SetSignal(statisticsValues.GetCOG().GetX(),"COGX");
     comVar.AddNamedSlot("COGY");
     comVar.SetSignal(statisticsValues.GetCOG().GetY(),"COGY");
     comVar.AddNamedSlot("COGZ");
     comVar.SetSignal(statisticsValues.GetCOG().GetZ(),"COGZ");
     comVar.AddNamedSlot("COGZSigma");
     comVar.SetSignal(statisticsValues.GetCOGZSigma(),"COGZSigma");
     comVar.AddNamedSlot("MinPulseTime");
     comVar.SetSignal(statisticsValues.GetMinPulseTime(),"MinPulseTime");
     comVar.AddNamedSlot("MaxPulseTime");
     comVar.SetSignal(statisticsValues.GetMaxPulseTime(),"MaxPulseTime");
     comVar.AddNamedSlot("QMaxDoms");
     comVar.SetSignal(statisticsValues.GetQMaxDoms(),"QMaxDoms");
     comVar.AddNamedSlot("QTotPulses");
     comVar.SetSignal(statisticsValues.GetQTotPulses(),"QTotPulses");
     comVar.AddNamedSlot("ZMin");
     comVar.SetSignal(statisticsValues.GetZMin(),"ZMin");
     comVar.AddNamedSlot("ZMax");
     comVar.SetSignal(statisticsValues.GetZMax(),"ZMax");
     comVar.AddNamedSlot("ZMean");
     comVar.SetSignal(statisticsValues.GetZMean(),"ZMean");
     comVar.AddNamedSlot("ZSigma");
     comVar.SetSignal(statisticsValues.GetZSigma(),"ZSigma");
     comVar.AddNamedSlot("ZTravel");
     comVar.SetSignal(statisticsValues.GetZTravel(),"ZTravel");
    }

    if (addcom) evt->AddDevice(comVar); // NvE
    
    
    // Retrieve the trigger info
    ReadTrigger(frame);
    
    // Retireve the waveform information
    WaveformExtractor(frame);
    
    // Retrieve the DOM Launch information
    DOMLaunchExtractor(frame);
    
    // Retrieve the hit information
    RecoPulseExtractor(frame);
    
    // Retrieve reconstructed and/or simulation tracks
    if(trackextractOn_) TrackExtractor(frame);
    
    // Write the event to the Icepack output Tree
    ttree_->Fill();
    
    log_debug("IcePacker: Finished event number : %i",count_);
    PushFrame(frame,"OutBox");
}

//____________________________________________________________________________
//Read trigger information
//____________________________________________________________________________
void IcePacker::ReadTrigger(I3FramePtr frame){
    log_debug("ABOUT TO DO TRIGGER");
    
    NcSignal s;
    Double_t trigtime;
    Double_t triglength;
    TString trigname;
    const char* trigsource;
    const char* trigtype;
    Int_t configid=0;
    
    trig.Reset( 1 ) ;
    
    I3TriggerHierarchyConstPtr triggerh = frame->Get<I3TriggerHierarchyConstPtr>(triggerName_);
    I3TriggerHierarchy::iterator iter;
    TriggerKey thiskey;
    
    if(triggerh){
        log_debug("Found the trigger hierarchy.");
        for(iter = triggerh->begin(); iter != triggerh->end(); iter++){
            trigtime = iter->GetTriggerTime();
            triglength = iter->GetTriggerLength();
            thiskey = iter->GetTriggerKey();
            trigsource = thiskey.GetSourceString();
            trigtype = thiskey.GetTypeString();
            configid=0;
            if (thiskey.CheckConfigID()) configid=thiskey.GetConfigID();
            
            trigname=trigsource;
            trigname+="/";
            trigname+=trigtype;
            trigname+="/ID";
            trigname+=configid;
            
            log_debug("Found trigger : %s at time: %d",trigname.Data(),trigtime);
            
            s.Reset(1);
            s.SetName(trigname);
            s.SetSlotName("trig_pulse_le",1);
            s.SetSignal(trigtime,1);
            s.SetSlotName("trig_pulse_tot",2);
            s.SetSignal(triglength,2);
            trig.AddHit(s);
        }
        evt->AddDevice(trig);
    }
}


//Perform the actual waveform extraction
//_________________________________________________________________________________________
void IcePacker::WaveformExtractor(I3FramePtr frame)
{
    log_debug("ABOUT TO DO WAVEFORMEXTRACTION");
    
    Int_t OM;
    Int_t String;
    Int_t omid;
    TString hname;
    TH1F histo;
    vector<double> waveform;
    Double_t starttimewaveform;
    Int_t numberbinswaveform;
    Double_t binwidth;
    Int_t ibin;
    Nc3Vector omr;
    Double_t ompos[3];
    
    // Reading the I3Waveforms from the frame Starting with the In-Ice part
    // First the InIce ATWD Waveforms
    I3WaveformSeriesMapConstPtr wvfrmmapInIceATWD = frame->Get<I3WaveformSeriesMapConstPtr>(wvfrmMapInIceATWDName_);
    if(wvfrmmapInIceATWD){
        // loop over the OMs
        for(I3WaveformSeriesMap::const_iterator map_iter = wvfrmmapInIceATWD->begin(); map_iter != wvfrmmapInIceATWD->end(); map_iter++){
            if(map_iter->second.size() == 0) continue;
            
            String = map_iter->first.GetString();
            OM = map_iter->first.GetOM();
            omid = icdom.GetOMId(String,OM);
            log_debug("ABOUT TO EXTRACT OM INFO of omid : %i",omid);
            
            //Initializing the OM device in case of absense
            om = (IceGOM*)evt->GetIdDevice(omid);
            if(!om){
                if(String>=79 && String<=90){
                    idom=&dcdom;
                }else if(String<=78){
                    idom=&icdom;
                }
                idom->Reset(1);
                idom->SetUniqueID(omid);
                
                const I3Geometry& geometry = frame->Get<I3Geometry>();
                const I3OMGeoMap& om_geo = geometry.omgeo;
                I3OMGeoMap::const_iterator geo_iter = om_geo.find(map_iter->first);
                assert(geo_iter != om_geo.end());
                const I3Position om_pos = (geo_iter->second).position;
                ompos[0] = (Double_t)om_pos.GetX();
                ompos[1] = (Double_t)om_pos.GetY();
                ompos[2] = (Double_t)om_pos.GetZ();
                omr.SetVector(ompos,"car");
                idom->SetPosition(omr);
                
                evt->AddDevice(idom);
                om = (IceGOM*)evt->GetIdDevice(omid);
            }
            if(!om) continue;
            
            //There can be multiple waveforms per OM
            for(I3WaveformSeries::const_iterator waveform_iter = map_iter->second.begin(); waveform_iter != map_iter->second.end(); waveform_iter++){
                
                starttimewaveform = (Double_t)waveform_iter->GetStartTime()/I3Units::nanosecond;
                log_debug("Starttime of waveform is ... ns: %d",starttimewaveform);
                binwidth = (Double_t)waveform_iter->GetBinWidth()/I3Units::nanosecond;
                waveform = (waveform_iter->GetWaveform());
                numberbinswaveform = waveform.size();
                
                if(numberbinswaveform>0){
                    
                    histo.Reset();
                    hname = "OM";
                    hname+=omid;
                    hname+="-CAL-ATWD";
                    histo.SetName(hname.Data());
                    histo.SetBins(numberbinswaveform,starttimewaveform,starttimewaveform+numberbinswaveform*binwidth);
                    
                    ibin = 1;
                    for(vector<double>::iterator wfiter = waveform.begin(); wfiter!=waveform.end(); wfiter++){
                        histo.SetBinContent(ibin,((Double_t)*wfiter)/I3Units::mV);
                        ibin++;
                    }
                    om->SetWaveform(&histo,om->GetNwaveforms()+1);
                }
                
            }// end for(waveform_iter)
            
            // **** Todo: the status information of the waveform, this goes per bin and .GetWaveformInformation() returns a vector
            
        }// end for(map_iter)
    }// end if(wvfrmmapInIceATWD)
    
    //The FADC waveforms
    I3WaveformSeriesMapConstPtr wvfrmmapInIceFADC = frame->Get<I3WaveformSeriesMapConstPtr>(wvfrmMapInIceFADCName_);
    if(wvfrmmapInIceFADC){
        // loop over the OMs
        for(I3WaveformSeriesMap::const_iterator map_iter = wvfrmmapInIceFADC->begin(); map_iter != wvfrmmapInIceFADC->end(); map_iter++){
            if(map_iter->second.size() == 0) continue;
            
            String = map_iter->first.GetString();
            OM = map_iter->first.GetOM();
            omid = icdom.GetOMId(String,OM);
            log_debug("ABOUT TO EXTRACT OM INFO of omid : %i",omid);
            
            //Initializing the OM device in case of absense
            om = (IceGOM*)evt->GetIdDevice(omid);
            if(!om){
                if(String>=79 && String<=90){
                    idom=&dcdom;
                }else if(String<=78){
                    idom=&icdom;
                }
                idom->Reset(1);
                idom->SetUniqueID(omid);
                
                const I3Geometry& geometry = frame->Get<I3Geometry>();
                const I3OMGeoMap& om_geo = geometry.omgeo;
                I3OMGeoMap::const_iterator geo_iter = om_geo.find(map_iter->first);
                assert(geo_iter != om_geo.end());
                const I3Position om_pos = (geo_iter->second).position;
                ompos[0] = (Double_t)om_pos.GetX();
                ompos[1] = (Double_t)om_pos.GetY();
                ompos[2] = (Double_t)om_pos.GetZ();
                omr.SetVector(ompos,"car");
                idom->SetPosition(omr);
                
                evt->AddDevice(idom);
                om = (IceGOM*)evt->GetIdDevice(omid);
            }
            if(!om) continue;
            
            //There can be multiple waveforms per OM
            for(I3WaveformSeries::const_iterator waveform_iter = map_iter->second.begin(); waveform_iter != map_iter->second.end(); waveform_iter++){
                
                starttimewaveform = (Double_t)waveform_iter->GetStartTime()/I3Units::nanosecond;
                log_debug("Starttime of waveform is ... ns: %d",starttimewaveform);
                binwidth = (Double_t)waveform_iter->GetBinWidth()/I3Units::nanosecond;
                waveform = (waveform_iter->GetWaveform());
                numberbinswaveform = waveform.size();
                
                if(numberbinswaveform>0){
                    
                    histo.Reset();
                    hname = "OM";
                    hname+=omid;
                    hname+="-CAL-FADC";
                    histo.SetName(hname.Data());
                    histo.SetBins(numberbinswaveform,starttimewaveform,starttimewaveform+numberbinswaveform*binwidth);
                    
                    ibin = 1;
                    for(vector<double>::iterator wfiter = waveform.begin(); wfiter!=waveform.end(); wfiter++){
                        histo.SetBinContent(ibin,((Double_t)*wfiter)/I3Units::mV);
                        ibin++;
                    }
                    om->SetWaveform(&histo,om->GetNwaveforms()+1);
                }
            }// end for(waveform_iter)
            
            // **** Todo: the status information of the waveform, this goes per bin and .GetWaveformInformation() returns a vector
            
        }// end for(map_iter)
    }//end if(wvfrmmapInIceFADC) loop
    
    //
    // The IceTop Waveforms
    //
    
    // First the IceTOP ATWD Waveforms
    I3WaveformSeriesMapConstPtr wvfrmmapIceTopATWD = frame->Get<I3WaveformSeriesMapConstPtr>(wvfrmMapIceTopATWDName_);
    if(wvfrmmapIceTopATWD){
        // loop over the OMs
        for(I3WaveformSeriesMap::const_iterator map_iter = wvfrmmapIceTopATWD->begin(); map_iter != wvfrmmapIceTopATWD->end(); map_iter++){
            if(map_iter->second.size() == 0) continue;
            
            String = map_iter->first.GetString();
            OM = map_iter->first.GetOM();
            omid = tdom.GetOMId(String,OM);
            log_debug("ABOUT TO EXTRACT OM INFO of omid : %i",omid);
            
            //Initializing the OM device in case of absense
            om = (IceGOM*)evt->GetIdDevice(omid);
            if(!om){
                tdom.Reset(1);
                tdom.SetUniqueID(omid);
                
                const I3Geometry& geometry = frame->Get<I3Geometry>();
                const I3OMGeoMap& om_geo = geometry.omgeo;
                I3OMGeoMap::const_iterator geo_iter = om_geo.find(map_iter->first);
                assert(geo_iter != om_geo.end());
                const I3Position om_pos = (geo_iter->second).position;
                ompos[0] = (Double_t)om_pos.GetX();
                ompos[1] = (Double_t)om_pos.GetY();
                ompos[2] = (Double_t)om_pos.GetZ();
                omr.SetVector(ompos,"car");
                tdom.SetPosition(omr);
                
                evt->AddDevice(tdom);
                om = (IceGOM*)evt->GetIdDevice(omid);
            }
            if(!om) continue;
            
            //There can be multiple waveforms per OM
            for(I3WaveformSeries::const_iterator waveform_iter = map_iter->second.begin(); waveform_iter != map_iter->second.end(); waveform_iter++){
                
                starttimewaveform = (Double_t)waveform_iter->GetStartTime()/I3Units::nanosecond;
                log_debug("Starttime of waveform is ... ns: %d",starttimewaveform);
                binwidth = (Double_t)waveform_iter->GetBinWidth()/I3Units::nanosecond;
                waveform = (waveform_iter->GetWaveform());
                numberbinswaveform = waveform.size();
                
                if(numberbinswaveform>0){
                    
                    histo.Reset();
                    hname = "OM";
                    hname+=omid;
                    hname+="-CAL-ATWD";
                    histo.SetName(hname.Data());
                    histo.SetBins(numberbinswaveform,starttimewaveform,starttimewaveform+numberbinswaveform*binwidth);
                    
                    ibin = 1;
                    for(vector<double>::iterator wfiter = waveform.begin(); wfiter!=waveform.end(); wfiter++){
                        histo.SetBinContent(ibin,((Double_t)*wfiter)/I3Units::mV);
                        ibin++;
                    }
                    om->SetWaveform(&histo,om->GetNwaveforms()+1);
                }
            }// end for(waveform_iter)
            
            // **** Todo: the status information of the waveform, this goes per bin and .GetWaveformInformation() returns a vector
            
        }// end for(map_iter)
    }// end if(wvfrmmapIceTopATWD)
    
    //The FADC waveforms
    I3WaveformSeriesMapConstPtr wvfrmmapIceTopFADC = frame->Get<I3WaveformSeriesMapConstPtr>(wvfrmMapIceTopFADCName_);
    if(wvfrmmapIceTopFADC){
        // loop over the OMs
        for(I3WaveformSeriesMap::const_iterator map_iter = wvfrmmapIceTopFADC->begin(); map_iter != wvfrmmapIceTopFADC->end(); map_iter++){
            if(map_iter->second.size() == 0) continue;
            
            String = map_iter->first.GetString();
            OM = map_iter->first.GetOM();
            omid = tdom.GetOMId(String,OM);
            log_debug("ABOUT TO EXTRACT OM INFO of omid : %i",omid);
            
            //Initializing the OM device in case of absense
            om = (IceGOM*)evt->GetIdDevice(omid);
            if(!om){
                tdom.Reset(1);
                tdom.SetUniqueID(omid);
                
                const I3Geometry& geometry = frame->Get<I3Geometry>();
                const I3OMGeoMap& om_geo = geometry.omgeo;
                I3OMGeoMap::const_iterator geo_iter = om_geo.find(map_iter->first);
                assert(geo_iter != om_geo.end());
                const I3Position om_pos = (geo_iter->second).position;
                ompos[0] = (Double_t)om_pos.GetX();
                ompos[1] = (Double_t)om_pos.GetY();
                ompos[2] = (Double_t)om_pos.GetZ();
                omr.SetVector(ompos,"car");
                tdom.SetPosition(omr);
                
                evt->AddDevice(tdom);
                om = (IceGOM*)evt->GetIdDevice(omid);
            }
            if(!om) continue;
            
            //There can be multiple waveforms per OM
            for(I3WaveformSeries::const_iterator waveform_iter = map_iter->second.begin(); waveform_iter != map_iter->second.end(); waveform_iter++){
                
                starttimewaveform = (Double_t)waveform_iter->GetStartTime()/I3Units::nanosecond;
                log_debug("Starttime of waveform is ... ns: %d",starttimewaveform);
                binwidth = (Double_t)waveform_iter->GetBinWidth()/I3Units::nanosecond;
                waveform = (waveform_iter->GetWaveform());
                numberbinswaveform = waveform.size();
                
                if(numberbinswaveform>0){
                    
                    histo.Reset();
                    hname = "OM";
                    hname+=omid;
                    hname+="-CAL-FADC";
                    histo.SetName(hname.Data());
                    histo.SetBins(numberbinswaveform,starttimewaveform,starttimewaveform+numberbinswaveform*binwidth);
                    
                    ibin = 1;
                    for(vector<double>::iterator wfiter = waveform.begin(); wfiter!=waveform.end(); wfiter++){
                        histo.SetBinContent(ibin,((Double_t)*wfiter)/I3Units::mV);
                        ibin++;
                    }
                    om->SetWaveform(&histo,om->GetNwaveforms()+1);
                }
            }// end for(waveform_iter)
            
            // **** Todo: the status information of the waveform, this goes per bin and .GetWaveformInformation() returns a vector
            
        }// end for(map_iter)
    }//end if(wvfrmmapIceTopFADC) loop
}//end of WaveformExtractor


//Perform the DOMLaunch Raw waveform extraction
//_________________________________________________________________________________________
void IcePacker::DOMLaunchExtractor(I3FramePtr frame)
{
    
    log_debug("ABOUT TO DO DOMLaunchEXTRACTION");
    
    Int_t OM;
    Int_t String;
    Int_t omid;
    TString hname;
    TH1F histo;
    vector<int> dml;
    Double_t starttimedml;
    Int_t numberbinsdml;
    Double_t binwidth = 3; //The sampling rate of ATWD channels is adjustable and taken from CAL-ATWD if available
    Double_t binwidthFADC = 25; //The sampling rate of the FADCs is always 40 MHz.
    Int_t ibin;
    I3DOMLaunchSeriesMapConstPtr domlmap;
    Short_t triggertype;
    Short_t triggermode;
    Short_t whichATWD;
    Nc3Vector omr;
    Double_t ompos[3];
    TH1F* histx=0; // Temp. histogram pointer
    
    //
    // First the In Ice DOMLaunches
    //
    // Reading the I3DOMLaunches from the frame
    domlmap = frame->Get<I3DOMLaunchSeriesMapConstPtr>(domlMapInIceName_);
    if(domlmap){
        // loop over the OMs
        for(I3DOMLaunchSeriesMap::const_iterator map_iter = domlmap->begin(); map_iter != domlmap->end(); map_iter++){
            //  cout<<"In the DOM Launch Serie Map loop"<<endl;
            if(map_iter->second.size() == 0) continue;
            
            log_debug("ABOUT TO EXTRACT OM INFO");
            
            String = map_iter->first.GetString();
            OM=map_iter->first.GetOM();
            omid = icdom.GetOMId(String,OM);
            log_debug("ABOUT TO EXTRACT OM INFO of omid : %i",omid);
            
            //Initializing the OM device in case of absense
            om = (IceGOM*)evt->GetIdDevice(omid);
            if(!om){
                if(String>=79 && String<=90){
                    idom=&dcdom;
                }else if(String<=78){
                    idom=&icdom;
                }
                idom->Reset(1);
                idom->SetUniqueID(omid);
                
                const I3Geometry& geometry = frame->Get<I3Geometry>();
                const I3OMGeoMap& om_geo = geometry.omgeo;
                I3OMGeoMap::const_iterator geo_iter = om_geo.find(map_iter->first);
                assert(geo_iter != om_geo.end());
                const I3Position om_pos = (geo_iter->second).position;
                ompos[0] = (Double_t)om_pos.GetX();
                ompos[1] = (Double_t)om_pos.GetY();
                ompos[2] = (Double_t)om_pos.GetZ();
                omr.SetVector(ompos,"car");
                idom->SetPosition(omr);
                
                evt->AddDevice(idom);
                om = (IceGOM*)evt->GetIdDevice(omid);
            }
            if(!om) continue;
            
            // Obtain binwidth (in ns) from the CAL-ATWD if available
            binwidth=3;
            histx=om->GetWaveform("CAL-ATWD");
            if (histx) binwidth=histx->GetBinWidth(1);
            
            //There can be multiple DOMLaunches per OM
            for(I3DOMLaunchSeries::const_iterator dml_iter = map_iter->second.begin(); dml_iter != map_iter->second.end(); dml_iter++){
                
                starttimedml = (Double_t)dml_iter->GetStartTime()/I3Units::nanosecond;
                log_debug("Starttime of the DOMLaunch is ... ns: %d",starttimedml);
                triggertype =(Short_t)(dml_iter->GetTriggerType());
                triggermode =(Short_t)(dml_iter->GetTriggerMode());
                whichATWD =(Short_t)(dml_iter->GetWhichATWD()); //There's no ATWD selector needed: per DOMLaunch it's either ATWD 1 (a) or ATWD 2 (b), never both!
                
                //**** Todo: add information about LC triggering of indiviual hit
                
                //Extracting the 4 ATWD channel's and the FADC's waveforms and putting it in histos:
                dml = (dml_iter->GetRawATWD(0)); //First the high ATWD
                numberbinsdml = 0;
                numberbinsdml = dml.size();
                if(numberbinsdml>0){
                    histo.Reset();
                    hname = "OM";
                    hname+=omid;
                    hname+="-RAW-ATWD0";
                    histo.SetName(hname.Data());
                    histo.SetBins(numberbinsdml,starttimedml,starttimedml+numberbinsdml*binwidth);
                    ibin = 1;
                    for(vector<int>::iterator wfiter = dml.begin(); wfiter!=dml.end(); wfiter++){
                        histo.SetBinContent(ibin,(Int_t)*wfiter);
                        ibin++;
                    }
                    om->SetWaveform(&histo,om->GetNwaveforms()+1);
                }
                
                dml = (dml_iter->GetRawATWD(1)); //Second the Medium ATWD
                numberbinsdml = 0;
                numberbinsdml = dml.size();
                if(numberbinsdml>0){
                    histo.Reset();
                    hname = "OM";
                    hname+=omid;
                    hname+="-RAW-ATWD1";
                    histo.SetName(hname.Data());
                    histo.SetBins(numberbinsdml,starttimedml,starttimedml+numberbinsdml*binwidth);
                    ibin = 1;
                    for(vector<int>::iterator wfiter = dml.begin(); wfiter!=dml.end(); wfiter++){
                        histo.SetBinContent(ibin,(Int_t)*wfiter);
                        ibin++;
                    }
                    om->SetWaveform(&histo,om->GetNwaveforms()+1);
                }
                
                dml = (dml_iter->GetRawATWD(2)); //Third the Low ATWD
                numberbinsdml = 0;
                numberbinsdml = dml.size();
                if(numberbinsdml>0){
                    histo.Reset();
                    hname = "OM";
                    hname+=omid;
                    hname+="-RAW-ATWD2";
                    histo.SetName(hname.Data());
                    histo.SetBins(numberbinsdml,starttimedml,starttimedml+numberbinsdml*binwidth);
                    ibin = 1;
                    for(vector<int>::iterator wfiter = dml.begin(); wfiter!=dml.end(); wfiter++){
                        histo.SetBinContent(ibin,(Int_t)*wfiter);
                        ibin++;
                    }
                    om->SetWaveform(&histo,om->GetNwaveforms()+1);
                }
                
                dml = (dml_iter->GetRawATWD(3)); //Fourth the Extra ATWD
                numberbinsdml = 0;
                numberbinsdml = dml.size();
                if(numberbinsdml>0){
                    histo.Reset();
                    hname = "OM";
                    hname+=omid;
                    hname+="-RAW-ATWD3";
                    histo.SetName(hname.Data());
                    histo.SetBins(numberbinsdml,starttimedml,starttimedml+numberbinsdml*binwidth);
                    ibin = 1;
                    for(vector<int>::iterator wfiter = dml.begin(); wfiter!=dml.end(); wfiter++){
                        histo.SetBinContent(ibin,(Int_t)*wfiter);
                        ibin++;
                    }
                    om->SetWaveform(&histo,om->GetNwaveforms()+1);
                }
                
                dml = (dml_iter->GetRawFADC()); //Last the FADC
                numberbinsdml = 0;
                numberbinsdml = dml.size();
                if(numberbinsdml>0){
                    histo.Reset();
                    hname = "OM";
                    hname+=omid;
                    hname+="-RAW-FADC";
                    histo.SetName(hname.Data());
                    histo.SetBins(numberbinsdml,starttimedml,starttimedml+numberbinsdml*binwidthFADC);
                    ibin = 1;
                    for(vector<int>::iterator wfiter = dml.begin(); wfiter!=dml.end(); wfiter++){
                        histo.SetBinContent(ibin,(Int_t)*wfiter);
                        ibin++;
                    }
                    om->SetWaveform(&histo,om->GetNwaveforms()+1);
                }
                
                dml = (dml_iter->GetRawChargeStamp()); //Get the RawChargeStamp info ONLY APPROPRIATE FOR IN-ICE DOMs
                numberbinsdml = 0;
                numberbinsdml = dml.size();
                if(numberbinsdml>0){
                    histo.Reset();
                    hname = "OM";
                    hname+=omid;
                    hname+="-RAWCHARGESTAMP";
                    histo.SetName(hname.Data());
                    histo.SetBins(numberbinsdml,starttimedml,starttimedml+numberbinsdml*binwidthFADC);
                    ibin = 1;
                    for(vector<int>::iterator wfiter = dml.begin(); wfiter!=dml.end(); wfiter++){
                        histo.SetBinContent(ibin,(Int_t)*wfiter);
                        ibin++;
                    }
                    om->SetWaveform(&histo,om->GetNwaveforms()+1);
                }
                
                
            }// for ( dml_iter)
        }// for(map_iter)
    }//  if(domlmap)
    
    
    //
    // Now The IceTOP DOMLaunches
    //
    // Reading the I3DOMLaunches from the frame
    domlmap = frame->Get<I3DOMLaunchSeriesMapConstPtr>(domlMapIceTopName_);
    if(domlmap){
        // loop over the OMs
        for(I3DOMLaunchSeriesMap::const_iterator map_iter = domlmap->begin(); map_iter != domlmap->end(); map_iter++){
            //  cout<<"In the DOM Launch Serie Map loop"<<endl;
            if(map_iter->second.size() == 0) continue;
            
            log_debug("ABOUT TO EXTRACT OM INFO");
            
            String = map_iter->first.GetString();
            OM=map_iter->first.GetOM();
            omid = tdom.GetOMId(String,OM);
            log_debug("ABOUT TO EXTRACT OM INFO of omid : %i",omid);
            
            //Initializing the OM device in case of absense
            om = (IceGOM*)evt->GetIdDevice(omid);
            if(!om){
                tdom.Reset(1);
                tdom.SetUniqueID(omid);
                
                const I3Geometry& geometry = frame->Get<I3Geometry>();
                const I3OMGeoMap& om_geo = geometry.omgeo;
                I3OMGeoMap::const_iterator geo_iter = om_geo.find(map_iter->first);
                assert(geo_iter != om_geo.end());
                const I3Position om_pos = (geo_iter->second).position;
                ompos[0] = (Double_t)om_pos.GetX();
                ompos[1] = (Double_t)om_pos.GetY();
                ompos[2] = (Double_t)om_pos.GetZ();
                omr.SetVector(ompos,"car");
                tdom.SetPosition(omr);
                
                evt->AddDevice(tdom);
                om = (IceGOM*)evt->GetIdDevice(omid);
            }
            if(!om) continue;
            
            // Obtain binwidth (in ns) from the CAL-ATWD if available
            binwidth=3;
            histx=om->GetWaveform("CAL-ATWD");
            if (histx) binwidth=histx->GetBinWidth(1);
            
            //There can be multiple DOMLaunches per OM
            for(I3DOMLaunchSeries::const_iterator dml_iter = map_iter->second.begin(); dml_iter != map_iter->second.end(); dml_iter++){
                
                starttimedml = (Double_t)dml_iter->GetStartTime()/I3Units::nanosecond;
                triggertype =(Short_t)(dml_iter->GetTriggerType());
                triggermode =(Short_t)(dml_iter->GetTriggerMode());
                whichATWD =(Short_t)(dml_iter->GetWhichATWD()); //There's no ATWD selector needed: per DOMLaunch it's either ATWD 1 (a) or ATWD 2 (b), never both!
                log_debug("Starttime of the DOMLaunch is ... ns: %d",starttimedml);
                
                //Extracting the 4 ATWD channel's and the FADC's waveforms and putting it in histos:
                dml = (dml_iter->GetRawATWD(0)); //First the high ATWD
                numberbinsdml = 0;
                numberbinsdml = dml.size();
                if(numberbinsdml>0){
                    histo.Reset();
                    hname = "OM";
                    hname+=omid;
                    hname+="-RAW-ATWD0";
                    histo.SetName(hname.Data());
                    histo.SetBins(numberbinsdml,starttimedml,starttimedml+numberbinsdml*binwidth);
                    ibin = 1;
                    for(vector<int>::iterator wfiter = dml.begin(); wfiter!=dml.end(); wfiter++){
                        histo.SetBinContent(ibin,(Int_t)*wfiter);
                        ibin++;
                    }
                    om->SetWaveform(&histo,om->GetNwaveforms()+1);
                }
                
                dml = (dml_iter->GetRawATWD(1)); //Second the Medium ATWD
                numberbinsdml = 0;
                numberbinsdml = dml.size();
                if(numberbinsdml>0){
                    histo.Reset();
                    hname = "OM";
                    hname+=omid;
                    hname+="-RAW-ATWD1";
                    histo.SetName(hname.Data());
                    histo.SetBins(numberbinsdml,starttimedml,starttimedml+numberbinsdml*binwidth);
                    ibin = 1;
                    for(vector<int>::iterator wfiter = dml.begin(); wfiter!=dml.end(); wfiter++){
                        histo.SetBinContent(ibin,(Int_t)*wfiter);
                        ibin++;
                    }
                    om->SetWaveform(&histo,om->GetNwaveforms()+1);
                }
                
                dml = (dml_iter->GetRawATWD(2)); //Third the Low ATWD
                numberbinsdml = 0;
                numberbinsdml = dml.size();
                if(numberbinsdml>0){
                    histo.Reset();
                    hname = "OM";
                    hname+=omid;
                    hname+="-RAW-ATWD2";
                    histo.SetName(hname.Data());
                    histo.SetBins(numberbinsdml,starttimedml,starttimedml+numberbinsdml*binwidth);
                    ibin = 1;
                    for(vector<int>::iterator wfiter = dml.begin(); wfiter!=dml.end(); wfiter++){
                        histo.SetBinContent(ibin,(Int_t)*wfiter);
                        ibin++;
                    }
                    om->SetWaveform(&histo,om->GetNwaveforms()+1);
                }
                
                dml = (dml_iter->GetRawATWD(3)); //Fourth the Extra ATWD
                numberbinsdml = 0;
                numberbinsdml = dml.size();
                if(numberbinsdml>0){
                    histo.Reset();
                    hname = "OM";
                    hname+=omid;
                    hname+="-RAW-ATWD3";
                    histo.SetName(hname.Data());
                    histo.SetBins(numberbinsdml,starttimedml,starttimedml+numberbinsdml*binwidth);
                    ibin = 1;
                    for(vector<int>::iterator wfiter = dml.begin(); wfiter!=dml.end(); wfiter++){
                        histo.SetBinContent(ibin,(Int_t)*wfiter);
                        ibin++;
                    }
                    om->SetWaveform(&histo,om->GetNwaveforms()+1);
                }
                
                dml = (dml_iter->GetRawFADC()); //Last the FADC
                numberbinsdml = 0;
                numberbinsdml = dml.size();
                if(numberbinsdml>0){
                    histo.Reset();
                    hname = "OM";
                    hname+=omid;
                    hname+="-RAW-FADC";
                    histo.SetName(hname.Data());
                    histo.SetBins(numberbinsdml,starttimedml,starttimedml+numberbinsdml*binwidthFADC);
                    ibin = 1;
                    for(vector<int>::iterator wfiter = dml.begin(); wfiter!=dml.end(); wfiter++){
                        histo.SetBinContent(ibin,(Int_t)*wfiter);
                        ibin++;
                    }
                    om->SetWaveform(&histo,om->GetNwaveforms()+1);
                }
                
            }// for ( dml_iter)
        }// for(map_iter)
    }//  if(domlmap)
}// end DOMLaunchExtractor


//Perform the RecoPulse extraction of calibrated and feature extracted hits
//_________________________________________________________________________________________
void IcePacker::RecoPulseExtractor(I3FramePtr frame)
{
    
    log_debug("ABOUT TO DO RecoPulseEXTRACTION");
    
    Int_t OM;
    Int_t String;
    Int_t omid;
    TString pulsesname;
    I3RecoPulseSeriesMapConstPtr recomap;
    Double_t pulsetime;
    Double_t pulsecharge;
    Double_t pulsewidth;
    NcSignal s;
    s.SetSlotName("ADC",1);
    s.SetSlotName("LE",2);
    s.SetSlotName("TOT",3);
    s.SetSlotName("SLC",4);
    s.SetSlotName("HiveSplit",5);//LB:include HiveSplit pulses
    Nc3Vector omr;
    Double_t ompos[3];
    
    vector<vector<Double_t> > pulsesMatrixHS;
    vector<Double_t> pulsesHS;
    
    //LB:First, store HS pulses in an array
    recomap = frame->Get<I3RecoPulseSeriesMapConstPtr>((string)"HiveSplit");
    if(recomap){
        for(I3RecoPulseSeriesMap::const_iterator map_iter = recomap->begin(); map_iter != recomap->end(); map_iter++){
            if(map_iter->second.size() == 0) continue;
            
            log_debug("ABOUT TO EXTRACT HS OM INFO");
            
            String = map_iter->first.GetString();
            OM=map_iter->first.GetOM();
            
            // Skip IceTop DOMs in case of merged pulse series
            if (OM>60) continue;
            
            omid = icdom.GetOMId(String,OM);
            log_debug("ABOUT TO EXTRACT HS OM INFO of omid : %i",omid);
            
            pulsesHS.push_back(omid);
            
            for(I3RecoPulseSeries::const_iterator reco_iter = map_iter->second.begin(); reco_iter != map_iter->second.end(); reco_iter++){
                pulsetime = (Double_t)reco_iter->GetTime()/I3Units::nanosecond;
                
                pulsesHS.push_back(pulsetime);
                
            }//end loop reco_iter
            
            pulsesMatrixHS.push_back(pulsesHS);
            pulsesHS.clear();
            
        }//end loop map_iter
    }// end if recomap
    //LB:--------------------------------------
    
    //
    // First the In Ice RecoPulses
    //
    if(InIcePulsesNames_.size() != 0){// Test whether multiple recopulseseries names are specified
        for ( std::vector<std::string>::iterator ipn = InIcePulsesNames_.begin(); ipn != InIcePulsesNames_.end(); ++ipn ){
            recomap = frame->Get<I3RecoPulseSeriesMapConstPtr>(*ipn);
            if(recomap){
                pulsesname=*ipn;
                // loop over the OMs
                for(I3RecoPulseSeriesMap::const_iterator map_iter = recomap->begin(); map_iter != recomap->end(); map_iter++){
                    if(map_iter->second.size() == 0) continue;
                    
                    log_debug("ABOUT TO EXTRACT OM INFO");
                    
                    String = map_iter->first.GetString();
                    OM=map_iter->first.GetOM();
                    
                    // Skip IceTop DOMs in case of merged pulse series
                    if (OM>60) continue;
                    
                    omid = icdom.GetOMId(String,OM);
                    log_debug("ABOUT TO EXTRACT OM INFO of omid : %i",omid);
                    
                    //Initializing the OM device in case of absense
                    om = (IceGOM*)evt->GetIdDevice(omid);
                    if(!om){
                        if(String>=79 && String<=90){
                            idom=&dcdom;
                        }else if(String<=78){
                            idom=&icdom;
                        }
                        idom->Reset(1);
                        idom->SetUniqueID(omid);
                        
                        const I3Geometry& geometry = frame->Get<I3Geometry>();
                        const I3OMGeoMap& om_geo = geometry.omgeo;
                        I3OMGeoMap::const_iterator geo_iter = om_geo.find(map_iter->first);
                        assert(geo_iter != om_geo.end());
                        const I3Position om_pos = (geo_iter->second).position;
                        ompos[0] = (Double_t)om_pos.GetX();
                        ompos[1] = (Double_t)om_pos.GetY();
                        ompos[2] = (Double_t)om_pos.GetZ();
                        omr.SetVector(ompos,"car");
                        idom->SetPosition(omr);
                        
                        evt->AddDevice(idom);
                        om = (IceGOM*)evt->GetIdDevice(omid);
                    }
                    if(!om) continue;
                    
                    //There can be multiple RecoPulses per OM
                    for(I3RecoPulseSeries::const_iterator reco_iter = map_iter->second.begin(); reco_iter != map_iter->second.end(); reco_iter++){
                        pulsetime = (Double_t)reco_iter->GetTime()/I3Units::nanosecond;
                        pulsecharge =(Double_t)(reco_iter->GetCharge()); // In pe
                        pulsewidth =(Double_t)(reco_iter->GetWidth()/I3Units::nanosecond);
                        log_debug("Time of the RecoPulse is ... ns: %d",pulsetime);
                        
                        s.Reset();
                        s.SetSignal(pulsecharge,1);
                        s.SetSignal(pulsetime,2);
                        s.SetSignal(pulsewidth,3);
                        //@@@	    if(pulsesname.Index("SLC")>=0) s.SetSignal(1,4); // Mark a SLC hit, which is recognised by naming convention of RecoPulseSeries
                        if (!(reco_iter->GetFlags() & reco_iter->LC)) s.SetSignal(1,4); // Mark a SLC hit, which is recognised by the missing HLC bit
                        
                        //LB: Loop over the HS pulses to check if this one is a HS or not
                        for ( std::vector<std::vector<Double_t> >::iterator domHS = pulsesMatrixHS.begin(); domHS != pulsesMatrixHS.end(); ++domHS ){
                            
                            if((Int_t)domHS->front()==omid){
                                
                                for ( std::vector<Double_t>::iterator pulsesHStime = domHS->begin(); pulsesHStime != domHS->end(); ++pulsesHStime ){
                                    
                                    if((Double_t)*pulsesHStime==pulsetime) {
                                        
                                        s.SetSignal(1,5);
                                        break;
                                    }
                                }
                                break;
                            }
                        }
                        //LB:------------------------------------------
                        
                        om->AddHit(s);
                    }// for (reco_iter)
                }// for (map_iter)
            }// if(recomap)
        }// for (ipn)
    }else{// Only one recopulseseries name is specified
        // Reading the I3RecoPulses from the frame
        recomap = frame->Get<I3RecoPulseSeriesMapConstPtr>(InIcePulsesName_);
        if(recomap){
            pulsesname=InIcePulsesName_;
            // loop over the OMs
            for(I3RecoPulseSeriesMap::const_iterator map_iter = recomap->begin(); map_iter != recomap->end(); map_iter++){
                if(map_iter->second.size() == 0) continue;
                
                log_debug("ABOUT TO EXTRACT OM INFO");
                
                String = map_iter->first.GetString();
                OM=map_iter->first.GetOM();
                
                // Skip IceTop DOMs in case of merged pulse series
                if (OM>60) continue;
                
                omid = icdom.GetOMId(String,OM);
                log_debug("ABOUT TO EXTRACT OM INFO of omid : %i",omid);
                
                //Initializing the OM device in case of absense
                om = (IceGOM*)evt->GetIdDevice(omid);
                if(!om){
                    if(String>=79 && String<=90){
                        idom=&dcdom;
                    }else if(String<=78){
                        idom=&icdom;
                    }
                    idom->Reset(1);
                    idom->SetUniqueID(omid);
                    
                    const I3Geometry& geometry = frame->Get<I3Geometry>();
                    const I3OMGeoMap& om_geo = geometry.omgeo;
                    I3OMGeoMap::const_iterator geo_iter = om_geo.find(map_iter->first);
                    assert(geo_iter != om_geo.end());
                    const I3Position om_pos = (geo_iter->second).position;
                    ompos[0] = (Double_t)om_pos.GetX();
                    ompos[1] = (Double_t)om_pos.GetY();
                    ompos[2] = (Double_t)om_pos.GetZ();
                    omr.SetVector(ompos,"car");
                    idom->SetPosition(omr);
                    
                    evt->AddDevice(idom);
                    om = (IceGOM*)evt->GetIdDevice(omid);
                }
                if(!om) continue;
                
                //There can be multiple RecoPulses per OM
                for(I3RecoPulseSeries::const_iterator reco_iter = map_iter->second.begin(); reco_iter != map_iter->second.end(); reco_iter++){
                    pulsetime = (Double_t)reco_iter->GetTime()/I3Units::nanosecond;
                    pulsecharge =(Double_t)(reco_iter->GetCharge()); // In pe
                    pulsewidth =(Double_t)(reco_iter->GetWidth()/I3Units::nanosecond);
                    log_debug("Time of the RecoPulse is ... ns: %d",pulsetime);
                    
                    s.Reset();
                    s.SetSignal(pulsecharge,1);
                    s.SetSignal(pulsetime,2);
                    s.SetSignal(pulsewidth,3);
                    //@@@	  if(pulsesname.Index("SLC")>=0) s.SetSignal(1,4); // Mark a SLC hit, which is recognised by naming convention of RecoPulseSeries
                    if (!(reco_iter->GetFlags() & reco_iter->LC)) s.SetSignal(1,4); // Mark a SLC hit, which is recognised by the missing HLC bit
                   
                    //LB: Loop over the HS pulses to check if this one is a HS or not
                    for ( std::vector<std::vector<Double_t> >::iterator domHS = pulsesMatrixHS.begin(); domHS != pulsesMatrixHS.end(); ++domHS ){
                        
                        if((Int_t)domHS->front()==omid){
                            
                            for ( std::vector<Double_t>::iterator pulsesHStime = domHS->begin(); pulsesHStime != domHS->end(); ++pulsesHStime ){
                                
                                if((Double_t)*pulsesHStime==pulsetime) {
                                    
                                    s.SetSignal(1,5);
                                    break;
                                }
                            }
                            break;
                        }
                    }
                    //LB:------------------------------------------
                    
                    om->AddHit(s);
                }// for (reco_iter)
            }// for (map_iter)
        }// if(recomap)
    }//if (InIcePulsesNames_ != 0)
    
    //
    // Now the IceTOP RecoPulses
    //
    if(IceTopPulsesNames_.size() != 0){// Test whether multiple recopulseseries names are specified
        for ( std::vector<std::string>::iterator ipn = IceTopPulsesNames_.begin(); ipn != IceTopPulsesNames_.end(); ++ipn ){
            recomap = frame->Get<I3RecoPulseSeriesMapConstPtr>(*ipn);
            if(recomap){
                pulsesname=*ipn;
                // loop over the OMs
                for(I3RecoPulseSeriesMap::const_iterator map_iter = recomap->begin(); map_iter != recomap->end(); map_iter++){
                    if(map_iter->second.size() == 0) continue;
                    
                    log_debug("ABOUT TO EXTRACT OM INFO");
                    
                    String = map_iter->first.GetString();
                    OM=map_iter->first.GetOM();
                    
                    // Skip InIce DOMs in case of merged pulse series
                    if (OM<61) continue;
                    
                    omid = tdom.GetOMId(String,OM);
                    log_debug("ABOUT TO EXTRACT OM INFO of omid : %i",omid);
                    
                    //Initializing the OM device in case of absense
                    om = (IceGOM*)evt->GetIdDevice(omid);
                    if(!om){
                        tdom.Reset(1);
                        tdom.SetUniqueID(omid);
                        
                        const I3Geometry& geometry = frame->Get<I3Geometry>();
                        const I3OMGeoMap& om_geo = geometry.omgeo;
                        I3OMGeoMap::const_iterator geo_iter = om_geo.find(map_iter->first);
                        assert(geo_iter != om_geo.end());
                        const I3Position om_pos = (geo_iter->second).position;
                        ompos[0] = (Double_t)om_pos.GetX();
                        ompos[1] = (Double_t)om_pos.GetY();
                        ompos[2] = (Double_t)om_pos.GetZ();
                        omr.SetVector(ompos,"car");
                        tdom.SetPosition(omr);
                        
                        evt->AddDevice(tdom);
                        om = (IceGOM*)evt->GetIdDevice(omid);
                    }
                    if(!om) continue;
                    
                    //There can be multiple RecoPulses per OM
                    for(I3RecoPulseSeries::const_iterator reco_iter = map_iter->second.begin(); reco_iter != map_iter->second.end(); reco_iter++){
                        pulsetime = (Double_t)reco_iter->GetTime()/I3Units::nanosecond;
                        pulsecharge =(Double_t)(reco_iter->GetCharge()); // In pe
                        pulsewidth =(Double_t)(reco_iter->GetWidth()/I3Units::nanosecond);
                        log_debug("Time of the RecoPulse is ... ns: %d",pulsetime);
                        
                        s.Reset();
                        s.SetSignal(pulsecharge,1);
                        s.SetSignal(pulsetime,2);
                        s.SetSignal(pulsewidth,3);
                        //@@@       if(pulsesname.Index("SLC")>=0) s.SetSignal(1,4); // Mark a SLC hit, which is recognised by naming convention of RecoPulseSeries
                        if (!(reco_iter->GetFlags() & reco_iter->LC)) s.SetSignal(1,4); // Mark a SLC hit, which is recognised by the missing HLC bit
                        
                        om->AddHit(s);
                    }// for (reco_iter)
                }// for (map_iter)
            }// if(recomap)
        }//for (ipn)
    }else{// Only one IceTOP recopulseseries name is specified
        // Reading the I3RecoPulses from the frame
        recomap = frame->Get<I3RecoPulseSeriesMapConstPtr>(IceTopPulsesName_);
        if(recomap){
            pulsesname=IceTopPulsesName_;
            // loop over the OMs
            for(I3RecoPulseSeriesMap::const_iterator map_iter = recomap->begin(); map_iter != recomap->end(); map_iter++){
                if(map_iter->second.size() == 0) continue;
                
                log_debug("ABOUT TO EXTRACT OM INFO");
                
                String = map_iter->first.GetString();
                OM=map_iter->first.GetOM();
                
                // Skip InIce DOMs in case of merged pulse series
                if (OM<61) continue;
                
                omid = tdom.GetOMId(String,OM);
                log_debug("ABOUT TO EXTRACT OM INFO of omid : %i",omid);
                
                //Initializing the OM device in case of absense
                om = (IceGOM*)evt->GetIdDevice(omid);
                if(!om){
                    tdom.Reset(1);
                    tdom.SetUniqueID(omid);
                    
                    const I3Geometry& geometry = frame->Get<I3Geometry>();
                    const I3OMGeoMap& om_geo = geometry.omgeo;
                    I3OMGeoMap::const_iterator geo_iter = om_geo.find(map_iter->first);
                    assert(geo_iter != om_geo.end());
                    const I3Position om_pos = (geo_iter->second).position;
                    ompos[0] = (Double_t)om_pos.GetX();
                    ompos[1] = (Double_t)om_pos.GetY();
                    ompos[2] = (Double_t)om_pos.GetZ();
                    omr.SetVector(ompos,"car");
                    tdom.SetPosition(omr);
                    
                    evt->AddDevice(tdom);
                    om = (IceGOM*)evt->GetIdDevice(omid);
                }
                if(!om) continue;
                
                //There can be multiple RecoPulses per OM
                for(I3RecoPulseSeries::const_iterator reco_iter = map_iter->second.begin(); reco_iter != map_iter->second.end(); reco_iter++){
                    pulsetime = (Double_t)reco_iter->GetTime()/I3Units::nanosecond;
                    pulsecharge =(Double_t)(reco_iter->GetCharge()); // In pe
                    pulsewidth =(Double_t)(reco_iter->GetWidth()/I3Units::nanosecond);
                    log_debug("Time of the RecoPulse is ... ns: %d",pulsetime);
                    
                    s.Reset();
                    s.SetSignal(pulsecharge,1);
                    s.SetSignal(pulsetime,2);
                    s.SetSignal(pulsewidth,3);
                    //@@@	  if(pulsesname.Index("SLC")>=0) s.SetSignal(1,4); // Mark a SLC hit, which is recognised by naming convention of RecoPulseSeries
                    if (!(reco_iter->GetFlags() & reco_iter->LC)) s.SetSignal(1,4); // Mark a SLC hit, which is recognised by the missing HLC bit
                    
                    om->AddHit(s);
                }// for (reco_iter)
            }// for (map_iter)
        }// if(recomap)
    }//if (IceTopPulsesNames_ != 0)
}// end RecoPulseExtractor


//Extracts the track information
//____________________________________________________________________________
void IcePacker::TrackExtractor(I3FramePtr frame)
{
    log_debug("ABOUT TO EXTRACT THE TRACKS");
    
    Double_t vec[3];
    NcPosition r;
    Nc3Vector p;
    NcSignal s;
    Float_t pi=acos(-1.);
    Double_t speed;
    Double_t beta;
    Double_t energy;
    
    I3Particle particle;
    I3VectorI3Particle pvec;
    TString name,paramName,paramName2; // NvE: Added paramName2 for flexibility
    I3LineFitParams parLF;
    I3LogLikelihoodFitParams par;
    I3ParaboloidFitParams param;
    //LB: add some parameter types
    CramerRaoParams parCR;
    I3FiniteCuts parFR;
    I3TrackCharacteristicsValues parTCV;
    I3DirectHitsValues parDHV;
    I3Double parDouble;
    vector<string> dHits(4);
    dHits[0]="A";
    dHits[1]="B";
    dHits[2]="C";
    dHits[3]="D";
    Double_t cramraophi=-1.;
    Double_t cramraotheta=-1.;
    Double_t dirLength=-1.;
    //LB:------------
    Int_t nrecotrack=0;
    Int_t nmctrack=0;
    I3FrameObjectConstPtr frame_object;
    
    // Extract the DST_RLogL value (if any) from the frame
    Double_t rlogl=-1.;
    Double_t logl=-1.;
    if (frame->Has("DST_RLogL")) rlogl=frame->Get<I3DoubleConstPtr>("DST_RLogL")->value;
    
    ////////////////////////////////////////
    // Extraction of reconstructed tracks //
    ////////////////////////////////////////
    
    for (I3Frame::typename_iterator iter=frame->typename_begin(); iter!=frame->typename_end(); iter++)
    {
        if (iter->second=="I3Particle")
        {
            name = (TString)iter->first;
            log_debug("Discovered a particle");
            particle = frame->Get<I3Particle>((string)name);
            
            t.Reset();
            
            if (particle.HasPosition())
            {
                vec[0]=particle.GetPos().GetX();
                vec[1]=particle.GetPos().GetY();
                vec[2]=particle.GetPos().GetZ();
                r.SetPosition(vec,"car");
                t.SetBeginPoint(r);
            }
            
            if (particle.HasDirection())
            {
                vec[0]=1; //Set default momentum to 1 GeV
                if (particle.HasEnergy())
                {
                    vec[0]=particle.GetEnergy();
                    if (vec[0]<=0) vec[0]=1;
                }
                vec[1]=fabs(pi-particle.GetDir().GetZenith());
                vec[2]=pi+particle.GetDir().GetAzimuth();
                if (vec[2]>=2.*pi) vec[2]-=2.*pi;
                p.SetVector(vec,"sph");
                t.Set3Momentum(p);
            }
            else if (particle.HasEnergy())
            {
                energy=particle.GetEnergy();
                if (energy>0) t.SetMass(energy);
            }
            
            nrecotrack++;
            t.SetId(nrecotrack);
            t.SetName(name);
            t.SetTitle(name);
            
            // Introduce fit parameters to the track data
            s.Reset(1);
            
            // Add global fit params to all DST tracks
            if (name.Contains("DST"))
            {
                if (rlogl>=0)
                {
                    s.AddNamedSlot("RLogL");
                    s.SetSignal(rlogl,"RLogL");
                }
                speed=particle.GetSpeed();
                if (!isnan(speed))
                {
                    s.AddNamedSlot("Speed");
                    s.SetSignal(speed,"Speed");
                    beta=speed/0.299792458;
                    s.AddNamedSlot("Beta");
                    s.SetSignal(beta,"Beta");
                }
                if (s.GetNslots()) t.SetFitDetails(s);
            }
            
            // Check for track specific fit parameters
            if (name.Contains("LineFit")) // Specific LineFit naming.
            {
                paramName=name+"Params"; //LB: Use paramName istead of changing name because sometimes several parameter keys are necessary.
            }
            else
            {
                paramName=name+"FitParams";
            }
            
            for (I3Frame::typename_iterator iter_param=frame->typename_begin(); iter_param!=frame->typename_end(); iter_param++)
            {
                if (paramName==(TString)iter_param->first && iter_param->second=="I3ParaboloidFitParams")
                {
                    param = frame->Get<I3ParaboloidFitParams>((string)paramName);
                    s.Reset(1);
                    s.SetSlotName("Sigma",1);
                    s.SetSignal(sqrt(param.pbfSigmaZen_*param.pbfSigmaZen_ + param.pbfSigmaAzi_*param.pbfSigmaAzi_),1);
                    s.SetSlotName("Sigma_Zenith",2);
                    s.SetSignal(param.pbfSigmaZen_,2);
                    s.SetSlotName("Sigma_Azimuth",3);
                    s.SetSignal(param.pbfSigmaAzi_,3);
                    s.SetSlotName("LogL",4);
                    s.SetSignal(param.pbfLlh_,4);
                    t.SetFitDetails(s);
                }
                else if (paramName==(TString)iter_param->first && iter_param->second=="I3LogLikelihoodFitParams")
                {
                    par = frame->Get<I3LogLikelihoodFitParams>((string)paramName);
                    s.Reset(1);
                    rlogl=-1.;
                    logl=-1.;
                    s.SetSlotName("RLogL",1);
                    if(!isnan(par.rlogl_)) rlogl=par.rlogl_;
                    s.SetSignal(rlogl,1);
                    s.SetSlotName("LogL",2);
                    if(!isnan(par.logl_)) logl=par.logl_;
                    s.SetSignal(logl,2);
                    s.SetSlotName("NDOF",3);
                    s.SetSignal(par.ndof_,3);
                    s.SetSlotName("Nmini",4);
                    s.SetSignal(par.nmini_,4);
                    //NvE: The following is done below in a generic way for all tracks that have the Cramer-Rao info
/*************
                    //LB: Include more track specific variables.
                    if (name=="MPEFit" || name=="SPEFit2")
                    {
                        paramName=name+"CramerRaoParams";
                        if (frame->Has((string)paramName)){
                            cramraotheta=-1.;
                            cramraophi=-1.;
                            parCR = frame->Get<CramerRaoParams>((string)paramName);
                            s.SetSlotName("CramerRao_theta",5);
                            if(!isnan(parCR.cramer_rao_theta))cramraotheta=parCR.cramer_rao_theta;
                            s.SetSignal(cramraotheta,5);
                            s.SetSlotName("CramerRao_phi",6);
                            if(!isnan(parCR.cramer_rao_phi))cramraophi=parCR.cramer_rao_phi;
                            s.SetSignal(cramraophi,6);
                        }
                    }
                    //LB:-------------------------------
***************/
                    t.SetFitDetails(s);
                }
                else if (paramName==(TString)iter_param->first && iter_param->second=="I3LineFitParams")
                {
                    parLF = frame->Get<I3LineFitParams>((string)paramName);
                    s.Reset(1);
                    s.SetSlotName("Speed",1);
                    s.SetSignal(parLF.LFVel,1);
                    s.SetSlotName("Nhits",2);
                    s.SetSignal(parLF.nHits,2);
                    beta=parLF.LFVel/0.299792458;
                    s.SetSlotName("Beta",3);
                    s.SetSignal(beta,3);
                    t.SetFitDetails(s);
                }
                
                //LB: The MuEx4MPE and FiniteReco cases are different.
                if (name=="FiniteRecoFit")
                {
                    if (frame->Has("FiniteRecoCuts")){
                        parFR = frame->Get<I3FiniteCuts>((string)"FiniteRecoCuts");
                        s.Reset(1);
                        s.SetSlotName("Length",1);
                        s.SetSignal(parFR.Length,1);
                        s.SetSlotName("Lend",2);
                        s.SetSignal(parFR.Lend,2);
                        s.SetSlotName("Lstart",3);
                        s.SetSignal(parFR.Lstart,3);
                        s.SetSlotName("Sdet",4);
                        s.SetSignal(parFR.Sdet,4);
                        s.SetSlotName("finiteCut",5);
                        s.SetSignal(parFR.finiteCut,5);
                        s.SetSlotName("DetectorLength",6);
                        s.SetSignal(parFR.DetectorLength,6);
                        t.SetFitDetails(s);
                    }
                }
               
                if (name=="MuEx4MPE_HiveSplit")
                {
                    s.Reset(1);
//NvE: The following is done in a more generic way for all tracks that contain this info below                    
/***********
                    if (frame->Has("MuEx4MPE_HiveSplitCharacteristics")){
                        parTCV = frame->Get<I3TrackCharacteristicsValues>((string)"MuEx4MPE_HiveSplitCharacteristics");
                        s.AddNamedSlot("AvgDomDistQTotDom");
                        s.SetSignal(parTCV.GetAvgDomDistQTotDom(),"AvgDomDistQTotDom");
                        s.AddNamedSlot("EmptyHitsTrackLength");
                        s.SetSignal(parTCV.GetEmptyHitsTrackLength(),"EmptyHitsTrackLength");
                        s.AddNamedSlot("TrackHitsSeparationLength");
                        s.SetSignal(parTCV.GetTrackHitsSeparationLength(),"TrackHitsSeparationLength");
                        s.AddNamedSlot("TrackHitsDistributionSmoothness");
                        s.SetSignal(parTCV.GetTrackHitsDistributionSmoothness(),"TrackHitsDistributionSmoothness");
                    }
*************/
                    for ( std::vector<string>::iterator dHits_iter = dHits.begin(); dHits_iter != dHits.end(); ++dHits_iter ){
                        
                        if (frame->Has("MuEx4MPE_HiveSplitDirectHits"+(*dHits_iter))){
                            dirLength=-1.;
                            parDHV = frame->Get<I3DirectHitsValues>((string)"MuEx4MPE_HiveSplitDirectHits"+(*dHits_iter));
                            s.AddNamedSlot("NDirStrings"+(*dHits_iter));
                            s.SetSignal(parDHV.GetNDirStrings(),"NDirStrings"+(*dHits_iter));
                            s.AddNamedSlot("NDirDoms"+(*dHits_iter));
                            s.SetSignal(parDHV.GetNDirDoms(),"NDirDoms"+(*dHits_iter));
                            s.AddNamedSlot("NDirPulses"+(*dHits_iter));
                            s.SetSignal(parDHV.GetNDirPulses(),"NDirPulses"+(*dHits_iter));
                            s.AddNamedSlot("QDirPulses"+(*dHits_iter));
                            s.SetSignal(parDHV.GetQDirPulses(),"QDirPulses"+(*dHits_iter));
                            s.AddNamedSlot("NEarlyStrings"+(*dHits_iter));
                            s.SetSignal(parDHV.GetNEarlyStrings(),"NEarlyStrings"+(*dHits_iter));
                            s.AddNamedSlot("NEarlyDoms"+(*dHits_iter));
                            s.SetSignal(parDHV.GetNEarlyDoms(),"NEarlyDoms"+(*dHits_iter));
                            s.AddNamedSlot("NEarlyPulses"+(*dHits_iter));
                            s.SetSignal(parDHV.GetNEarlyPulses(),"NEarlyPulses"+(*dHits_iter));
                            s.AddNamedSlot("QEarlyPulses"+(*dHits_iter));
                            s.SetSignal(parDHV.GetQEarlyPulses(),"QEarlyPulses"+(*dHits_iter));
                            s.AddNamedSlot("NLateStrings"+(*dHits_iter));
                            s.SetSignal(parDHV.GetNLateStrings(),"NLateStrings"+(*dHits_iter));
                            s.AddNamedSlot("NLateDoms"+(*dHits_iter));
                            s.SetSignal(parDHV.GetNLateDoms(),"NLateDoms"+(*dHits_iter));
                            s.AddNamedSlot("NLatePulses"+(*dHits_iter));
                            s.SetSignal(parDHV.GetNLatePulses(),"NLatePulses"+(*dHits_iter));
                            s.AddNamedSlot("QLatePulses"+(*dHits_iter));
                            s.SetSignal(parDHV.GetQLatePulses(),"QLatePulses"+(*dHits_iter));
                            s.AddNamedSlot("DirTrackLength"+(*dHits_iter));
                            if(!isnan(parDHV.GetDirTrackLength())) dirLength=parDHV.GetDirTrackLength();
                            s.SetSignal(dirLength,"DirTrackLength"+(*dHits_iter));
                            s.AddNamedSlot("DirTrackHitDistributionSmoothness"+(*dHits_iter));
                            s.SetSignal(parDHV.GetDirTrackHitDistributionSmoothness(),"DirTrackHitDistributionSmoothness"+(*dHits_iter));
                        }
                    }
                    if (frame->Has("MuEx4MPE_HiveSplit_EnUnc")){
                        parDouble = frame->Get<I3Double>((string)"MuEx4MPE_HiveSplit_EnUnc");
                        s.AddNamedSlot("EnUnc");
                        s.SetSignal(parDouble.value,"EnUnc");
                    }
                    if (frame->Has("MuEx4MPE_HiveSplit_Sigma")){
                        parDouble = frame->Get<I3Double>((string)"MuEx4MPE_HiveSplit_Sigma");
                        s.AddNamedSlot("Sigma");
                        s.SetSignal(parDouble.value,"Sigma");
                    }
                    if (frame->Has("MuEx4MPE_HiveSplit_r")){
                        parDouble = frame->Get<I3Double>((string)"MuEx4MPE_HiveSplit_r");
                        s.AddNamedSlot("r");
                        s.SetSignal(parDouble.value,"r");
                    }
                    if (frame->Has("MuEx4MPE_HiveSplit_rlle")){
                        parDouble = frame->Get<I3Double>((string)"MuEx4MPE_HiveSplit_rlle");
                        s.AddNamedSlot("rlle");
                        s.SetSignal(parDouble.value,"rlle");
                    }
                    if (frame->Has("MuEx4MPE_HiveSplit_rllt")){
                        parDouble = frame->Get<I3Double>((string)"MuEx4MPE_HiveSplit_rllt");
                        s.AddNamedSlot("rllt");
                        s.SetSignal(parDouble.value,"rllt");
                    }
                    t.SetFitDetails(s);
                }
                //LB: ---------------------------------

                //NvE: Add Cramer-Rao parameters for the current track whenever they are available
                paramName2=name+"CramerRaoParams";
                if (frame->Has((string)paramName2))
                {
                 cramraotheta=-1.;
                 cramraophi=-1.;
                 parCR=frame->Get<CramerRaoParams>((string)paramName2);
                 s.AddNamedSlot("CramerRao_theta");
                 if(!isnan(parCR.cramer_rao_theta))cramraotheta=parCR.cramer_rao_theta;
                 s.SetSignal(cramraotheta,"CramerRao_theta");
                 s.AddNamedSlot("CramerRao_phi");
                 if(!isnan(parCR.cramer_rao_phi))cramraophi=parCR.cramer_rao_phi;
                 s.SetSignal(cramraophi,"CramerRao_phi");
                 t.SetFitDetails(s);
                }

                //NvE: Add track characteristics for the current track whenever they are available
                paramName2=name+"Characteristics";
                if (frame->Has((string)paramName2))
                {
                 parTCV=frame->Get<I3TrackCharacteristicsValues>((string)paramName2);
                 s.AddNamedSlot("AvgDomDistQTotDom");
                 s.SetSignal(parTCV.GetAvgDomDistQTotDom(),"AvgDomDistQTotDom");
                 s.AddNamedSlot("EmptyHitsTrackLength");
                 s.SetSignal(parTCV.GetEmptyHitsTrackLength(),"EmptyHitsTrackLength");
                 s.AddNamedSlot("TrackHitsSeparationLength");
                 s.SetSignal(parTCV.GetTrackHitsSeparationLength(),"TrackHitsSeparationLength");
                 s.AddNamedSlot("TrackHitsDistributionSmoothness");
                 s.SetSignal(parTCV.GetTrackHitsDistributionSmoothness(),"TrackHitsDistributionSmoothness");
                 t.SetFitDetails(s);
                }

            } //End of iter_param loop
            
            evt->AddTrack(t);
        }
        else if(iter->second=="I3Vector<I3Particle>")
        {
            name = (TString)iter->first;
            log_debug("Discovered a particle vector.");
            
            if (name.Index("IceDwalk")<0 && name.Index("I3DirectWalk")<0) continue; //Only the DirectWalk tracks!!!
            
            pvec = frame->Get<I3VectorI3Particle>((string)name);
            for(I3VectorI3Particle::const_iterator veciter=pvec.begin(); veciter!=pvec.end(); veciter++)
            {
                particle=*veciter;
                
                t.Reset();
                
                if (particle.HasPosition())
                {
                    vec[0] = particle.GetPos().GetX();
                    vec[1] = particle.GetPos().GetY();
                    vec[2] = particle.GetPos().GetZ();
                    r.SetPosition(vec,"car");
                    t.SetBeginPoint(r);
                }
                
                if (particle.HasDirection())
                {
                    vec[0]=1; //Set default momentum to 1 GeV
                    if (particle.HasEnergy())
                    {
                        vec[0]=particle.GetEnergy();
                        if (vec[0]<=0) vec[0]=1;
                    }
                    vec[1]=fabs(pi-particle.GetDir().GetZenith());
                    vec[2]=pi+particle.GetDir().GetAzimuth();
                    if (vec[2]>=2.*pi) vec[2]-=2.*pi;
                    p.SetVector(vec,"sph");
                    t.Set3Momentum(p);
                }
                else if (particle.HasEnergy())
                {
                    energy=particle.GetEnergy();
                    if (energy>0) t.SetMass(energy);
                }
                
                nrecotrack++;
                t.SetId(nrecotrack);
                t.SetName(name);
                t.SetTitle(name);
                
                evt->AddTrack(t);
            } // End of veciter loop
        }
    } //End of iter loop
    
    log_debug("The total amount of reco tracks equals : %i",nrecotrack);

/////****
//@@@ Temporarily de-activate extraction of simulation tracks because of technical problems
//@@@ due to changed and incompatible coding conventions from icerec-V04 to icerec-V05

    
    /////////////////////////////////////
    // Extraction of simulation tracks //
    /////////////////////////////////////
    
//@@@@    I3MCTreeConstPtr mctree = I3MCTreeUtils::Get(frame,mctreeName_);
    I3MCTreeConstPtr mctree = I3MCTreeUtils::Get(*frame,mctreeName_);
    if (mctree)
    {
        // Extracting the neutrino data in case of nugen files
//@@@@        for(I3MCTree::iterator mciter = mctree->begin(); mciter!=mctree->end(); mciter++)
        for(I3MCTree::const_iterator mciter = mctree->begin(); mciter!=mctree->end(); mciter++)
        {
            name = (TString)mciter->GetTypeString();
            
            if (multimctrack_ || name.Index("NuMu")>=0)
            {
                log_debug("Extracting the Neutrino info");
                
                t.Reset();
                
                if (mciter->HasPosition())
                {
                    vec[0]=mciter->GetPos().GetX();
                    vec[1]=mciter->GetPos().GetY();
                    vec[2]=mciter->GetPos().GetZ();
                    r.SetPosition(vec,"car");
                    t.SetBeginPoint(r);
                }
                
                if (mciter->HasDirection())
                {
                    vec[0]=1; //Set default momentum to 1 GeV
                    if (mciter->HasEnergy())
                    {
                        vec[0]=mciter->GetEnergy();
                        if (vec[0]<=0) vec[0]=1;
                    }
                    vec[1]=fabs(pi-mciter->GetDir().GetZenith());
                    vec[2]=pi+mciter->GetDir().GetAzimuth();
                    if (vec[2]>=2.*pi) vec[2]-=2.*pi;
                    p.SetVector(vec,"sph");
                    t.Set3Momentum(p);
                }
                else if (mciter->HasEnergy())
                {
                    energy=mciter->GetEnergy();
                    if (energy>0) t.SetMass(energy);
                }
                
                nmctrack++;
                t.SetId(-nmctrack); //MC tracks get negative index
                t.SetName(name);
                t.SetTitle(name);
                
                evt->AddTrack(t);
            }
        } // End of mciter loop
        
        if (!multimctrack_) // Only some of the MC tracks are extracted if !multimctrack_
        {
            
            ///////////////////////////////////////////////////////
            // Extracting the most energetic MC primary particle //
            ///////////////////////////////////////////////////////
            
//@@@            if (I3MCTreePhysicsLibrary::GetMostEnergeticPrimary(mctree)!=mctree->end())
            if (I3MCTreePhysicsLibrary::GetMostEnergeticPrimary(mctree))
            {
                log_debug("There is a MostEnergeticPrimary in the MCTree");
                particle = *(I3MCTreePhysicsLibrary::GetMostEnergeticPrimary(mctree));
                
//@@@@                if (particle.IsPrimary() && particle.HasDirection() && particle.HasPosition())
                if (!mctree->parent(particle) && particle.HasDirection() && particle.HasPosition())
                {
                    log_debug("Extracting the MostEnergeticPrimary info");
                    name = "MostEnergeticPrimary";
                    
                    t.Reset();
                    
                    // Get position
                    vec[0]=particle.GetPos().GetX();
                    vec[1]=particle.GetPos().GetY();
                    vec[2]=particle.GetPos().GetZ();
                    
                    r.SetPosition(vec,"car");
                    t.SetBeginPoint(r);
                    
                    // Get momentum (direction)
                    vec[0]=1; //Set default momentum to 1 GeV
                    if (particle.HasEnergy())
                    {
                        vec[0]=particle.GetEnergy();
                        if (vec[0]<=0) vec[0]=1;
                    }
                    vec[1]=fabs(pi-particle.GetDir().GetZenith());
                    vec[2]=pi+particle.GetDir().GetAzimuth();
                    if (vec[2]>=2.*pi) vec[2]-=2.*pi;
                    
                    p.SetVector(vec,"sph");
                    t.Set3Momentum(p);
                    
                    nmctrack++;
                    t.SetId(-nmctrack); //MC tracks get negative index
                    t.SetName(name);
                    t.SetTitle(name);
                    
                    evt->AddTrack(t);
                }
            }
            
            /////////////////////////////////////////
            // Extracting the most energetic track //
            /////////////////////////////////////////
            
//@@@            if (I3MCTreePhysicsLibrary::GetMostEnergeticTrack(mctree)!=mctree->end())
            if (I3MCTreePhysicsLibrary::GetMostEnergeticTrack(mctree))
            {
                log_debug("There is a MostEnergeticTrack in the MCTree");
                particle = *(I3MCTreePhysicsLibrary::GetMostEnergeticTrack(mctree));
                if (particle.IsTrack() && particle.HasDirection() && particle.HasPosition())
                {
                    log_debug("Extracting the MostEnergeticTrack info");
                    name = "MostEnergeticTrack";
                    
                    t.Reset();
                    
                    // Get position
                    vec[0]=particle.GetPos().GetX();
                    vec[1]=particle.GetPos().GetY();
                    vec[2]=particle.GetPos().GetZ();
                    
                    r.SetPosition(vec,"car");
                    t.SetBeginPoint(r);
                    
                    // Get momentum (direction)
                    vec[0]=1; //Set default momentum to 1 GeV
                    if (particle.HasEnergy())
                    {
                        vec[0]=particle.GetEnergy();
                        if (vec[0]<=0) vec[0]=1;
                    }
                    vec[1]=fabs(pi-particle.GetDir().GetZenith());
                    vec[2]=pi+particle.GetDir().GetAzimuth();
                    if (vec[2]>=2.*pi) vec[2]-=2.*pi;
                    
                    p.SetVector(vec,"sph");
                    t.Set3Momentum(p);
                    
                    nmctrack++;
                    t.SetId(-nmctrack); //MC tracks get negative index
                    t.SetName(name);
                    t.SetTitle(name);
                    
                    evt->AddTrack(t);
                }
            }
            
            ///////////////////////////////////////////
            // Extracting the most energetic cascade //
            ///////////////////////////////////////////
            
//@@@@            if (I3MCTreePhysicsLibrary::GetMostEnergeticCascade(mctree)!=mctree->end())
            if (I3MCTreePhysicsLibrary::GetMostEnergeticCascade(mctree))
            {
                log_debug("There is a MostEnergeticCascade in the MCTree");
                particle = *(I3MCTreePhysicsLibrary::GetMostEnergeticCascade(mctree));
                if (particle.IsCascade() && particle.HasDirection() && particle.HasPosition())
                {
                    log_debug("Extracting the MostEnergeticCascade info");
                    name = "MostEnergeticCascade";
                    
                    t.Reset();
                    
                    // Get position
                    vec[0]=particle.GetPos().GetX();
                    vec[1]=particle.GetPos().GetY();
                    vec[2]=particle.GetPos().GetZ();
                    
                    r.SetPosition(vec,"car");
                    t.SetBeginPoint(r);
                    
                    // Get momentum (direction)
                    vec[0]=1; //Set default momentum to 1 GeV
                    if (particle.HasEnergy())
                    {
                        vec[0]=particle.GetEnergy();
                        if (vec[0]<=0) vec[0]=1;
                    }
                    vec[1]=fabs(pi-particle.GetDir().GetZenith());
                    vec[2]=pi+particle.GetDir().GetAzimuth();
                    if (vec[2]>=2.*pi) vec[2]-=2.*pi;
                    
                    p.SetVector(vec,"sph");
                    t.Set3Momentum(p);
                    
                    nmctrack++;
                    t.SetId(-nmctrack); //MC tracks get negative index
                    t.SetName(name);
                    t.SetTitle(name);
                    
                    evt->AddTrack(t);
                }
            }
            
            ////////////////////////////////////////////
            // Extracting the most energetic Neutrino //
            ////////////////////////////////////////////
            
//@@@            if (I3MCTreePhysicsLibrary::GetMostEnergeticNeutrino(mctree)!=mctree->end())
            if (I3MCTreePhysicsLibrary::GetMostEnergeticNeutrino(mctree))
            {
                log_debug("There is a MostEnergeticNeutrino in the MCTree");
                particle = *(I3MCTreePhysicsLibrary::GetMostEnergeticNeutrino(mctree));
                if (particle.IsNeutrino() && particle.HasDirection() && particle.HasPosition())
                {
                    log_debug("Extracting the MostEnergeticNeutrino info");
                    name = "MostEnergeticNeutrino";
                    
                    t.Reset();
                    
                    // Get position
                    vec[0]=particle.GetPos().GetX();
                    vec[1]=particle.GetPos().GetY();
                    vec[2]=particle.GetPos().GetZ();
                    
                    r.SetPosition(vec,"car");
                    t.SetBeginPoint(r);
                    
                    // Get momentum (direction)
                    vec[0]=1; //Set default momentum to 1 GeV
                    if (particle.HasEnergy())
                    {
                        vec[0]=particle.GetEnergy();
                        if (vec[0]<=0) vec[0]=1;
                    }
                    vec[1]=fabs(pi-particle.GetDir().GetZenith());
                    vec[2]=pi+particle.GetDir().GetAzimuth();
                    if (vec[2]>=2.*pi) vec[2]-=2.*pi;
                    
                    p.SetVector(vec,"sph");
                    t.Set3Momentum(p);
                    
                    nmctrack++;
                    t.SetId(-nmctrack); //MC tracks get negative index
                    t.SetName(name);
                    t.SetTitle(name);
                    
                    evt->AddTrack(t);
                }
            }
            
            ////////////////////////////////////////
            // Extracting the most energetic Muon //
            ////////////////////////////////////////
            
//@@@@            if (I3MCTreePhysicsLibrary::GetMostEnergeticMuon(mctree)!=mctree->end())
            if (I3MCTreePhysicsLibrary::GetMostEnergeticMuon(mctree))
            {
                log_debug("There is a MostEnergeticMuon in the MCTree");
                particle = *(I3MCTreePhysicsLibrary::GetMostEnergeticMuon(mctree));
                if (particle.IsTrack() && particle.HasDirection() && particle.HasPosition())
                {
                    log_debug("Extracting the MostEnergeticMuon info");
                    name = "MostEnergeticMuon";
                    
                    t.Reset();
                    
                    // Get position
                    vec[0]=particle.GetPos().GetX();
                    vec[1]=particle.GetPos().GetY();
                    vec[2]=particle.GetPos().GetZ();
                    
                    r.SetPosition(vec,"car");
                    t.SetBeginPoint(r);
                    
                    // Get momentum (direction)
                    vec[0]=1; //Set default momentum to 1 GeV
                    if (particle.HasEnergy())
                    {
                        vec[0]=particle.GetEnergy();
                        if (vec[0]<=0) vec[0]=1;
                    }
                    vec[1]=fabs(pi-particle.GetDir().GetZenith());
                    vec[2]=pi+particle.GetDir().GetAzimuth();
                    if (vec[2]>=2.*pi) vec[2]-=2.*pi;
                    
                    p.SetVector(vec,"sph");
                    t.Set3Momentum(p);
                    
                    nmctrack++;
                    t.SetId(-nmctrack); //MC tracks get negative index
                    t.SetName(name);
                    t.SetTitle(name);
                    
                    evt->AddTrack(t);
                }
            }
            
            ///////////////////////////////////////////
            // Extracting the most energetic Nucleus //
            ///////////////////////////////////////////
            
//@@@@            if (I3MCTreePhysicsLibrary::GetMostEnergeticNucleus(mctree)!=mctree->end())
            if (I3MCTreePhysicsLibrary::GetMostEnergeticNucleus(mctree))
            {
                log_debug("There is a MostEnergeticNucleus in the MCTree");
                particle = *(I3MCTreePhysicsLibrary::GetMostEnergeticNucleus(mctree));
                if (particle.IsNucleus() && particle.HasDirection() && particle.HasPosition())
                {
                    log_debug("Extracting the MostEnergeticNucleus info");
                    name = "MostEnergeticNucleus";
                    
                    t.Reset();
                    
                    // Get position
                    vec[0]=particle.GetPos().GetX();
                    vec[1]=particle.GetPos().GetY();
                    vec[2]=particle.GetPos().GetZ();
                    
                    r.SetPosition(vec,"car");
                    t.SetBeginPoint(r);
                    
                    // Get momentum (direction)
                    vec[0]=1; //Set default momentum to 1 GeV
                    if (particle.HasEnergy())
                    {
                        vec[0]=particle.GetEnergy();
                        if (vec[0]<=0) vec[0]=1;
                    }
                    vec[1]=fabs(pi-particle.GetDir().GetZenith());
                    vec[2]=pi+particle.GetDir().GetAzimuth();
                    if (vec[2]>=2.*pi) vec[2]-=2.*pi;
                    
                    p.SetVector(vec,"sph");
                    t.Set3Momentum(p);
                    
                    nmctrack++;
                    t.SetId(-nmctrack); //MC tracks get negative index
                    t.SetName(name);
                    t.SetTitle(name);
                    
                    evt->AddTrack(t);
                }
            }
        }
    }
    
    log_debug("The total amount of MC tracks equals : %i",nmctrack);

//@@@ End of the temporarily de-activated code
//****/
    
}// End of IcePacker::TrackExtractor

//------------------------------------------------------------------
void IcePacker::Finish()
{
    log_debug("FINISHING");
    
    tfile_->cd();
    ttree_->AutoSave();
    //  tfile_->Write();
    tfile_->Close();
    delete tfile_;
    
    log_info("Wrote %i events to file: %s.",count_,rootfilename_.c_str());
    log_debug("done");
}



