#ifndef ICEPACKER_H
#define ICEPACKER_H
/**
 * class: IcePacker
 *
 * Version $Id: $
 *
 * date: $Date: $
 *
 * @author duvoort
 *
 * (c) 2009 IceCube Collaboration
 */

// header files
#include "../IcePackBase/IPBHeaders.h"
#include "icetray/I3Module.h"
#include "icetray/I3Frame.h"
#include <TTree.h>
#include <TFile.h>
#include <vector>

using namespace std;

typedef struct{
  int pretrig;
}ITrigg;

///
/// The module that writes the trees
///
class IcePacker : public I3Module
{
public:
  IcePacker();
  IcePacker(const IcePacker&);
  IcePacker(const I3Context&);
  IcePacker& operator=(const IcePacker&);
  ~IcePacker() {};

  void Configure();
  void Reconfigure();
  void Physics(I3FramePtr);
  void Finish();

private:
  TTree* ttree_;
  TFile* tfile_;
  IceEvent* evt;
  NcDevice daq;
  NcDevice mcinfo;
  NcDevice trig;
  IceICDOM icdom;
  IceDCDOM dcdom;
  IceIDOM* idom;
  IceTDOM tdom;
  IceGOM* om;
  NcTrack t;
  string rootfilename_;
  string treename_;
  string domlMapInIceName_;
  string domlMapIceTopName_;
  string wvfrmMapInIceATWDName_;
  string wvfrmMapInIceFADCName_;
  string wvfrmMapIceTopATWDName_;
  string wvfrmMapIceTopFADCName_;
  string InIcePulsesName_;
  string IceTopPulsesName_;
  std::vector<std::string> InIcePulsesNames_;
  std::vector<std::string> IceTopPulsesNames_;
  std::vector<std::string> fSubeventSels;
  std::vector<std::string> fFilterSels;
  Bool_t fFilterPrescale;
  Bool_t fSelNamePattern;
  string eventheaderName_;
  string triggerName_;
  string filterName_;
  string mctreeName_;
  string mcweightinfo_;
  Int_t count_;
  Bool_t trackextractOn_;
  Bool_t multimctrack_;
  long maxtreesize;

 protected:
 
  void ReadTrigger(I3FramePtr);
  void WaveformExtractor(I3FramePtr);
  void DOMLaunchExtractor(I3FramePtr);
  void RecoPulseExtractor(I3FramePtr);
  void TrackExtractor(I3FramePtr);

  SET_LOGGER("IcePacker");

};


#endif //ICEPACKER_H
