// eutelescope inlcudes
#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelExceptions.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelTrackerDataInterfacerImpl.h"

// eutelescope geometry
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTelGenericPixGeoDescr.h"

#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <UTIL/CellIDDecoder.h>

#include <algorithm>
#include "marlin/Processor.h"
#include <iostream>

#include "EVENT/TrackerPulse.h"

using namespace eutelescope;
using namespace std;





#include "marlin/Processor.h"

// system includes <>
#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include "EUTelTrack.h"
#include "EUTelHit.h"
#include "EUTelReaderGenericLCIO.h"


class output;
class GBL_trackOutput;
namespace eutelescope {


  class EUTelOutputTTree : public  marlin::Processor{

  public:
    virtual Processor*  newProcessor() { return new EUTelOutputTTree; }

    EUTelOutputTTree();
    ~EUTelOutputTTree(){ delete m_gbl; }
    virtual void init();
    virtual void processRunHeader(LCRunHeader* run);
    virtual void processEvent(LCEvent * evt);
    virtual void check(LCEvent * /*evt*/){ ; };
    virtual void end();

  protected:
    //TbTrack additions
    void prepareTree();
    void clear();

    bool readZsHits(std::string colName, LCEvent* event);
    bool readTracks(LCEvent* event);
    bool readHits(std::string hitColName, LCEvent* event);
    bool first;
    typedef std::map<std::string, output*> output_map_t;
    output_map_t m_out;
    GBL_trackOutput* m_gbl;
    TTree *m_tree;
    TFile *m_file;
    std::string path;
  private:

    

  };
}
class ignoreNames{
    public:
    static bool isIgnored(const std::string & name){
    if (name == "TrackFORselTrack"){
        return true;
    }
    if (name == "StatesFORselTrack")
    {
        return true;
    }
    if (name == "StateHitFORselTrack"){
        return true;
    }
    if(name.find("FOR") != std::string::npos){
        return true;
    }
    return false;
                                                }
};

TFile *gFile_ = NULL;
std::string gStupitNameForShittyROOTFile;
void print_vec(const std::vector<double>& vec, std::ostream& out){

  for (size_t i = 0; i < vec.size(); ++i)
  {
    out << vec[i];
    if (i != vec.size())
    {
      out << ", ";
    }
  }

}

class output{
public:
  virtual void pushCollection(const EVENT::LCCollection* col)=0;

  virtual void newEvent(int eventNR){};
  virtual void FinnishEvent() {};
  virtual void eventEnd();

  virtual void eventStart(int eventNR);
  output(const std::string& name, const std::string& type);
  virtual ~output();
  int m_event_nr;
  TTree *m_tree;
  std::string m_name, m_type;
  std::vector<double> m_x, m_y, m_id;
};

void output::eventEnd()
{

  FinnishEvent();
  if (m_tree){
    m_tree->Fill();
  }
}

void output::eventStart(int eventNR)
{
  m_event_nr = eventNR;
  m_id.clear();
  m_x.clear();
  m_y.clear();
  newEvent(eventNR);
}

output::output(const std::string& name, const std::string& type) :m_event_nr(0), m_name(name), m_type(type)
{
  if (!gFile_)
  {
  //  std::cout << "create file with name: " << gStupitNameForShittyROOTFile << std::endl;
    gFile_ = new TFile(gStupitNameForShittyROOTFile.c_str(), "RECREATE");
  }


  m_tree = new TTree(name.c_str(), name.c_str());
  m_tree->SetDirectory(gFile_->GetDirectory("/"));

  m_tree->Branch("ID", &m_id);
  m_tree->Branch("x", &m_x);
  m_tree->Branch("y", &m_y);
  m_tree->Branch("event_nr", &m_event_nr);
}

output::~output()
{
  if (m_tree)
  {
  //  std::cout << "Write!!! " << m_name << " " << m_type  << std::endl;
    m_tree->Write();
  }
}

class Tracks_output : public output{
public:
  typedef const EVENT::Track* data_t;
  static const char* TypeName(){
    return LCIO::TRACK;
  }
  Tracks_output(const std::string& name, const std::string& type) :output(name, type), warning(true) {}
  virtual void pushCollection(const EVENT::LCCollection* col) {

    if (warning)
    {
      std::cout << "not implemented type: " << TypeName() << std::endl;
      warning = false;
    }
  };
  bool warning;
};

class GBL_trackOutput {
public:
  typedef EUTelTrack * data_t;
  GBL_trackOutput(const std::string& name )  {
    if (!gFile_)
    {
        std::cout << "create file with name: " << gStupitNameForShittyROOTFile << std::endl;
        gFile_ = new TFile(gStupitNameForShittyROOTFile.c_str(), "RECREATE");
    }
    m_tree = new TTree(name.c_str(), name.c_str());
    m_tree->SetDirectory(gFile_->GetDirectory("/"));

    m_tree->Branch("ID", &m_id);
    m_tree->Branch("x", &m_x);
    m_tree->Branch("y", &m_y);
    m_tree->Branch("x_glo", &m_Glox);
    m_tree->Branch("y_glo", &m_Gloy);
    m_tree->Branch("z_glo", &m_Gloz);


    m_tree->Branch("chi2", &m_chi2);
    m_tree->Branch("ndf", &m_ndf);
    m_tree->Branch("phi", &m_phi);
    m_tree->Branch("theta", &m_theta);

    m_tree->Branch("event_nr", &m_event_nr);
    m_tree->Branch("TDC", &m_ptdc);
    m_tree->Branch("Fx_glo", &m_GloFocus_x);
    m_tree->Branch("Fy_glo", &m_GloFocus_y);
    m_tree->Branch("Fz_glo", &m_GloFocus_z);
  }
  static bool hasCollection(EVENT::LCEvent* evt){

    const std::vector<std::string>*  names = evt->getCollectionNames();
    std::vector< std::string >::const_iterator it= std::find(names->begin(), names->end(), "GBL");

    return true;
  }
  ~GBL_trackOutput(){
    if (m_tree){
//        std::cout<<"Write GBL!! " << std::endl; 
      m_tree->Write();
    }
  }
  virtual void pushCollection( EVENT::LCEvent* ev) {
    beginEvent();
    try{
        std::vector<EUTelTrack> tr = lc_reader.getTracks(ev, "tracks");
//        std::cout<< "Track size: " <<tr.size() << std::endl;
        TString strptdc=(ev->getParameters()).getStringVal("PTDC.BIT");
        m_ptdc=strptdc.Atoi();
        for (size_t i = 0; i < tr.size();++i)
        {
          tr[i].print();
          processTrack(tr[i]);
        
        }
    }catch (DataNotAvailableException e) {
		streamlog_out(MESSAGE0) << " collection not available" << std::endl;
	}

    endEvent();
  };
  void beginEvent(){
    m_x.clear();
    m_y.clear();
    m_Glox.clear();
    m_Gloy.clear();
    m_Gloz.clear();
    m_id.clear();
    m_phi.clear();
    m_theta.clear();
    m_ndf.clear();
    m_chi2.clear();
    m_GloFocus_x.clear();
    m_GloFocus_y.clear();
    m_GloFocus_z.clear();
    ++m_event_nr;

  }
  void endEvent(){
    //  std::cout<<"Fill GBL"<<std::endl;
    m_tree->Fill();
  }
  void processTrack(EUTelTrack& trc){
//      std::cout<<"Track" <<std::endl;
   
    std::vector<EUTelState>& planes = trc.getStates();

    chi2 = trc.getChi2();
    ndf = trc.getNdf();
//    std::cout<<"Planes size " << planes.size() << std::endl;
    for (size_t i = 0; i < planes.size(); ++i){
//     std::cout<<"Plane " << i  <<std::endl;
   
      processPlanes(planes[i]);
    }
  }
  void processPlanes(EUTelState& pln){
  
   x = pln.getPosition()[0];
   y = pln.getPosition()[1];
 
   double gloPosition[3];
   geo::gGeometry().local2Master( pln.getLocation(),pln.getPosition(), gloPosition );
   Glox = gloPosition[0];
   Gloy = gloPosition[1];
   Gloz = gloPosition[2];

 
   ID = static_cast<double>(pln.getLocation());
//   std::cout<<"ID " << ID <<std::endl;

   phi = pln.getDirLocalX() / pln.getDirLocalZ();
   theta = pln.getDirLocalY() / pln.getDirLocalZ();
//  std::cout<<"Push back" <<std::endl;
    
   pushHit();
  
    //get the global coordinate of the focus point  
    std::map<int, geo::EUTelAnnulusGear>::iterator mapIt = geo::gGeometry()._AnnulusGearMap.find(pln.getLocation() );
    if( mapIt != geo::gGeometry()._AnnulusGearMap.end() ) {
        geo::EUTelAnnulusGear para = mapIt->second;
        const double local[] = {para.Fx,para.Fy,0};
        double global[3];
        geo::gGeometry().local2Master( pln.getLocation(),local, global );
        m_GloFocus_x.push_back(global[0]);
        m_GloFocus_y.push_back(global[1]);
        m_GloFocus_z.push_back(global[2]);
    } 
    else{
        m_GloFocus_x.push_back(-1);
        m_GloFocus_y.push_back(-1);
        m_GloFocus_z.push_back(-1);

    }

  }

  void pushHit(){
  
    m_x.push_back(x);
    m_y.push_back(y);
    m_Glox.push_back(Glox);
    m_Gloy.push_back(Gloy);
    m_Gloz.push_back(Gloz);
    m_id.push_back(ID);
    m_phi.push_back(phi);
    m_theta.push_back(theta);

    m_ndf.push_back(ndf);
    m_chi2.push_back(chi2);

  }


  double x;
  double y;
  double Glox;
  double Gloy;
  double Gloz;

  double ID;

  double phi;
  double theta;



  double chi2;
  double ndf;
  EUTelReaderGenericLCIO lc_reader;
  int m_event_nr;
  TTree *m_tree;
  std::string m_name, m_type;
  std::vector<double> m_x, m_y, m_Glox, m_Gloy, m_Gloz, m_id,
    m_chi2,m_ndf,m_phi,m_theta;
  int m_ptdc;
  std::vector<double> m_GloFocus_x, m_GloFocus_y, m_GloFocus_z;
};


class LCGenericObject_output : public output{
public:
  typedef const void* data_t;
  static const char* TypeName(){
    return "LCGenericObject";
  }
  LCGenericObject_output(const std::string& name, const std::string& type) :output(name, type), warning(true) {}
  virtual void pushCollection(const EVENT::LCCollection* col) {

    if (warning)
    {
      std::cout << "not implemented type: " << TypeName() << std::endl;
      warning = false;
    }
  };
  
  bool warning;
};


class Track_output :public output{
public:
  typedef const EVENT::Track* data_t;
  static const char* TypeName(){
    return LCIO::TRACK;
  }
  virtual void newEvent(int eventNR){

    m_z.clear();
    m_D0.clear();
    m_phi.clear();
    m_omega.clear();
  }
  Track_output(const std::string& name, const std::string& type) :output(name, type) {
    m_tree->Branch("z", &m_z);
    m_tree->Branch("D0", &m_D0);
    m_tree->Branch("phi", &m_phi);
    m_tree->Branch("omega", &m_omega);

  }
  virtual void pushCollection(const EVENT::LCCollection* col) {

    for (int i = 0; i < col->getNumberOfElements(); ++i)
    {

      data_t  trk = dynamic_cast<data_t>(col->getElementAt(i));

      if (!trk)
      {
        std::cout << "unable to cast the pointer in " << m_type << ": " << m_name << std::endl;
        return;
      }
      m_id.push_back(0);
      m_x.push_back(trk->getReferencePoint()[0]);
      m_y.push_back(trk->getReferencePoint()[1]);
      m_z.push_back(trk->getReferencePoint()[2]);
      m_D0.push_back(trk->getD0());
      m_phi.push_back(trk->getPhi());
      m_omega.push_back(trk->getOmega());
    }
  }

  std::vector<double> m_z, m_D0, m_phi, m_omega;
};


class TrackerData_output :public output{
public:
  typedef const EVENT::TrackerData data_t;
  static const char* TypeName(){
    return LCIO::TRACKERDATA;
  }
  virtual void pushCollection(const EVENT::LCCollection* col) {
  
    for (size_t i = 0; i < col->getNumberOfElements(); ++i){



      data_t* hit = dynamic_cast<data_t*> (col->getElementAt(i));
      if (!hit)
      {
        std::cout << "unable to cast the pointer in " << m_type << ": " << m_name << std::endl;
        return;
      }

      for (size_t i = 0; i < hit->getChargeValues().size() / 4; i++)
      {

        UTIL::CellIDDecoder<data_t> hitDecoder(col);
        int sensorID = hitDecoder(hit)["sensorID"];
        m_id.push_back(sensorID);
        m_x.push_back(hit->getChargeValues().at(i * 4));
        m_y.push_back(hit->getChargeValues().at(i * 4 + 1));
      }

    }
  }

  TrackerData_output(const std::string& name, const std::string& type) :output(name,type){ }
};


class TrackerPulse_Output :public output{
public:
  typedef const EVENT::TrackerPulse data_t;
  static const char* TypeName(){
    return LCIO::TRACKERPULSE;
  }

  virtual void pushCollection(const EVENT::LCCollection* col){

    for (size_t i = 0; i < col->getNumberOfElements(); ++i){



      data_t* hit = dynamic_cast<data_t*> (col->getElementAt(i));
      if (!hit)
      {
        std::cout << "unable to cast the pointer" << std::endl;
        return;
      }

      for (size_t i = 0; i < hit->getTrackerData()->getChargeValues().size() / 4; i++)
      {

        UTIL::CellIDDecoder<data_t> hitDecoder(col);
        int sensorID = hitDecoder(hit)["sensorID"];
        m_id.push_back(sensorID);

        m_x.push_back(hit->getTrackerData()->getChargeValues().at(i * 4));
        m_y.push_back(hit->getTrackerData()->getChargeValues().at(i * 4 + 1));
      }

    }
  }
  TrackerPulse_Output(const std::string& name, const std::string& type):output(name,type){
  }
  virtual ~TrackerPulse_Output(){}

};

class TrackerHit_output :public output{
public:
  typedef const EVENT::TrackerHit data_t;
  static const char* TypeName(){
    return LCIO::TRACKERHIT;
  }

  virtual void newEvent(int eventNR){
     m_id_2.clear();
     m_x_2.clear();
     m_y_2.clear();
     m_Fx_2.clear();
     m_Fy_2.clear();
  }

  virtual void pushCollection(const EVENT::LCCollection* col){

    data_t* hit = NULL;
    for (size_t i = 0; i < col->getNumberOfElements(); ++i){

      hit = dynamic_cast<data_t*>(col->getElementAt(i));
      if (!hit){
        std::cout << "unable to cast the pointer" << std::endl;
        continue;

      }
      UTIL::CellIDDecoder<data_t> hitDecoder(col);
      int sensorID = hitDecoder(hit)["sensorID"];
      m_id.push_back(sensorID);
      m_x.push_back(hit->getPosition()[0]);
      m_y.push_back(hit->getPosition()[1]);
      if(m_name == "local_hit")
      { 
          double id_2 = sensorID;
          if(sensorID>=10&&sensorID<=21){
             id_2 = sensorID+12;
          } else if (sensorID>=22&&sensorID<=33){
             id_2 = sensorID-12;
          }
                    
          double gloPosition[3];
          geo::gGeometry().local2Master( sensorID,hit->getPosition(), gloPosition );
          double localPosition[2];
          geo::gGeometry().master2Local( id_2, gloPosition, localPosition );
          
          m_id_2.push_back(id_2); 
          m_x_2.push_back(localPosition[0]); 
          m_y_2.push_back(localPosition[1]);

          std::map<int, geo::EUTelAnnulusGear>::iterator mapIt = geo::gGeometry()._AnnulusGearMap.find(sensorID );
          if( mapIt != geo::gGeometry()._AnnulusGearMap.end() ) {
             geo::EUTelAnnulusGear para = mapIt->second;
             const double local[] = {para.Fx,para.Fy,0};
             double glofocus[3], localfocus_2[2];
             geo::gGeometry().local2Master( sensorID,local, glofocus );
             geo::gGeometry().master2Local( id_2, glofocus, localfocus_2);
             m_Fx_2.push_back(localfocus_2[0]);
             m_Fy_2.push_back(localfocus_2[1]);
          } else{
             m_Fx_2.push_back(-1);
             m_Fy_2.push_back(-1);
          }
 
      }
    }
  }

  TrackerHit_output(const std::string& name, const std::string& type) :output(name, type){
     if(name =="local_hit"){
        m_tree->Branch("ID_2", &m_id_2);
        m_tree->Branch("x_2", &m_x_2);
        m_tree->Branch("y_2", &m_y_2);
        m_tree->Branch("Fx_2", &m_Fx_2);
        m_tree->Branch("Fy_2", &m_Fy_2);
     } 
  }
  
  std::vector<double> m_id_2, m_x_2, m_y_2, m_Fx_2, m_Fy_2;
};

output* createOutput(const std::string& name, const std::string& type){
 
  if (ignoreNames::isIgnored(name))
  {
    return NULL;
  }

  if (type == Track_output::TypeName()){
    return new Track_output(name, type);
  }

  if (type == TrackerPulse_Output::TypeName()){

    return new TrackerPulse_Output(name, type);
  }
  if (type == TrackerData_output::TypeName())
  {
    return new  TrackerData_output(name, type);
  }
  if (type==TrackerHit_output::TypeName())
  {
    return new TrackerHit_output(name, type);
  }
  if (type==LCGenericObject_output::TypeName())
  {
    return new LCGenericObject_output(name, type);
  }
// if (type == TrackerHit_output::TypeName())
// {
//    return new TrackerHit_output(name, type);
//  }

  std::cout << "unsupported type  " << type << std::endl;
  return NULL;
}


EUTelOutputTTree::EUTelOutputTTree() :Processor("EUTelOutputTTree"), m_file(NULL), m_tree(NULL) , m_gbl(NULL){

  
  registerProcessorParameter("OutputPath", "Path/File where root-file should be stored",
    path, std::string("NTuple.root"));
}
void EUTelOutputTTree::init()
{
 
  std::string name("test.root");
  geo::gGeometry().initializeTGeoDescription(name, true);
  first = true;


  gStupitNameForShittyROOTFile = path;




  std::cout << "EUTelOutputTTree::init" << std::endl;
}

void EUTelOutputTTree::processRunHeader(LCRunHeader* run)
{

  std::cout << "EUTelOutputTTree::processRunHeader" << std::endl;
}

void EUTelOutputTTree::processEvent(LCEvent * evt)
{
  const std::vector<std::string>*  names = evt->getCollectionNames();
  std::vector< std::string >::const_iterator name;
    for (name = names->begin(); name != names->end(); name++){
//          std::cout<< "Collection name: "<< *name <<std::endl;
      if (ignoreNames::isIgnored(*name))
      {
//          std::cout<< "Collection name Ignored: "<< *name <<std::endl;

        continue;
      }
      LCCollection* col = evt->getCollection(*name);
      if (!col)
      {
        continue;
      }
      output_map_t::const_iterator  it = m_out.find(*name);
      if (it == m_out.end()){
        m_out[*name] = createOutput(*name, col->getTypeName());

      }
      m_out[*name]->eventStart(evt->getEventNumber());

      m_out[*name]->pushCollection(col);

    }
    for (output_map_t::const_iterator it = m_out.begin(); it != m_out.end(); ++it){
      it->second->eventEnd();
    }

    if (GBL_trackOutput::hasCollection(evt))
    {
      if (!m_gbl)
      {
        m_gbl = new GBL_trackOutput("GBL_tracks");
      }
      
      m_gbl->pushCollection(evt);

    }

}

void EUTelOutputTTree::end()
{
  for (output_map_t::const_iterator it = m_out.begin(); it != m_out.end(); ++it){
    delete it->second;
  }
  m_out.clear();
  if (m_gbl)
  {
    delete m_gbl;
    m_gbl = NULL;
    std::cout << "closing  GBL " << std::endl;

  }
  if (gFile_)
  {
    std::cout << "closing  file " << std::endl;
    gFile_->Write();
    gFile_->Close();
    delete gFile_;
    gFile_ = NULL;
  }


  std::cout << "EUTelOutputTTree::end" << std::endl;
}


EUTelOutputTTree gEUTelOutputTTree;
