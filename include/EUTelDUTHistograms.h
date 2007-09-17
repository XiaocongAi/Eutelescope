// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTelDUTHistograms_h
#define EUTelDUTHistograms_h 1

#include "marlin/Processor.h"

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// lcio includes <.h> 
#include "lcio.h"

// AIDA includes <.h>
#ifdef MARLIN_USE_AIDA
#include <AIDA/IBaseHistogram.h>
#include <AIDA/IHistogram1D.h>
#endif

// system includes <>
#include <string>
#include <vector>
#include <map>

namespace eutelescope {


  //! DUT analysis processor for EUDET Telescope
  /*! This processor was designed for analysis of DUT performancs
   *  based on the analytic track fitting results.
   *
   * \par Geometry description
   * Geometry information is taken from GEAR.
   *
   * \par Input  
   * \c Track collection with fit results and \c TrackerHit collection
   * with DUT hits are taken as an input. 
   * 
   * \param InputTrackCollectionName  Name of the input Track collection
   *
   * \param InputHitCollectionName  Name of the input TrackerHit collection,
   *  from which DUT hits are taken 
   *
   * \param DUTalignment Alignment corrections for DUT: shift in X, Y
   *                     and rotation around Z
   *
   * \param DistMax Maximum allowed distance between fit and matched
   *                 DUT hit.
   *
   * \param HistoInfoFileName Name of the histogram information file.
   *        Using this file histogram parameters can be changed without 
   *        recompiling the code.
   *
   * \param DebugEventCount      Print out debug and information
   * messages only for one out of given number of events. If zero, no
   * debug information is printed. 
   *
   * \author A.F.Zarnecki, University of Warsaw
   * @version $Id: EUTelDUTHistograms.h,v 1.0 2007-09-17 22:27:25 zarnecki Exp $
   * \date 2007.09.15
   *
   */ 


  class EUTelDUTHistograms : public marlin::Processor {
  
  public:

  
     
    //! Returns a new instance of EUTelDUTHistograms
    /*! This method returns an new instance of the this processor.  It
     *  is called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *  
     *  @return a new EUTelDUTHistograms
     */
    virtual Processor*  newProcessor() { return new EUTelDUTHistograms ; }
  
    //! Default constructor 
    EUTelDUTHistograms() ;
  
    //! Called at the job beginning.
    /*! This is executed only once in the whole execution. 
     *  
     */
    virtual void init() ;
  
    //! Called for every run.
    /*!
     * @param run the LCRunHeader of the current run
     */
    virtual void processRunHeader( LCRunHeader* run ) ;
  
    //! Called every event
    /*! This is called for each event in the file.
     * 
     *  @param evt the current LCEvent event 
     */
    virtual void processEvent( LCEvent * evt ) ; 
  
    //! Check event method
    /*! This method is called by the Marlin execution framework as
     *  soon as the processEvent is over. It can be used to fill check
     *  plots. For the time being there is nothing to check and do in
     *  this slot.
     * 
     *  @param evt The LCEvent event as passed by the ProcessMgr
     */
    virtual void check( LCEvent * evt ) ; 
  
  
    //! Book histograms
    /*! This method is used to books all required
     *  histograms. Histogram pointers are stored into
     *  _aidaHistoMap so that they can be recalled and filled 
     * from anywhere in the code.  
     */
    void bookHistos();


    //! Called after data processing for clean up.
    /*! Used to release memory allocated in init() step
     */
    virtual void end() ;
  
  protected:

    //! Silicon planes parameters as described in GEAR
    /*! This structure actually contains the following:
     *  @li A reference to the telescope geoemtry and layout
     *  @li An integer number saying if the telescope is w/ or w/o DUT
     *  @li An integer number saying the number of planes in the
     *  telescope.
     *
     *  This object is provided by GEAR during the init() phase and
     *  stored here for local use.
     */ 
    gear::SiPlanesParameters * _siPlanesParameters;

    //! Silicon plane layer layout
    /*! This is the real geoemetry description. For each layer
     *  composing the telescope the relevant information are
     *  available.
     *  
     *  This object is taken from the _siPlanesParameters during the
     *  init() phase and stored for local use
     */ 
    gear::SiPlanesLayerLayout * _siPlanesLayerLayout;

    //! The histogram information file
    /*! This string contain the name of the histogram information
     *  file. This is selected by the user in the steering file.
     * 
     *  @see eutelescope::EUTelHistogramManager
     *  @see eutelescope::EUTelHistogramInfo
     */ 
    std::string _histoInfoFileName;


    //! Input \c Track collection name
    std::string _inputTrackColName ;

    //! Input \c TrackerHit collection name
    std::string _inputHitColName ;

    //!  Debug print out for one out of given number of events.
    int _debugCount ;



    // Internal processor variables
    // ----------------------------


    int _nRun ;
    int _nEvt ;

    int _iDUT;
    double _zDUT;  
    double _distMax;

    std::vector<double> _measuredX;
    std::vector<double> _measuredY;

    std::vector<double> _fittedX;
    std::vector<double> _fittedY;

    std::vector<float > _DUTalign;


 #ifdef MARLIN_USE_AIDA
    //! AIDA histogram map
    /*! Used to refer to histograms by their names, i.e. to recall 
     *  a histogram pointer using histogram name.
     */

    std::map<std::string , AIDA::IBaseHistogram * > _aidaHistoMap;

    static std::string _MeasuredXHistoName;
    static std::string _MeasuredYHistoName;
    static std::string _MeasuredXYHistoName;

    static std::string _FittedXHistoName;
    static std::string _FittedYHistoName;
    static std::string _FittedXYHistoName;

    static std::string _EfficiencyXHistoName;
    static std::string _EfficiencyYHistoName;
    static std::string _EfficiencyXYHistoName;

    static std::string _ShiftXHistoName;
    static std::string _ShiftYHistoName;
    static std::string _ShiftXYHistoName;
    static std::string _ShiftXvsYHistoName;
    static std::string _ShiftYvsXHistoName;

#endif 

 } ;

  
  //! A global instance of the processor.
  EUTelDUTHistograms aEUTelDUTHistograms ;


}

#endif



