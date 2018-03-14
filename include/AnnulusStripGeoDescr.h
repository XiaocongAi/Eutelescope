#ifndef ANNULUSSTRIPGE0DESCR_H
#define ANNULUSSTRIPGE0DESCR_H

  /** @class AnnulusStripGeoDescr 
    */

#include <string> //std::string
#include <utility> //std::pair
 
 //EUTELESCOPE
#include "EUTelGenericPixGeoDescr.h"
 
 //ROOT
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"



namespace eutelescope {
namespace geo {

class AnnulusStripGeoDescr: public EUTelGenericPixGeoDescr {

        public:
               AnnulusStripGeoDescr( int xPixel, int yPixel, double xSize, double ySize, double zSize, double pitchPhi, double stereoAngle, double rmin, double rmax, double rCentre, int order, double radLength);
               ~AnnulusStripGeoDescr();
        
               void createRootDescr(char const *);
               std::string getPixName(int, int);
               std::pair<int, int> getPixIndex(char const *);

              void getPitchesAndStereo(double& pitchPhi,  double& stereoAngle){
                      pitchPhi= _pitchPhi;
                      stereoAngle= _stereoAngle;
              }


              void getRadii(double& rMin, double& rMax, double& rCentre){
                      rMin= _rMin;
                      rMax= _rMax;
                      rCentre= _rCentre;
              }

             void getstripLength(double & stripLength){
                     stripLength = _stripLength;
             }

             void getstripOrder(int & stripOrder){
                     stripOrder = _stripOrder;
             }
 
 
        protected:

                TGeoMaterial* matSi;
		TGeoMedium* Si;
		TGeoVolume* plane,*rowstrip;
                double _pitchPhi, _stereoAngle;
                double _rMin, _rMax, _rCentre;
                double _stripLength;
                int _stripOrder;

};

}
}
#endif
