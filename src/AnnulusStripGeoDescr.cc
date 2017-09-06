#include "AnnulusStripGeoDescr.h"
#include <cmath>
#include <iostream>


namespace eutelescope {
 namespace geo {


//AnnulusStripGeoDescr::AnnulusStripGeoDescr( 1025, 1153, 3, 97.779, 106.417, 0.00031, 0.0001932745,0.0001718368,0.02,384.5, 403.481, 427.462, 456.442, 488.423, 438.614,18.981, 23.981, 28.98, 31.981, 9.3660734): 
//        EUTelGenericPixGeoDescr( 97.779, 106.417, 0.00031,//size X, Y, Z (size in mm of sensor) (taken using points Bx & Cx, Ay and rmax-rCentre)
//                                 0, 1025,0,1153, 0, 3,       //minX1, maxX1, minX2, maxX2, minY, maxY number of pixels in x and y
//                                 0.0001932745,0.0001718368,0.02, //pitchPhi1, pitchPhi2, stereo angle (in rads)
//                                 384.5, 403.481, 427.462, 456.442, 488.423, 438.614, //rmin, r1, r2, r3, rmax, rCentre (in mm)
//                                 18.981, 23.981, 28.98, 31.981,// the 4 striplengths (in mm)
//                                 9.3660734 )         //rad length

AnnulusStripGeoDescr::AnnulusStripGeoDescr( int xPixel, int yPixel, double xSize, double ySize, double zSize, double pitchPhi, double stereoAngle, double rmin, double rmax, double rCentre, double radLength):
       EUTelGenericPixGeoDescr( xSize, ySize, zSize,//size X, Y, Z (size in mm of sensor) (taken using points Bx & Cx, Ay and rmax-rCentre)
                                0, xPixel-1, 0, yPixel-1,       //minX1, maxX1, minX2, maxX2, minY, maxY number of pixels in x and y
                                radLength)         //rad length



{
      Double_t PI=3.14159265358979,deg=180/PI;
 
      Double_t phi_i,b,r,c,x,y,gradient,theta;
      //Create the material for the sensor
      matSi = new TGeoMaterial( "Si", 28.0855 , 14.0, 2.33, -_radLength, 45.753206 );
      Si = new TGeoMedium("MimosaSilicon",1, matSi);

      plane = _tGeoManager->MakeBox("sensarea_box",Si,150,500,1);
      //Define the four strip volumes
      rowstrip = _tGeoManager->MakeTubs( "sensarea_strip" , Si, rmin, rmax, zSize/2, 90+(-pitchPhi/2)*deg,90+(pitchPhi/2)*deg);

      //The formula used for calculating strip position is defined in the ATLAS12EC Technical Specs
      theta=pitchPhi*((xPixel-1)/2.0-0.5)+stereoAngle+PI/2;  

    
      if(xPixel<1090){  
      //placement of first two rows
        for( int i = xPixel-1; i >=0; i-- ){
          TGeoCombiTrans* transform=new TGeoCombiTrans(0,0,0,new TGeoRotation("rot",0,0,0));
          //get position of each strip for first row
          phi_i=(i -(xPixel-1)/2.0)*pitchPhi;  
          b=-2*(2*rCentre*sin(stereoAngle/2))*sin(stereoAngle/2+phi_i);
          c=pow((2*rCentre*sin(stereoAngle/2)),2)-pow(rmin,2);
          r=0.5*(-b+sqrt(pow(b,2)-4*c));
          y=r*cos(phi_i+stereoAngle) - rCentre*cos(stereoAngle);
          x=-r*sin(phi_i+stereoAngle) + rCentre*sin(stereoAngle);
          //create first transform
          //rotate to get correct angle of strip
          transform->RotateZ(theta*deg-90);
          transform->SetTranslation(x-rmin*cos(theta),y-rmin*sin(theta),0);  //The R0 center has the coordinate (0, 0)
          //add the nodes
          plane->AddNode(rowstrip,i+1,transform);
          //get angle of next strip
          theta-=pitchPhi;
       }
      }
      else{
      //placement of last two rows
       for( int i = 0; i <=xPixel-1; i++ ){
          TGeoCombiTrans* transform=new TGeoCombiTrans(0,0,0,new TGeoRotation("rot",0,0,0));
          //get position of each strip for second row
          phi_i=((xPixel -1 - i) - (xPixel-1)/2.0)*pitchPhi;  
          b=-2*(2*rCentre*sin(stereoAngle/2))*sin(stereoAngle/2+phi_i);
          c=pow((2*rCentre*sin(stereoAngle/2)),2)-pow(rmin,2);
          r=0.5*(-b+sqrt(pow(b,2)-4*c));
          y=r*cos(phi_i+stereoAngle) - rCentre*cos(stereoAngle);
          x=-r*sin(phi_i+stereoAngle) + rCentre*sin(stereoAngle);
          //create first transform
          //rotate to get correct angle of strip
          transform->RotateZ(theta*deg-90);
          transform->SetTranslation(x-rmin*cos(theta),y-rmin*sin(theta),0);  //The R0 center has the coordinate (0, 0)
          //add the nodes
          plane->AddNode(rowstrip,i+1,transform);
          //get angle of next strip
          theta-=pitchPhi;
       }
 
      }

    }
  
    AnnulusStripGeoDescr::~ AnnulusStripGeoDescr()
    {
      //It appears that ROOT will take ownership and delete that stuff! 
      //delete matSi,
      //delete Si;
    }
  
    void  AnnulusStripGeoDescr::createRootDescr(char const * planeVolume)
    {
      //Get the plane as provided by the EUTelGeometryTelescopeGeoDescription
      TGeoVolume* topplane =_tGeoManager->GetVolume(planeVolume);
      //Add the sensitive area to the plane
      topplane->AddNode(plane, 1,new TGeoTranslation(0,0,0) );
 
    }
  
    std::string AnnulusStripGeoDescr::getPixName(int x , int y){  //needs to be changed to include the information of rows and indice of pixels in x/y direction?
      char buffer [100];
      //return path to the pixel, don't forget to shift indices by +1+
      if(x<512) 
      snprintf( buffer, 100, "/sensarea_box_1/sensarea_strip_%d", x+1);
      else if(x>=512 && x<=_maxIndexX + 256 )
      snprintf( buffer, 100, "/sensarea_box_1/sensarea_strip_%d", x+1-256);
      else
      snprintf( buffer, 100, "invalid index");
      return std::string( buffer ); 
    }
  
  
  
    /*TODO*/ std::pair<int, int>  AnnulusStripGeoDescr::getPixIndex(char const*){return std::make_pair(0,0); }
  
    //EUTelGenericPixGeoDescr* maker()
  //  {
  //    AnnulusStripGeoDescr* mPixGeoDescr = new AnnulusStripGeoDescr();
  //    return dynamic_cast<EUTelGenericPixGeoDescr*>(mPixGeoDescr);
  //  }
  
  } //namespace geo
} //namespace eutelescope



