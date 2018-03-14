#include "AnnulusStripGeoDescr.h"
#include <cmath>
#include <iostream>


namespace eutelescope {
 namespace geo {

AnnulusStripGeoDescr::AnnulusStripGeoDescr( int xPixel, int yPixel, double xSize, double ySize, double zSize, double pitchPhi, double stereoAngle, double rmin, double rmax, double rCentre, int order, double radLength):
       EUTelGenericPixGeoDescr( xSize, ySize, zSize,
                                0, xPixel-1, 0, yPixel-1,    
                                radLength), _pitchPhi(pitchPhi), _stereoAngle(stereoAngle), _rMin(rmin), _rMax(rmax), _rCentre(rCentre), _stripOrder(order), _stripLength(rmax-rmin)      
     {
      
      Double_t PI=3.14159265358979,deg=180/PI;
  
      Double_t phi_i,b,r,c,x,y,gradient,theta;
      
      //Create the material for the sensor
      matSi = new TGeoMaterial( "Si", 28.0855 , 14.0, 2.33, -_radLength, 45.753206 );
      Si = new TGeoMedium("MimosaSilicon",1, matSi);

      plane = _tGeoManager->MakeBox("sensarea_box",Si,150,500,1);
      rowstrip = _tGeoManager->MakeTubs( "sensarea_strip" , Si, rmin, rmax, zSize/2, 90+(-pitchPhi/2)*deg,90+(pitchPhi/2)*deg);

      //The formula used for calculating strip position is defined in the ATLAS12EC Technical Specs
      theta=pitchPhi*((xPixel-1)/2.0)+stereoAngle+PI/2;  
    
      if(order==-1){  
         //strip index counting from right to left
         for( int i = xPixel-1; i >=0; i-- ){
             TGeoCombiTrans* transform=new TGeoCombiTrans(0,0,0,new TGeoRotation("rot",0,0,0));
             //get position of each strip
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
      else if (order ==1) {
         //strip index counting from left to right
         for( int i = 0; i <=xPixel-1; i++ ){
            TGeoCombiTrans* transform=new TGeoCombiTrans(0,0,0,new TGeoRotation("rot",0,0,0));
            //get position of each strip
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
      else std::cout<<"Strip index order must be 1 or -1! Please check!"<<std::endl;
      

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
  
  } //namespace geo
} //namespace eutelescope



