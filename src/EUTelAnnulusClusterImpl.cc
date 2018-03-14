#include "EUTelAnnulusClusterImpl.h"
#include "EUTelAnnulusPixel.h"

using namespace eutelescope;

EUTelAnnulusClusterImpl::EUTelAnnulusClusterImpl(IMPL::TrackerDataImpl* data): EUTelGenericSparseClusterImpl<EUTelAnnulusPixel>(data)
{
/*nothing else to do*/
} 

EUTelAnnulusClusterImpl::~EUTelAnnulusClusterImpl()
{

}

void EUTelAnnulusClusterImpl::getClusterGeomInfo(float& FxCoord, float& FyCoord, float& xPos, float& yPos, float& rPos, float& phiPos, float& rSize, float& phiSize) const  //get the information of the cluster with several annulus pixels, just for plotting? 
{
	float rMin = std::numeric_limits<float>::max(); 	//stores the largest possible value every pixel will be lower, 
	float phiMin = rMin;					//so its OK for max 
	float rMax = -rMin;								
	float phiMax = rMax;
	float rMinBoundary = 0;
	float rMaxBoundary = 0;
	float phiMinBoundary = 0;
	float phiMaxBoundary = 0;

	EUTelAnnulusPixel* pixel = new EUTelAnnulusPixel;
	//Loop over all the pixels in the cluster
	for( unsigned int index = 0; index < size() ; index++ ) 
	{
		//Get the pixel
		getSparsePixelAt( index , pixel);

		//And its position
                float rlength = pixel->getr();
                float ang = pixel->getang();

		if( ang < phiMin )
		{
			phiMin = ang;
			phiMinBoundary = pixel->getdphi();
		}
		if ( ang > phiMax ) 
		{
			phiMax = ang;
			phiMaxBoundary = pixel->getdphi();
		}
		if ( rlength < rMin )
		{
			rMin = rlength;
			rMinBoundary = pixel->getRmax()-pixel->getRmin();
		}
		if ( rlength > rMax ) 
		{
			rMax = rlength;
			rMaxBoundary = pixel->getRmax()-pixel->getRmin();
		}
	}
        rSize = rMax - rMin + rMinBoundary*0.5 + rMaxBoundary*0.5;  //the r size of the cluster
	phiSize = phiMax - phiMin + phiMaxBoundary*0.5 + phiMinBoundary*0.5; //the width of the angle range of the cluster

	rPos = rMax + rMaxBoundary*0.5 - rSize*0.5;   //the radius of the cluster center
	phiPos = phiMax + phiMaxBoundary*0.5 - phiSize*0.5; // the angle of the cluster center

        xPos = rPos*cos(phiPos) + FxCoord;
        yPos = rPos*sin(phiPos) + FyCoord; 
delete pixel;
}

void EUTelAnnulusClusterImpl::getGeometricCenterOfGravity(float& xCoG, float& yCoG) const
{
	xCoG = 0;
	yCoG = 0;
	
	double totalCharge = 0;
	
	EUTelAnnulusPixel* pixel = new EUTelAnnulusPixel;
	for( unsigned int index = 0; index < size() ; index++ ) 
	{
		getSparsePixelAt( index , pixel);

		double curSignal = pixel->getSignal();
		xCoG += (pixel->getPosX())*curSignal;  //the position in the cartesian coordinate originated in the R0 center
		yCoG += (pixel->getPosY())*curSignal;
		totalCharge += curSignal;
	} 

	xCoG /= totalCharge;
	yCoG /= totalCharge;
	delete pixel;
}

void EUTelAnnulusClusterImpl::getAnnulusCenterOfGravity(float& FxCoord, float& FyCoord, float& xCoG, float& yCoG) const
{
        xCoG = 0;
        yCoG = 0;
        double  rCoG=0;
        double  phiCoG=0;
        double totalCharge = 0;

        EUTelAnnulusPixel* pixel = new EUTelAnnulusPixel;
        for( unsigned int index = 0; index < size() ; index++ )
        {
                getSparsePixelAt( index , pixel);
        
                double curSignal = pixel->getSignal(); 
                rCoG += (pixel->getr())*curSignal;  //the position in the cartesian coordinate originated in the R0 center
                phiCoG += (pixel->getang())*curSignal;
                totalCharge += curSignal;
        }

        rCoG /= totalCharge;
        phiCoG /= totalCharge;
        xCoG = rCoG*cos(phiCoG) + FxCoord;
        yCoG = rCoG*sin(phiCoG) + FyCoord;

 //       streamlog_out(MESSAGE5)<<"cl_size = "<<size()<<" Hit position is "<<xCoG<<" "<<yCoG<<" with r = "<<rCoG<<" and phi = "<<phiCoG<<std::endl;
        delete pixel;
}

void EUTelAnnulusClusterImpl::getHitCovMatrix(double & Phires, double & Rres, double (&cov)[4] ) const
{
      double xCoG = 0;
      double yCoG = 0;

      double totalCharge = 0;
      
      EUTelAnnulusPixel* pixel = new EUTelAnnulusPixel;
      for( unsigned int index = 0; index < size() ; index++ )
      {
          getSparsePixelAt( index , pixel);

          double curSignal = pixel->getSignal();
          double sintheta = sin(pixel->getang());
          double costheta = cos(pixel->getang());
          double radius = pixel->getr();
      	  xCoG += radius*costheta*curSignal;
	  yCoG += radius*sintheta*curSignal;
          cov[0] += Rres*Rres*costheta*costheta*curSignal*curSignal + Phires*Phires*sintheta*sintheta*radius*radius*curSignal*curSignal ;
          cov[1] += Rres*Rres*sintheta*sintheta*curSignal*curSignal + Phires*Phires*costheta*costheta*radius*radius*curSignal*curSignal ;
         // streamlog_out(MESSAGE5)<<"R = "<<radius<<" costheta = "<<costheta<<" sintheta = "<<sintheta<<std::endl; 
	  totalCharge += curSignal;
      }

           
         cov[0] = sqrt(cov[0])/totalCharge; //the resolution of x0
         cov[1] = sqrt(cov[1])/totalCharge; //the resolution of y0
	xCoG /= totalCharge;
	yCoG /= totalCharge;

      for( unsigned int index = 0; index < size() ; index++ )
       {
           getSparsePixelAt( index , pixel);
     
           double curSignal = pixel->getSignal();
           double sintheta = sin(pixel->getang());
           double costheta = cos(pixel->getang());
           double radius = pixel->getr();
           cov[2] += Rres*Rres*(xCoG*costheta+yCoG*sintheta)*(xCoG*costheta+yCoG*sintheta)*curSignal*curSignal + Phires*Phires*(xCoG*sintheta-yCoG*costheta)*(xCoG*sintheta-yCoG*costheta)*radius*radius*curSignal*curSignal ;
           cov[3] += Rres*Rres*(xCoG*sintheta-yCoG*costheta)*(xCoG*sintheta-yCoG*costheta)*curSignal*curSignal + Phires*Phires*(xCoG*costheta+yCoG*sintheta)*(xCoG*costheta+yCoG*sintheta)*radius*radius*curSignal*curSignal ;
       }
 
         cov[2] = sqrt(cov[2])/sqrt(xCoG*xCoG + yCoG*yCoG)/totalCharge;  //the resolution of r
         cov[3] = sqrt(cov[3])/(xCoG*xCoG + yCoG*yCoG)/totalCharge;  //the resolution of phi 
        
       //  streamlog_out(MESSAGE5)<<"InAnnlusCluster res"<<Phires<<",  "<<Rres<<"   "<<cov[0]<<", "<<cov[1]<<", "<<cov[2]<<", "<<cov[3]<<" total charge = "<<totalCharge<<std::endl; 
}
