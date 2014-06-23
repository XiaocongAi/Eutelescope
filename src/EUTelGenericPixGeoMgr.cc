//STL
#include <map>
#include <iostream>
#include <stdexcept>

//EUTelescope
#include "EUTelGenericPixGeoMgr.h"

// MARLIN
#include "marlin/Global.h"
#include "marlin/VerbosityLevels.h"

//Geometry implementations
#include "Mimosa26GeoDescr.h"
#include "FEI4Single.h"
#include "FEI4Double.h"
#include "FEI4FourChip.h"

//ROOT includes
#include "TGeoBBox.h"
#include "TObjArray.h"

using namespace eutelescope;
using namespace geo;

//Default constructor
EUTelGenericPixGeoMgr::EUTelGenericPixGeoMgr() {}

//Destructor
EUTelGenericPixGeoMgr::~EUTelGenericPixGeoMgr() 
{
	//Clean up by deleting all EUTelGenericPixGeoDescr 
	std::map<std::string, EUTelGenericPixGeoDescr*>::iterator it;
	for(it = _pixelDescriptions.begin(); it != _pixelDescriptions.end(); ++it )
	{
		streamlog_out( MESSAGE3 ) << "Deleting " << (*it).first << std::endl;
		delete (*it).second;
	}
}


void EUTelGenericPixGeoMgr::addPlane(int planeID, std::string geoName, std::string planeVolume)
{
	EUTelGenericPixGeoDescr* pixgeodescrptr = NULL;
	std::map<std::string, EUTelGenericPixGeoDescr*>::iterator it;
	
	//Check if the geoemtry is loaded already (e.g. for some other plane)
	it = _pixelDescriptions.find(geoName);	
	if( it != _pixelDescriptions.end() )
	{
		//if it is, use it!
		streamlog_out( MESSAGE3 )  << "Found " << geoName << ", using it" << std::endl;
		pixgeodescrptr = (*it).second;
		streamlog_out( MESSAGE3 ) << "Inserting " << planeID << " into map" << std::endl;
		_geoDescriptions.insert( std::make_pair(planeID, pixgeodescrptr) );
	}	
	//Otherwise we load it
	else
	{
		streamlog_out( MESSAGE3 ) << "Didnt find " << geoName << " yet, thus creating" << std::endl;

		//MIMOSA26
		if(geoName == "Mimosa26.so")
		{
			pixgeodescrptr = new Mimosa26GeoDescr();
			streamlog_out( MESSAGE3 ) << "Inserting " << planeID << " into map" << std::endl;
			_geoDescriptions.insert( std::make_pair(planeID, pixgeodescrptr) );
			_pixelDescriptions.insert ( std::make_pair(geoName, pixgeodescrptr) );
		}		
		//FE-I4 Single
		else if(geoName == "FEI4Single.so")
		{
			pixgeodescrptr = new FEI4Single();
			streamlog_out( MESSAGE3 ) << "Inserting " << planeID << " into map" << std::endl;
			_geoDescriptions.insert( std::make_pair(planeID, pixgeodescrptr) );
			_pixelDescriptions.insert ( std::make_pair(geoName, pixgeodescrptr) );
		}
		//FE-I4 Double
		else if(geoName == "FEI4Double.so")
		{
			pixgeodescrptr = new FEI4Double();
			streamlog_out( MESSAGE3 ) << "Inserting " << planeID << " into map" << std::endl;
			_geoDescriptions.insert( std::make_pair(planeID, pixgeodescrptr) );
			_pixelDescriptions.insert ( std::make_pair(geoName, pixgeodescrptr) );
		}
		//FE-I4 FourChip
		else if(geoName == "FEI4FourChip.so")
		{
			pixgeodescrptr = new FEI4FourChip();
			streamlog_out( MESSAGE3 ) << "Inserting " << planeID << " into map" << std::endl;
			_geoDescriptions.insert( std::make_pair(planeID, pixgeodescrptr) );
			_pixelDescriptions.insert ( std::make_pair(geoName, pixgeodescrptr) );
		}
		//Unknown, TERMINATE!
		else
		{
			streamlog_out( ERROR3 ) << "EUTelGenericPixGeoMgr: Unknown geometry library: " << geoName << " Terminating!" << std::endl;
			//While runtime errors are not processed in Marlin, they will cause a controlled crash!
			throw std::runtime_error("Unknown geometry library!");
		}
	}

	streamlog_out( MESSAGE3 )  << "Adding plane: " << planeID << " with geoLibName: " << geoName << " in volume " << planeVolume << std::endl;
	
	//Call the factory method to actually load the geoemtry!
	pixgeodescrptr->createRootDescr(planeVolume);
}


EUTelGenericPixGeoDescr* EUTelGenericPixGeoMgr::getPixGeoDescr(int planeID)
{
	return static_cast< EUTelGenericPixGeoDescr* > ( _geoDescriptions.find(planeID)->second );
}
