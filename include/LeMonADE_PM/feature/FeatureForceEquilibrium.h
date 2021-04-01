/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2015 by
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers (see AUTHORS)
    ooo                        |
----------------------------------------------------------------------------------

This file is part of LeMonADE.

LeMonADE is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

LeMonADE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with LeMonADE.  If not, see <http://www.gnu.org/licenses/>.

--------------------------------------------------------------------------------*/


#ifndef LEMONADE_PM_FEATURE_FEATUREFORCEEQUILIBRIUM_H
#define LEMONADE_PM_FEATURE_FEATUREFORCEEQUILIBRIUM_H

#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/feature/FeatureConnectionSc.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE_PM/updater/moves/MoveForceEquilibrium.h>


/*****************************************************************************/
/**
 * @file
 * @date   2019/03/05
 * @author Toni
 *
 * @class FeatureForceEquilibrium
 * 
 * */

///////////////////////////////////////////////////////////////////////////////
//DEFINITION OF THE CLASS TEMPLATE   	                                ///////
//Implementation of the members below					///////
///////////////////////////////////////////////////////////////////////////////

class FeatureForceEquilibrium : public Feature {
  
public:
	FeatureForceEquilibrium():BuildTable(false){};
	
  	//! This Feature requires a monomer_extensions.
	typedef LOKI_TYPELIST_1(MonomerReactivity) monomer_extensions;
  
	//! Export the relevant functionality for reading bfm-files to the responsible reader object
	template<class IngredientsType>
	void exportRead(FileImport<IngredientsType>& fileReader);

	//! Export the relevant functionality for writing bfm-files to the responsible writer object
	template<class IngredientsType>
	void exportWrite(AnalyzerWriteBfmFile<IngredientsType>& fileWriter) const;
	
	
	//! check bas connect move - always true 
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients, const MoveBase& move) const { return true;};
	//! check bas connect move - always true 
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients, const MoveForceEquilibrium& move) const;
	//! Synchronize with system
	template<class IngredientsType>
	void synchronize(IngredientsType& ingredients)
	{
	  std::cout << "FeatureForceEquilibrium::synchronize() -> fillTables() "<<std::endl;
	  fillTables(ingredients);
	  std::cout << "done."<<std::endl;
	};
	//!getter function for the neighboring crosslinks
	std::vector<uint32_t> getCrossLinkNeighborIDs(uint32_t CrossLinkID) const 
	{
	  if (CrossLinkNeighbors.find(CrossLinkID) == CrossLinkNeighbors.end())
	  {
	    std::stringstream errormessage;
	    errormessage << "FeatureForceEquilibrium::getCrossLinkNeighborIDs Cross Link ID " << CrossLinkID <<" does not exist.";
	    throw std::runtime_error(errormessage.str());
	  }
	  
	  return CrossLinkNeighbors.at(CrossLinkID);
	};
	//! returns the chain Id for a pair of cross links 
	uint32_t getChainIDByPair(uint32_t MonID1, uint32_t MonID2) const 
	{
	  std::pair<uint32_t, uint32_t> CrosslinkPair=std::make_pair(std::min( MonID1,MonID2),std::max( MonID1,MonID2));
	  if ( CrossLinkPairChainTable.find(CrosslinkPair) == CrossLinkPairChainTable.end())
	  {
	    std::stringstream errormessage;
	    errormessage << "FeatureForceEquilibrium::getChainIDByPair Cross link pair " << CrosslinkPair.first << " and " << CrosslinkPair.second  <<" does not exist.";
	    throw std::runtime_error(errormessage.str());	    
	  }
	  return CrossLinkPairChainTable.at( CrosslinkPair );
	  
	}
private:
  
  //!
  bool BuildTable;

  //! convinience function to fill all tables 
  template<class IngredientsType>
  void fillTables(IngredientsType& ingredients);
  
  //!key CrossLink ID; value neighboring cross links
  std::map<uint32_t,std::vector<uint32_t> > CrossLinkNeighbors;
  
  //! key: pair of cross link IDs, value: chain ID 
  std::map<std::pair<uint32_t,uint32_t>,uint32_t> CrossLinkPairChainTable;
  
  //!
  std::vector<uint32_t> MonomerChainID;
};

template<class IngredientsType>
void FeatureForceEquilibrium::fillTables(IngredientsType& ingredients)
{

  uint32_t ChainID(1); 
  std::cout << "Fill tables and scan "<<  ingredients.getMolecules().size() << " monomers."<<std::endl;
  for (size_t i =0 ; i < ingredients.getMolecules().size(); i++ )
  {
    MonomerChainID.push_back(0);
      if (  ingredients.getMolecules()[i].isReactive() 
	&& (ingredients.getMolecules()[i].getNumMaxLinks() == 2)
	)
      {
	MonomerChainID.back()=ChainID;
	if (ChainID < 10 )
	  std::cout << i << " "<< ingredients.getMolecules().getNumLinks(i)  << " "  << MonomerChainID.back()  <<std::endl;
	if ( (i % 2) == 1  )
	  ChainID++;
      }
  }
  std::cout << "FeatureForceEquilibrium::fillTables Found " << ChainID-1 << " number of chains." <<std::endl; 
  for (uint32_t i = 0 ;i < ingredients.getMolecules().size();i++)
  {
    if( ingredients.getMolecules()[i].isReactive()  && ingredients.getMolecules()[i].getNumMaxLinks() > 2 )
    {
      std::vector<uint32_t> NeighborIDs;
      for (size_t j = 0 ; j < ingredients.getMolecules().getNumLinks(i); j++)
      {
	uint32_t tail(i);
	uint32_t head(ingredients.getMolecules().getNeighborIdx(i,j));
	uint32_t FirstChainMonomer(head);
	bool FoundCrossLink(false);
	while( ingredients.getMolecules().getNumLinks(head) == 2   && !FoundCrossLink  )
	{
	  for (size_t k = 0 ; k < ingredients.getMolecules().getNumLinks(head); k++)
	  {
	    uint32_t NextMonomer( ingredients.getMolecules().getNeighborIdx(head,k));
	    if ( NextMonomer != tail ) 
	    {
	      tail=head;
	      head=NextMonomer; 
	      break;
	    }
	  }
	  if (ingredients.getMolecules()[head].getNumMaxLinks() > 2) 
	      FoundCrossLink=true;
	}
	if( ingredients.getMolecules().getNumLinks(head) > 1 )
	{
	  NeighborIDs.push_back(head);
// 	  std::cout << "Pair of chain monomers: " <<std::min( i,head)<< " " << std::max( i,head) <<std::endl;
	  CrossLinkPairChainTable[std::make_pair(std::min( i,head ),std::max( i,head ) )]=MonomerChainID.at(FirstChainMonomer);
	}
      }
      CrossLinkNeighbors[i]=NeighborIDs;
    }
  }
  
}


/**
 * @details The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type FileImport. The export of the Reads is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Read-In Commands:
 * * !reactivity
 *
 * @param fileReader File importer for the bfm-file
 * @param destination List of Feature to write-in from the read values.
 * @tparam IngredientsType Features used in the system. See Ingredients.
 **/
template<class IngredientsType>
void FeatureForceEquilibrium::exportRead(FileImport< IngredientsType >& fileReader)
{
  fileReader.registerRead("!reactivity",new ReadReactivity<IngredientsType>(fileReader.getDestination()));
}


/**
 * The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type AnalyzerWriteBfmFile. The export of the Writes is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Write-Out Commands:
 * * !reactivity
 *
 * @param fileWriter File writer for the bfm-file.
 */
template<class IngredientsType>
void FeatureForceEquilibrium::exportWrite(AnalyzerWriteBfmFile< IngredientsType >& fileWriter) const
{
  fileWriter.registerWrite("!reactivity",new WriteReactivity<IngredientsType>(fileWriter.getIngredients_()));
}

template<class IngredientsType>
bool FeatureForceEquilibrium::checkMove(const IngredientsType& ingredients, const MoveForceEquilibrium& move) const
{

//   if(ingredients.getMolecules().getNumLinks(move.getIndex()) < 2 ) return false ;
  return true; 
}

#endif /*LEMONADE_PM_FEATURE_FEATUREFORCEEQUILIBRIUM_H*/