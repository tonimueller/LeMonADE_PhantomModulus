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

#ifndef LEMONADE_PM_ANALYZER_ANALYZERWRITEBFMFILESUBGROUPTM_H
#define LEMONADE_PM_ANALYZER_ANALYZERWRITEBFMFILESUBGROUPTM_H

#include <set>
#include <map>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <ctime>

#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include <LeMonADE/utility/DepthIterator.h>
#include <LeMonADE/utility/MonomerGroup.h>


/***********************************************************************/
/**
 * @file
 *
 * @class AnalyzerWriteBfmFileSubGroupTendomer
 *
 * @brief Analyzer writing the configurations into a given bfm-file.
 *
 * @details The output is appended to the file, if the file already exists.
 * If it does not exist, a new file is created and the header information is written
 * at the beginning
 *
 * @tparam IngredientsType Ingredients class storing all system information( e.g. monomers, bonds, etc).
 *
 * @todo rename to WriteBFMFile or similar.
 **/
template <class IngredientsType>
class AnalyzerWriteBfmFileSubGroupTendomer: public AnalyzerWriteBfmFile<IngredientsType>
{
public:
	//! Standard Constructor. Default: all data will be appended to the file
	AnalyzerWriteBfmFileSubGroupTendomer(const std::string& filename, const IngredientsType& ing,int writeType=AnalyzerWriteBfmFile<IngredientsType>::APPEND, std::string activeComponent_="active_00001.dat");

	//!Standard destructor closing the file stream.
	virtual ~AnalyzerWriteBfmFileSubGroupTendomer(){}


	//! Triggers writing of the Writes that need to be written for every step (i.e. excluding header)
	virtual bool execute();

	//! Initializes the writing. Opens the file and check for IO.
	virtual void initialize();


	void setAlwaysReindex(bool value){alwaysReindexGroups=value;}
private:

	//! Storage for data that are processed to file (mostly Ingredients).
	IngredientsType ingredients;

	const IngredientsType& ingredientsAllData;
	MonomerGroup<typename IngredientsType::molecules_type> subgroup;
	
    //! file with labels to be sol-0,gel-1,active-2
    std::string activeComponent;

	//!flag used in execute()
	bool alwaysReindexGroups;
};

/***********************************************************************/
//constructor
/***********************************************************************/
/**
 * @details It initialize all internal values and passes all information to the corresponding classes.
 * It register all known commands from the Feature.
 *
 *
 * @param filename Name of the file to write-out
 * @param ing Class holding all information of the system (mainly Ingredients )
 * @param writeType ENUM-type BFM_WRITE_TYPE to specify the write-out
 * @param activeComponent_  encodes the active components of the system 
 * 
 * @todo Did I understand that correctly that the first !mcs is read not the last in the file?
 */
template <class IngredientsType>
AnalyzerWriteBfmFileSubGroupTendomer<IngredientsType>::
AnalyzerWriteBfmFileSubGroupTendomer(const std::string& filename, const IngredientsType& ing, int writeType,std::string activeComponent_)
:AnalyzerWriteBfmFile<IngredientsType>(filename,ing,writeType)
,ingredientsAllData(ing), activeComponent(activeComponent_),subgroup(ing.getMolecules())
,alwaysReindexGroups(false)
{

}


/***********************************************************************/
//bool execute
/***********************************************************************/
/**
 * @details Write-out of the next configuration of all Feature.
 *
 * @return True if everthing is alrigth. False if something goes wrong.
 */
template <class IngredientsType>
bool AnalyzerWriteBfmFileSubGroupTendomer<IngredientsType>::execute()
{
	//copies all ingredients information again
	if(alwaysReindexGroups==true)
		initialize();

	//copies only the monomers
	else
		ingredients.modifyMolecules()=subgroup.copyGroup();

	AnalyzerWriteBfmFile<IngredientsType>::execute();

}


/***********************************************************************/
//bool execute
/***********************************************************************/
/**
 * @details Handling the creation and openiong of file.
 * * If APPEND is specified and file exists all output is append to the existing file.
 * * If APPEND is specified and file don't exists, a new file is created and the header is written.
 * * If NEWFILE is specified and file don't exists don't exists, a new file is created and the header is written.
 * * If NEWFILE is specified and file exists, an exception is thrown.
 *
 * @throw <std::runtime_error> IO-error or Writing is not allowed or if unknown ENUM-type is used.
 */
template <class IngredientsType>
void AnalyzerWriteBfmFileSubGroupTendomer<IngredientsType>::initialize()
{
    	//make a local copy of ingredients.
	ingredients=ingredientsAllData;  // implicit calls synchronize
 // activeObjects are first all cross links and afterwards all chain: 0-sol, 1-gel/pending, 2-active 
        std::ifstream in(activeComponent);
        std::string line;

		auto nCrossLinks(ingredients.getNumCrossLinkers());
		auto nChains(ingredients.getNumTendomers());
		auto nSegments(ingredients.getNumMonomersPerChain()*2);
		auto nChainMonomers(nChains*nSegments);
		auto objectID(0);
		auto nActiveCrossLinks(0);
		auto nActiveTendomers(0);
        std::cout << "nCrossLinks   =" << nCrossLinks << std::endl;
        std::cout << "nChains       =" << nChains << std::endl;
        std::cout << "nSegments     =" << nSegments << std::endl;
        std::cout << "nChainMonomers=" << nChainMonomers << std::endl;
		std::cout <<"Filestart:\n" ;
        while (in.good() &&  ! in.eof() ){
            getline(in,line);
            std::stringstream ss;
            uint32_t activeTag;
            ss << line;
			// std::cout << line <<"\n";
			if (! line.empty()) {
				ss >>activeTag;
				auto tmpAttribute(0);
				if(objectID < nCrossLinks ){
					auto IdX(nChainMonomers+objectID);
					#ifdef DEBUG
						std::cout <<  objectID << " "  << activeTag <<  " " << IdX <<  " " << IdX << std::endl;
					#endif
					if (activeTag == 2 ){
						tmpAttribute=1;	
						nActiveCrossLinks++;
					}
					ingredients.modifyMolecules()[nChainMonomers+objectID].setAttributeTag(tmpAttribute);
				}else{
					auto IdCStart((objectID-nCrossLinks)*nSegments);
					auto IdCEnd( (objectID-nCrossLinks+1)*nSegments-1);
					#ifdef DEBUG
						std::cout << objectID << " "  << activeTag <<" " << IdCStart << " " << IdCEnd << std::endl;
					#endif
					if (activeTag == 2 ) {
						nActiveTendomers++;
						tmpAttribute=1;
					}
					for (auto i =IdCStart ; i <= IdCEnd; i++ )
						ingredients.modifyMolecules()[i].setAttributeTag(tmpAttribute);
				}
				objectID++;
			}
        }
		std::cout <<"Fileend:\n" ;
		ingredients.setNumTendomers           (nActiveTendomers);
		ingredients.setNumCrossLinkers        (nActiveCrossLinks);
		ingredients.setNumMonomersPerChain    (ingredients.getNumMonomersPerChain());
		ingredients.setNumLabelsPerTendomerArm(ingredients.getNumLabelsPerTendomerArm());
		ingredients.synchronize();
	
		// MonomerGroup<typename Ing::molecules_type> subgroup(ingredients.getMolecules());
		subgroup.clear();
		hasThisType<1> predicate;
		for(size_t n=0;n<ingredients.getMolecules().size();n++)
		{
			if(predicate(ingredients.getMolecules(),n)==true) subgroup.push_back(n);
		}
		ingredients.modifyMolecules()=subgroup.copyGroup();

		for(size_t n=0;n<ingredients.getMolecules().size();n++)
		{
			ingredients.modifyMolecules()[n].setReactive(false);
			ingredients.modifyMolecules()[n].setNumMaxLinks(0);
			if( n < nActiveTendomers*nSegments ){
				if (n%(nSegments/2) == 0 ) {
					ingredients.modifyMolecules()[n].setReactive(true);
					ingredients.modifyMolecules()[n].setNumMaxLinks(2);
				}
			}else {
				ingredients.modifyMolecules()[n].setReactive(true);
				ingredients.modifyMolecules()[n].setNumMaxLinks(4);
			}
		}		
		ingredients.synchronize();




	ingredients.setName(this->getFilename());

	// //find the desired subgroup of monomers
	// subgroup.clear();
	// for(size_t n=0;n<ingredientsAllData.getMolecules().size();n++)
	// {
	// 	if(predicate(ingredientsAllData.getMolecules(),n)==true) subgroup.push_back(n);
	// }


	std::cout << "size of selected subgroup:" << subgroup.size() << std::endl;

	ingredients.modifyMolecules()=subgroup.copyGroup();
	ingredients.clearCompressedOutputIndices(); //necessary because the indices are now different

	std::cout << "sizeMolecules:" << this->ingredients.getMolecules().size() << std::endl;
	AnalyzerWriteBfmFile<IngredientsType>::initialize();

}




#endif
