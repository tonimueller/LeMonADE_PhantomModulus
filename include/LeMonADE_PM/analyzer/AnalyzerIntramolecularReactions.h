/****************************************************************************** 
 * based on LeMonADE: https://github.com/LeMonADE-project/LeMonADE/
 * author: Toni Müller
 * email: mueller-toni@ipfdd.de
 * project:Phantom Modulus
 *****************************************************************************/

#ifndef LEMONADE_PM_ANALYZER_ANALYZERINTRAMOLECULARREACTIONS_H
#define LEMONADE_PM_ANALYZER_ANALYZERINTRAMOLECULARREACTIONS_H

#include <string>

#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/analyzer/AbstractAnalyzer.h>
#include <LeMonADE/utility/ResultFormattingTools.h>
#include <LeMonADE/utility/MonomerGroup.h>
#include <LeMonADE/utility/DistanceCalculation.h>

/*************************************************************************
 * definition of AnalyzerIntramolecularReactions class
 * ***********************************************************************/

/**
 * @file
 *
 * @class AnalyzerIntramolecularReactions
 *
 * @brief Analyzer counting the intramolecular reactions
 *
 * @tparam IngredientsType Ingredients class storing all system information( e.g. monomers, bonds, etc).
 */
template < class IngredientsType > class AnalyzerIntramolecularReactions : public AbstractAnalyzer
{
private:
	//! typedef for the underlying container holding the monomers
	typedef typename IngredientsType::molecules_type molecules_type;
	//! reference to the complete system
	const IngredientsType& ingredients;
	//!cotnainer holding the information about conversion nBonds nIntraBonds nBondsGel nIntraBondsGels
	std::vector<std::vector<dobule> > data;
	//! name of output file for the average is outputFilePrefix_averages.dat  
	std::string outputFilename;
	//! save the collected data to a file 
	void dumpData();
	//! iterates over the structure and searches for the intramolecular reactions
	void getBiggestCluster();
public:
	//! constructor
	AnalyzerIntramolecularReactions(const IngredientsType& ingredients_, std::string outputFilename_);
	//! destructor. does nothing
	virtual ~AnalyzerIntramolecularReactions(){}
	//! nothing needs  to be initialized
	virtual void initialize();
	//! Runs over the structure and obtains the conversion nBonds nIntraBonds nBondsGel nIntraBondsGels.
	virtual bool execute();
	//! Writes the final results to file
	virtual void cleanup();	
};

/*************************************************************************
 * implementation of memebers
 * ***********************************************************************/

/**
 * @param ing reference to the object holding all information of the system
 * */
template<class IngredientsType>
AnalyzerIntramolecularReactions<IngredientsType>::AnalyzerIntramolecularReactions(
	const IngredientsType& ingredients_, std::string outputFilename_)
:ingredients(ingredients_)
,outputFilename(outputFilename_)
{}
////////////////////////////////////////////////////////////////////////////////
template< class IngredientsType >
void AnalyzerIntramolecularReactions<IngredientsType>::initialize(){}

template< class IngredientsType >
bool AnalyzerIntramolecularReactions<IngredientsType>::execute()
{
	
 	return true;
}
template<class IngredientsType>
void AnalyzerIntramolecularReactions<IngredientsType>::cleanup()
{
  	dumpData();
}
template<class IngredientsType>
void AnalyzerIntramolecularReactions<IngredientsType>::dumpData()
{  
	//output for the equilibrated positions 
	std::stringstream commentAveragePosition;
	commentAveragePosition<<"Created by AnalyzerIntramolecularReactions\n";
	commentAveragePosition<<"nTotalBonds="<<myIngredients.getNumCrossLinkers()*4<<"\n";
	commentAveragePosition<<"conversion nBonds nIntraBonds nBondsGel nIntraBondsGels \n";
	ResultFormattingTools::writeResultFile(
		outAvPos.str(),
		ingredients,
		data,
		commentAveragePosition.str()
	);
}


template <class IngredientsType>
void AnalyzerIntramolecularReactions<IngredientsType>::getBiggestCluster () {
	uint32_t nMolecules(0);
	std::vector<uint32_t> LargestCluster;
	const typename IngredientsType::molecules_type& getMolies = ingredients.getMolecules(); 
	//search largest cluster: use the attribute tag for that
	//at first store the initial attributes to reset them after the search
	//then color the whole graph. connected parts have the same color
	auto nMonomers(getMolies.size());
	std::vector<uint32_t> Tag(nMonomers,0);
	//color graph
	int freeColor(1);
	for (auto i=0 ; i < nMonomers;i++){
		std::vector<size_t> rememberBranch;
		rememberBranch.push_back(i);
		if(Tag[i] == 0) {
			while ( rememberBranch.size() > 0 ){
				auto ID(rememberBranch.back());
				auto color(Tag[ID]);
				rememberBranch.pop_back();
				if ( color == 0){
					Tag[ID]=freeColor;
					auto nLinks(getMolies.getNumLinks(ID));
					if(nLinks != 0 ){
						for(auto j=0; j < nLinks ; j++){
							auto Neighbor(getMolies.getNeighborIdx(ID,j));
							if( Tag[Neighbor] == 0 ) rememberBranch.push_back(Neighbor);
						}
					}
				}
			}
			freeColor++;
		}
	} 
	freeColor--;  
	nMolecules=freeColor;
	std::cout << "#Molecules=" << nMolecules <<std::endl;
	std::vector<std::vector<uint32_t> > ColoredGraphIDs(freeColor+1,std::vector<uint32_t>(0));
	for (auto i=0 ; i < nMonomers;i++){
		uint32_t atti(Tag[i]);
		if (atti==0 ){//this should never happen because all monomers should be covered above and colored 
			std::stringstream error_message;
			error_message << "AnalyzerActiveChainDensity::execute(): Found monomer which is not colored! ID is " << i << "\n";
			throw std::runtime_error(error_message.str());
		}
		ColoredGraphIDs[atti].push_back(i);
	}
	auto biggestClusterID(0);
	auto biggestClusterSize(ColoredGraphIDs[biggestClusterID].size());
	for (auto i=0; i < freeColor+1 ; i++){
		std::cout << "ClusterID  "  << i << " cluster size  " << ColoredGraphIDs[i].size() << std::endl; 
		if ( biggestClusterSize < ColoredGraphIDs[i].size() ){ 
			biggestClusterSize=ColoredGraphIDs[i].size();
			biggestClusterID=i;
		}
	}
	std::cout << "biggest ClusterID  "  << biggestClusterID << " biggestClusterSize " << biggestClusterSize << std::endl; 
	LargestCluster=ColoredGraphIDs[biggestClusterID];
	// return LargestCluster; 
}

#endif /*LEMONADE_PM_ANALYZER_ANALYZERINTRAMOLECULARREACTIONS_H*/


