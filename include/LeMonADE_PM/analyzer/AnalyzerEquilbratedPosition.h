/****************************************************************************** 
 * based on LeMonADE: https://github.com/LeMonADE-project/LeMonADE/
 * author: Toni Müller
 * email: mueller-toni@ipfdd.de
 * project: topological effects 
 *****************************************************************************/

#ifndef LEMONADE_PM_ANALYZER_ANALYZEREQUILIBRATEPOSITON_H
#define LEMONADE_PM_ANALYZER_ANALYZEREQUILIBRATEPOSITON_H

#include <string>

#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/analyzer/AbstractAnalyzer.h>
#include <LeMonADE/utility/ResultFormattingTools.h>
#include <LeMonADE/utility/MonomerGroup.h>
#include <LeMonADE/utility/DistanceCalculation.h>

/*************************************************************************
 * definition of AnalyzerEquilbratedPosition class
 * ***********************************************************************/

/**
 * @file
 *
 * @class AnalyzerEquilbratedPosition
 *
 * @brief Analyzer for evaluating ...
 *
 * @tparam IngredientsType Ingredients class storing all system information( e.g. monomers, bonds, etc).
 *
 * @details 
 */
template < class IngredientsType > class AnalyzerEquilbratedPosition : public AbstractAnalyzer
{

private:
	//! typedef for the underlying container holding the monomers
	typedef typename IngredientsType::molecules_type molecules_type;

	//! reference to the complete system
	const IngredientsType& ingredients;
public:
	//! constructor
	AnalyzerEquilbratedPosition(const IngredientsType& ingredients_, std::string outAvPosBasename_, std::string outDistBasename_);

	//! destructor. does nothing
	virtual ~AnalyzerEquilbratedPosition(){}

	//! Initializes data structures. Called by TaskManager::initialize()
	virtual void initialize();

	//! Calculates the Rg2 for the current timestep. Called by TaskManager::execute()
	virtual bool execute();

	//! Writes the final results to file
	virtual void cleanup();
	
	//! name of output file for the average is outputFilePrefix_averages.dat  
	std::string outAvPosBasename;
	
	//!name of the output file for the distribution is outputFilePrefix_timeseries.dat
	std::string  outDistBasename;

	// //! key: pair of cross link IDs, value: chain ID 
  	// std::map<std::pair<uint32_t,uint32_t>,uint32_t> CrossLinkPairChainTable;
	
	//! count up for each frame and give a unique number for the output files 
	uint32_t nExecutions;

	//! save the current values in Rg2TimeSeriesX, etc., to disk
	void dumpData();
	//! calculates the distance between crosslinks and stores IDs, distance vector and chainID
	std::vector< std::vector<double> >  CalculateDistance();
	//! just collects the id and the position for the cross links 
	std::vector<std::vector<double> > CollectAveragePositions();
	//! returns the chain Id for a pair of cross links 
	uint32_t getChainIDByPair(uint32_t MonID1, uint32_t MonID2) const ;
};

/*************************************************************************
 * implementation of memebers
 * ***********************************************************************/

/**
 * @param ing reference to the object holding all information of the system
 * */
template<class IngredientsType>
AnalyzerEquilbratedPosition<IngredientsType>::AnalyzerEquilbratedPosition(
	const IngredientsType& ingredients_, std::string outAvPosBasename_, std::string outDistBasename_)
:ingredients(ingredients_)
,outAvPosBasename(outAvPosBasename_)
,outDistBasename(outDistBasename_)
,nExecutions(0)
{}
////////////////////////////////////////////////////////////////////////////////
template< class IngredientsType >
std::vector< std::vector<double> >  AnalyzerEquilbratedPosition<IngredientsType>::CalculateDistance(){
	std::vector< std::vector<double> >  dist(7,std::vector<double>());
	auto crosslinkID(ingredients.getCrosslinkIDs());
	for (size_t i = 0 ; i < crosslinkID.size(); i++){
		auto IDx(crosslinkID[i]);
		std::vector<uint32_t> neighbors(ingredients.getCrossLinkNeighborIDs(IDx));
		for (size_t j=0; j < neighbors.size() ;j++){
			VectorDouble3 vec(ingredients.getMolecules()[IDx].getVector3D()-ing.getMolecules()[neighbors[j].ID].getVector3D() -neighbors[j].jump);
			dist[0].push_back(IDx);
			dist[1].push_back(neighbors[j]);
			dist[2].push_back(vec.getX());
			dist[3].push_back(vec.getY());
			dist[4].push_back(vec.getZ());
			dist[5].push_back(vec.getLength());
			// dist[6].push_back(getChainIDByPair(IDx,neighbors[j]));
			dist[6].push_back(neighbors[j].segDistance);
		}
	}
	return dist;
}
////////////////////////////////////////////////////////////////////////////////
template< class IngredientsType >
std::vector< std::vector<double> >  AnalyzerEquilbratedPosition<IngredientsType>::CollectAveragePositions(){
	std::vector<std::vector<double> > AveragePosition(4, std::vector<double>());
	auto crosslinkID(ingredients.getCrosslinkIDs());
	for (size_t i = 0 ; i < crosslinkID.size(); i++){
		auto IDx(crosslinkID[i]);
		AveragePosition[0].push_back(IDx);
		AveragePosition[1].push_back(ingredients.getMolecules()[IDx].getX() );
		AveragePosition[2].push_back(ingredients.getMolecules()[IDx].getY() );
		AveragePosition[3].push_back(ingredients.getMolecules()[IDx].getZ() );
	}
	return AveragePosition;
}
////////////////////////////////////////////////////////////////////////////////
template< class IngredientsType >
uint32_t  AnalyzerEquilbratedPosition<IngredientsType>::getChainIDByPair(uint32_t MonID1, uint32_t MonID2) const {
	std::pair<uint32_t, uint32_t> CrosslinkPair=std::make_pair(std::min( MonID1,MonID2),std::max( MonID1,MonID2));
	if ( CrossLinkPairChainTable.find(CrosslinkPair) == CrossLinkPairChainTable.end()){
		std::stringstream errormessage;
		errormessage << "AnalyzerEquilbratedPosition::getChainIDByPair Cross link pair " << CrosslinkPair.first << " and " << CrosslinkPair.second  <<" does not exist.";
		throw std::runtime_error(errormessage.str());	    
	}
	return CrossLinkPairChainTable.at( CrosslinkPair );
}
////////////////////////////////////////////////////////////////////////////////
/**
 * @brief calculates the average distances between monomers and their distribution 
 * */
template< class IngredientsType >
void AnalyzerEquilbratedPosition<IngredientsType>::initialize(){}
/**
 * @details 
 * */
template< class IngredientsType >
bool AnalyzerEquilbratedPosition<IngredientsType>::execute()
{
  dumpData();
  return true;
}
template<class IngredientsType>
void AnalyzerEquilbratedPosition<IngredientsType>::cleanup()
{
  dumpData();
}


template<class IngredientsType>
void AnalyzerEquilbratedPosition<IngredientsType>::dumpData()
{
  	double conversion, NReactedSites(0.0), NReactiveSites(0.0);
	auto crosslinkID(ingredients.getCrosslinkIDs());
	for (size_t i = 0 ; i < crosslinkID.size(); i++){
		auto IDx(crosslinkID[i]);			
		uint32_t NLinks(ingredients.getMolecules().getNumLinks(IDx));
		uint32_t nIrreversibleBonds=0;
		for (uint32_t n = 0 ; n < NLinks ;n++){
			uint32_t neighbor(ingredients.getMolecules().getNeighborIdx(IDx,n));
			if( ingredients.getMolecules()[neighbor].isReactive() )
				NReactedSites++;
			else
				nIrreversibleBonds++;
		}
		NReactiveSites+=(ingredients.getMolecules()[IDx].getNumMaxLinks()-nIrreversibleBonds);
	}
  
  conversion=NReactedSites/NReactiveSites;
  std::vector< std::vector<double> > CrossLinkPositions=CollectAveragePositions() ;
 
  //output for the equilibrated positions 
  std::stringstream commentAveragePosition;
  commentAveragePosition<<"Created by AnalyzerEquilbratedPosition\n";
  commentAveragePosition<<"ID's start at 0 \n";
  commentAveragePosition<<"ID equilibrated position\n";
  std::stringstream outAvPos;
  outAvPos<<   std::setw(6) << std::setfill('0') << nExecutions;
  outAvPos << "_" << outAvPosBasename;
  

  ResultFormattingTools::writeResultFile(
	  outAvPos.str(),
	  ingredients,
	  CrossLinkPositions,
	  commentAveragePosition.str()
  );
		
  // chain stretching distribution 
  std::vector< std::vector<double> >  dist=CalculateDistance();
  
  std::stringstream commentDistribution;
  commentDistribution<<"Created by AnalyzerEquilbratedPosition\n";
  commentDistribution<<"Monomer ID's start at 0 \n";
  commentDistribution<<"Chain ID's start at 1 \n";
  commentDistribution<<"ID1 ID2 vector length ChainID \n";
  std::stringstream outDist;
  outDist<<   std::setw(6) << std::setfill('0') << nExecutions;
  outDist << "_" << outDistBasename;

  ResultFormattingTools::writeResultFile(
	  outDist.str(),
	  ingredients,
	  dist,
	  commentDistribution.str()
  );
  
  nExecutions++;
}

#endif /*LEMONADE_PM_ANALYZER_ANALYZEREQUILIBRATEPOSITON_H*/


