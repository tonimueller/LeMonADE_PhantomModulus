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
	
	//! name of output files are outputFilePrefix_averages.dat and outputFilePrefix_timeseries.dat
	std::string outAvPosBasename;
	//!
	std::string  outDistBasename;
	//! key: pair of cross link IDs, value: chain ID 
  	std::map<std::pair<uint32_t,uint32_t>,uint32_t> CrossLinkPairChainTable;
  
	//! monomer to chain ID 
	std::vector<uint32_t> MonomerChainID;

	//! save the current values in Rg2TimeSeriesX, etc., to disk
	void dumpData();
	//! calculated the distances between pairs of cross links 
	std::vector< std::vector<double> >  CalculateDistance(){
	    std::vector< std::vector<double> >  dist(7,std::vector<double>());
	    for (size_t i = 0 ; i < ingredients.getMolecules().size(); i++){
	      	if ( ingredients.getMolecules()[i].isReactive() && ingredients.getMolecules()[i].getNumMaxLinks() > 2 ){ 
				std::vector<uint32_t> neighbors(ingredients.getCrossLinkNeighborIDs(i));
				for (size_t j=0; j < neighbors.size() ;j++){
					VectorDouble3 vec(LemonadeDistCalcs::MinImageVector( ingredients.getMolecules()[i].getVector3D(),ingredients.getMolecules()[neighbors[j]].getVector3D(),ingredients ));
					dist[0].push_back(i);
					dist[1].push_back(neighbors[j]);
					dist[2].push_back(vec.getX());
					dist[3].push_back(vec.getY());
					dist[4].push_back(vec.getZ());
					dist[5].push_back(vec.getLength());
					dist[6].push_back(ingredients.getChainIDByPair(i,neighbors[j]));
				}
			}
	    }
	    return dist;
	}
	//! just collects the id and the position for the cross links 
	std::vector<std::vector<double> > CollectAveragePositions(){
	  	std::vector<std::vector<double> > AveragePosition(4, std::vector<double>());
	    for (size_t i = 0 ; i < ingredients.getMolecules().size(); i++){
			if ( ingredients.getMolecules()[i].isReactive() && ingredients.getMolecules()[i].getNumMaxLinks() > 2 ){ 
					AveragePosition[0].push_back(i);
					AveragePosition[1].push_back(ingredients.getMolecules()[i].getX() );
					AveragePosition[2].push_back(ingredients.getMolecules()[i].getY() );
					AveragePosition[3].push_back(ingredients.getMolecules()[i].getZ() );
			}
	    }
	    return AveragePosition;
	}
	//! count up for each frame and give a unique number for the output files 
	uint32_t filenumber;
	//! returns the chain Id for a pair of cross links 
	uint32_t getChainIDByPair(uint32_t MonID1, uint32_t MonID2) const 
	{
		std::pair<uint32_t, uint32_t> CrosslinkPair=std::make_pair(std::min( MonID1,MonID2),std::max( MonID1,MonID2));
		if ( CrossLinkPairChainTable.find(CrosslinkPair) == CrossLinkPairChainTable.end())
		{
			std::stringstream errormessage;
			errormessage << "AnalyzerEquilbratedPosition::getChainIDByPair Cross link pair " << CrosslinkPair.first << " and " << CrosslinkPair.second  <<" does not exist.";
			throw std::runtime_error(errormessage.str());	    
		}
		return CrossLinkPairChainTable.at( CrosslinkPair );
	  
	}
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
,filenumber(0)
{}
/**
 * @details 
 * */
template< class IngredientsType >
void AnalyzerEquilbratedPosition<IngredientsType>::initialize(){
	// uint32_t ChainID(1); 
	// std::cout << "Fill tables and scan "<<  ingredients.getMolecules().size() << " monomers."<<std::endl;
	// for (size_t i =0 ; i < ingredients.getMolecules().size(); i++ ){
	// 		MonomerChainID.push_back(0);
	// 		if (  ingredients.getMolecules()[i].isReactive() && (ingredients.getMolecules()[i].getNumMaxLinks() == 2) ){
	// 			MonomerChainID.back()=ChainID;
	// 			if (ChainID < 10 )
	// 				std::cout << i << " "<< ingredients.getMolecules().getNumLinks(i)  << " "  << MonomerChainID.back()  <<std::endl;
	// 			if ( (i % 2) == 1  )
	// 				ChainID++;
	// 		}
	// }
	// std::cout << "AnalyzerEquilbratedPosition::initialize Found " << ChainID-1 << " number of chains." <<std::endl; 
	// for (uint32_t i = 0 ;i < ingredients.getMolecules().size();i++){
	// 		if( ingredients.getMolecules()[i].isReactive()  && ingredients.getMolecules()[i].getNumMaxLinks() > 2 ){
	// 			std::vector<uint32_t> NeighborIDs;
	// 			for (size_t j = 0 ; j < ingredients.getMolecules().getNumLinks(i); j++){
	// 				uint32_t tail(i);
	// 				uint32_t head(ingredients.getMolecules().getNeighborIdx(i,j));
	// 				uint32_t FirstChainMonomer(head);
	// 				bool FoundCrossLink(false);
	// 				while( ingredients.getMolecules().getNumLinks(head) == 2   && !FoundCrossLink  ){
	// 					for (size_t k = 0 ; k < ingredients.getMolecules().getNumLinks(head); k++){
	// 						uint32_t NextMonomer( ingredients.getMolecules().getNeighborIdx(head,k));
	// 						if ( NextMonomer != tail ) {
	// 							tail=head;
	// 							head=NextMonomer; 
	// 							break;
	// 						}
	// 					}
	// 					if (ingredients.getMolecules()[head].getNumMaxLinks() > 2) 
	// 						FoundCrossLink=true;
	// 				}
	// 				if( ingredients.getMolecules().getNumLinks(head) > 1 ){
	// 					NeighborIDs.push_back(head);
	// 				// 	  std::cout << "Pair of chain monomers: " <<std::min( i,head)<< " " << std::max( i,head) <<std::endl;
	// 					CrossLinkPairChainTable[std::make_pair(std::min( i,head ),std::max( i,head ) )]=MonomerChainID.at(FirstChainMonomer);
	// 				}
	// 			}
	// 		CrossLinkNeighbors[i]=NeighborIDs;
	// 		}
	// }
}
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
  
//   double conversion, NReactedSites(0.0), NReactiveSites(0.0);
//   for(size_t i = 0 ; i < ingredients.getMolecules().size(); i++ ){
// 		if ( ingredients.getMolecules()[i].isReactive() ){
// 			uint32_t NLinks(ingredients.getMolecules().getNumLinks(i));
// 			uint32_t nIrreversibleBonds=0;
// 			for (uint32_t n = 0 ; n < NLinks ;n++){
// 				uint32_t neighbor(ingredients.getMolecules().getNeighborIdx(i,n));
// 				if( ingredients.getMolecules()[neighbor].isReactive() )
// 					NReactedSites++;
// 				else
// 					nIrreversibleBonds++;
// 			}
// 			NReactiveSites+=(ingredients.getMolecules()[i].getNumMaxLinks()-nIrreversibleBonds);
// 		}
//   }
  
//   conversion=NReactedSites/NReactiveSites;
//   std::vector< std::vector<double> > CrossLinkPositions=CollectAveragePositions() ;
 
  //output for the equilibrated positions 
  std::stringstream commentAveragePosition;
  commentAveragePosition<<"Created by AnalyzerEquilbratedPosition\n";
  commentAveragePosition<<"ID's start at 0 \n";
//   commentAveragePosition<<"Conversion = " << conversion <<" \n";
  commentAveragePosition<<"ID equilibrated position\n";
  std::stringstream outAvPos;
  outAvPos<<   std::setw(6) << std::setfill('0') << filenumber;
  outAvPos << "_" << outAvPosBasename;
  

//   ResultFormattingTools::writeResultFile(
// 	  outAvPos.str(),
// 	  ingredients,
// 	  CrossLinkPositions,
// 	  commentAveragePosition.str()
//   );
		
  // chain stretching distribution 
  std::vector< std::vector<double> >  dist=CalculateDistance();
  
  std::stringstream commentDistribution;
  commentDistribution<<"Created by AnalyzerEquilbratedPosition\n";
  commentDistribution<<"Monomer ID's start at 0 \n";
  commentDistribution<<"Chain ID's start at 1 \n";
//   commentDistribution<<"Conversion = " << conversion <<" \n";
  commentDistribution<<"ID1 ID2 vector length ChainID \n";
  std::stringstream outDist;
  outDist<<   std::setw(6) << std::setfill('0') << filenumber;
  outDist << "_" << outDistBasename;

  ResultFormattingTools::writeResultFile(
	  outDist.str(),
	  ingredients,
	  dist,
	  commentDistribution.str()
  );
  
  filenumber++;
}

#endif /*LEMONADE_PM_ANALYZER_ANALYZEREQUILIBRATEPOSITON_H*/


