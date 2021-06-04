/****************************************************************************** 
 * based on LeMonADE: https://github.com/LeMonADE-project/LeMonADE/
 * author: Toni Müller
 * email: mueller-toni@ipfdd.de
 * project: Phantom Modulus
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
 * @brief Analyzes the connection distribution of chains and cross linker for the whole system, 
 * the gel and the active part of the gel.
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
public:
	//! constructor
	AnalyzerIntramolecularReactions(const IngredientsType& ingredients_, std::string outputFilename_);

	//! destructor. does nothing
	virtual ~AnalyzerIntramolecularReactions(){}

	//! Initializes data structures. Called by TaskManager::initialize()
	virtual void initialize();

	//! Calculates the Rg2 for the current timestep. Called by TaskManager::execute()
	virtual bool execute();

	//! Writes the final results to file
	virtual void cleanup();
	
	//! name of output file for the average is outputFilePrefix_averages.dat  
	std::string outputFilename;

	//! save the current values in Rg2TimeSeriesX, etc., to disk
	void dumpData();

	//! calculates the distance between crosslinks and stores IDs, distance vector and chainID
	std::vector< std::vector<double> >  CalculateDistance();

	//! just collects the id and the position for the cross links 
	std::vector<std::vector<int> > CollectAveragePositions();
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
/**
 * @brief calculates the average distances between monomers and their distribution 
 * */
template< class IngredientsType >
void AnalyzerIntramolecularReactions<IngredientsType>::initialize(){}
/**
 * @details 
 * */
template< class IngredientsType >
bool AnalyzerIntramolecularReactions<IngredientsType>::execute()
{
  dumpData();
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
	std::vector< std::vector<int> > CrossLinkPositions=CollectAveragePositions() ;
	std::stringstream commentAveragePosition;
	commentAveragePosition<<"Created by AnalyzerIntramolecularReactions\n";
	commentAveragePosition<<"conversion nBonds nIntraBonds nBondsGel nIntraBondsGels \n";

	ResultFormattingTools::writeResultFile(
		outputFilename,
		ingredients,
		CrossLinkPositions,
		commentAveragePosition.str()
	);
			
	
}

#endif /*LEMONADE_PM_ANALYZER_ANALYZERINTRAMOLECULARREACTIONS_H*/


