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
#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <math.h>  
#include <LeMonADE/updater/AbstractUpdater.h>

/**
 * @class UpdaterReadCrosslinkConnections
 */

template <class IngredientsType>
class UpdaterReadCrosslinkConnections: public AbstractUpdater
{
  
  
public:
  
      UpdaterReadCrosslinkConnections(IngredientsType& ing_, const std::string input_, const double stepwidth_,const double minConversion_):ing(ing_), initialIng(ing_),input(input_), stepwidth(stepwidth_),minConversion(minConversion_){};
      virtual void initialize();
      virtual bool execute();
      virtual void cleanup(){};
  
private:

      //! container storing system information about monomers 
      IngredientsType& ing;
      
      //! container storing system information about monomers 
      IngredientsType initialIng;
      
      //! input of the cross link positions and so on..
      const std::string input;
      
      //! 
      const double stepwidth;
      
      //! minimal conversion to be read in before further analyse 
      double minConversion;
      
      //! 
      double conversion;
      
      //!
      uint32_t ReadNLines;
      
    //   //! read the next line 
    //   void readline();
      
      //! 
      uint32_t NMaxConnection;
      
      //!
      uint32_t NOfConnections;
      
      //! stream reading input file
      std::ifstream stream;
      
      //!
      uint32_t NConnections;
      
      //!
      uint32_t NMonomerPerChain;
    void ConnectCrossLinkToChain(uint32_t MonID, uint32_t chainID )
    {
      uint32_t ChainMonomer1( (chainID -1) *NMonomerPerChain  );
      uint32_t ChainMonomer2( (chainID   ) *NMonomerPerChain-1);
//       std::cout << "ConnectCrossLinkToChain: " << ChainMonomer1 <<" " << ChainMonomer2 << " " << chainID <<std::endl;
      if( !ing.getMolecules().areConnected(MonID,ChainMonomer1) 
	&& ing.getMolecules().getNumLinks(ChainMonomer1) < ing.getMolecules()[ChainMonomer1].getNumMaxLinks() 
	)
	ing.modifyMolecules().connect(MonID, ChainMonomer1 );
      else if ( !ing.getMolecules().areConnected(MonID,ChainMonomer2) 
	      && ing.getMolecules().getNumLinks(ChainMonomer2) < ing.getMolecules()[ChainMonomer2].getNumMaxLinks()
	      )
	ing.modifyMolecules().connect(MonID, ChainMonomer2 );
      else 
      {
	std::stringstream errormessage;
	for ( uint32_t i = 0 ; i < ing.getMolecules().getNumLinks(ChainMonomer1); i++ )
	  std::cout << "NeighborID  " <<  ing.getMolecules().getNeighborIdx(ChainMonomer1,i) << " for ChainMonomer1 " << ChainMonomer1  <<std::endl;
	for ( uint32_t i = 0 ; i < ing.getMolecules().getNumLinks(ChainMonomer2); i++ )
	  std::cout << "NeighborID  " <<  ing.getMolecules().getNeighborIdx(ChainMonomer2,i) << " for ChainMonomer2 " << ChainMonomer2  <<std::endl;
	errormessage << "Wrong monomer ID 1 for connection of cross link " << MonID << " with chain  " << ChainMonomer1 <<"\n";  
	errormessage << "Wrong monomer ID 2 for connection of cross link " << MonID << " with chain  " << ChainMonomer2 <<"\n";
	errormessage << "for chain ID " << chainID<<"\n";
	throw std::runtime_error(errormessage.str());
      }
    } 
};
template <class IngredientsType>
void UpdaterReadCrosslinkConnections<IngredientsType>::initialize()
{
	NMaxConnection=0;
	for (uint32_t i =0; i < ing.getMolecules().size();i++)
	{
	  if(ing.getMolecules()[i].isReactive())
	    NMaxConnection+=(ing.getMolecules()[i].getNumMaxLinks()-ing.getMolecules().getNumLinks(i));
	}
	std::cout <<"Number of reactive sites is : " << NMaxConnection << std::endl;
	if ( (NMaxConnection % 2) == 1 ) throw std::runtime_error("Wrong number of maximum connections!");
	NMaxConnection/=2;
	std::cout <<"Number of maximum connection: " << NMaxConnection << std::endl;
	
// 	stream.open(input);
// 	if(stream.fail()) throw std::runtime_error(std::string("error opening input file ")+input+std::string("\n"));
// 	bool findStartofData(false);
// 	//set head to beginning of the data block 
// 	while (!findStartofData)
// 	{
// 	  std::string line;
// 	  getline(stream,line);
// 	  if ( line.empty() )
// 	    findStartofData=true;
// 	  else if ( line.at(0) == '#') //go to next line 
// 	    continue;
// 	  else throw std::runtime_error("Wrong input format!");
// 	}
	NMonomerPerChain=ing.getChainLength();
	conversion=minConversion-stepwidth/100.;
	
}

template <class IngredientsType>
bool UpdaterReadCrosslinkConnections<IngredientsType>::execute()
{
  
  ing=initialIng;
  stream.open(input);
  if(stream.fail()) throw std::runtime_error(std::string("error opening input file ")+input+std::string("\n"));
  bool findStartofData(false);
  //set head to beginning of the data block 
  while (!findStartofData)
  {
    std::string line;
    getline(stream,line);
    if ( line.empty() )
      findStartofData=true;
    else if ( line.at(0) == '#') //go to next line 
      continue;
    else throw std::runtime_error("Wrong input format!");
  }
  conversion+=stepwidth/100.;
//   uint32_t ReadNLines=(floor(NMaxConnection*conversion)-NConnections);
  uint32_t ReadNLines=(floor(NMaxConnection*conversion));
  uint32_t NewConnections(0);
  std::cout <<  "Start reading "   <<ReadNLines <<" number of lines " <<std::endl; 
  while (NewConnections <= ReadNLines && stream.good() )
  { 
      std::string line;
      getline(stream, line);
      if (line.empty()) break;
      std::stringstream ss;
      uint32_t Time, ChainID, MonID1, P1X,P1Y,P1Z, MonID2, P2X, P2Y, P2Z;
      ss << line;
      ss >> Time >> ChainID >> MonID1 >> P1X >> P1Y >>P1Z >>MonID2 >> P2X >> P2Y >>P2Z ; 
      ing.modifyMolecules().setAge(Time);
      ing.modifyMolecules()[MonID1].modifyVector3D().setAllCoordinates(P1X,P1Y,P1Z);
      ConnectCrossLinkToChain(MonID1,ChainID);
      if(MonID2 >0)
	ing.modifyMolecules()[MonID2].modifyVector3D().setAllCoordinates(P2X,P2Y,P2Z);
      NewConnections++;
  }
  NConnections+=NewConnections;
  std::cout << "Read and add " << NewConnections << " connections to "
	    << NConnections << "/" << NMaxConnection  << " to the system at time " 
	    << ing.getMolecules().getAge() << std::endl;   
  ing.synchronize();
  if(stream.eof()) {stream.close();return false ; }
  else {stream.close();return true;}
}


