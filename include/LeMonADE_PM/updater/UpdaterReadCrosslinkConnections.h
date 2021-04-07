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
 * @brief reads in connections for a system made of linear chains and crosslinker 
 * 
 * @details 
 *   -the input file needs to have a format like:
 *     #Time, ChainID, MonID1, P1X, P1Y, P1Z, MonID2, P2X, P2Y, P2Z
 *   -the chains are before the crosslinks in the bfm file
 *   -the chain length must be at least 1 
 * 
 */

template <class IngredientsType>
class UpdaterReadCrosslinkConnections : public AbstractUpdater
{
public:
    UpdaterReadCrosslinkConnections(
        IngredientsType &ing_, 
        const std::string input_, 
        const double stepwidth_, 
        const double minConversion_): 
        ing(ing_), 
        initialIng(ing_), 
        input(input_), 
        stepwidth(stepwidth_), 
        minConversion(minConversion_),
        nExecutions(0){};
    virtual void initialize();
    virtual bool execute();
    virtual void cleanup(){};

private:
  //! container storing system information about monomers
  IngredientsType &ing;

  //! container storing system information about monomers
  IngredientsType initialIng;

  //! input of the cross link positions and so on..
  const std::string input;

  //!conversion step used in the execute function to increase the conversion 
  const double stepwidth;

  //! minimal conversion to be read in before further analyse
  double minConversion;

  //!number of maximum connections for a cross link
  uint32_t NMaxConnection;

  //!numbe of monomers per chain 
  uint32_t NMonomerPerChain;

  //!number of executions;
  uint32_t nExecutions;

  //! connects the cross link to the chain
  void ConnectCrossLinkToChain(uint32_t MonID, uint32_t chainID);
  
};


template <class IngredientsType>
void UpdaterReadCrosslinkConnections<IngredientsType>::ConnectCrossLinkToChain(uint32_t MonID, uint32_t chainID){
    uint32_t ChainMonomer1((chainID - 1) * NMonomerPerChain);
    uint32_t ChainMonomer2((chainID)*NMonomerPerChain - 1);
    if (!ing.getMolecules().areConnected(MonID, ChainMonomer1) && ing.getMolecules().getNumLinks(ChainMonomer1) < 2)
      ing.modifyMolecules().connect(MonID, ChainMonomer1);
    else if (!ing.getMolecules().areConnected(MonID, ChainMonomer2) && ing.getMolecules().getNumLinks(ChainMonomer2) < 2)
       ing.modifyMolecules().connect(MonID, ChainMonomer2);
    else{
        std::stringstream errormessage;
        for (uint32_t i = 0; i < ing.getMolecules().getNumLinks(ChainMonomer1); i++)
          std::cout << "NeighborID  " << ing.getMolecules().getNeighborIdx(ChainMonomer1, i) << " for ChainMonomer1 " << ChainMonomer1 << std::endl;
        for (uint32_t i = 0; i < ing.getMolecules().getNumLinks(ChainMonomer2); i++)
          std::cout << "NeighborID  " << ing.getMolecules().getNeighborIdx(ChainMonomer2, i) << " for ChainMonomer2 " << ChainMonomer2 << std::endl;
        errormessage << "Wrong monomer ID 1 for connection of cross link " << MonID << " with chain  " << ChainMonomer1 << "\n";
        errormessage << "Wrong monomer ID 2 for connection of cross link " << MonID << " with chain  " << ChainMonomer2 << "\n";
        errormessage << "for chain ID " << chainID << "\n";
        throw std::runtime_error(errormessage.str());
    }
}

/**
 * @brief read in connections up to the minimum conversion given 
 * */
template <class IngredientsType>
void UpdaterReadCrosslinkConnections<IngredientsType>::initialize(){
    //assume a stochiometric mixture
    NMaxConnection=ing.getFunctionality()*ing.getNumOfCrosslinks();
    NMonomerPerChain = ing.getNumOfMonomersPerChain();
    std::cout << "Number of maximum connection: " << NMaxConnection << std::endl;
}
/**
 * @brief read in connections up tp the next step 
 * @details Copies the initial ingredients to the current one and adds connections up to the current conversion. 
 * */
template <class IngredientsType>
bool UpdaterReadCrosslinkConnections<IngredientsType>::execute(){
    //reset the ingredients container to the inital one
    ing = initialIng;
    //open input file to the connection table 
    //! stream reading input file
    std::ifstream stream;
    stream.open(input);
    if (stream.fail())
      throw std::runtime_error(std::string("error opening input file ") + input + std::string("\n"));
    bool findStartofData(false);
    //set head to beginning of the data block
    while (!findStartofData){
        std::string line;
        getline(stream, line);
        if (line.empty())
            findStartofData = true;
        else if (line.at(0) == '#') //go to next line
            continue;
        else
            throw std::runtime_error("Wrong input format!");
    }
    //current conversion 
    auto conversion = minConversion + static_cast<double>(nExecutions) * stepwidth;
    std::cout << "Current conversion is " <<conversion <<std::endl;
    //read in number of lines to reach the current conversion 
    //(always start from the initial state and thus read in connectiosn from the beginning)
    uint32_t ReadNLines = (floor(NMaxConnection * conversion));
    //counter for the new connections
    uint32_t NewConnections(0);
    std::cout << "Start reading " << ReadNLines << " number of lines " << std::endl;
    while (NewConnections < ReadNLines && stream.good()){
        std::string line;
        getline(stream, line);
        if (line.empty())
            break;
        std::stringstream ss;
        uint32_t Time, ChainID, MonID1, P1X, P1Y, P1Z, MonID2, P2X, P2Y, P2Z;
        ss << line;
        ss >> Time >> ChainID >> MonID1 >> P1X >> P1Y >> P1Z >> MonID2 >> P2X >> P2Y >> P2Z;
        // std::cout <<  Time << " " << ChainID<< " "
        //           << MonID1<< " " << P1X << " " << P1Y << " " << P1Z << " " 
        //           << MonID2<< " " << P2X << " " << P2Y << " " << P2Z << "\n";
        ing.modifyMolecules().setAge(Time);
        ing.modifyMolecules()[MonID1].modifyVector3D().setAllCoordinates(P1X, P1Y, P1Z);
        ConnectCrossLinkToChain(MonID1, ChainID);
        if (MonID2 > 0)
            ing.modifyMolecules()[MonID2].modifyVector3D().setAllCoordinates(P2X, P2Y, P2Z);
        NewConnections++;
    }
    std::cout << "Read and add " << NewConnections << "/" << NMaxConnection 
              << " to the system at time "
              << ing.getMolecules().getAge() << std::endl;
    ing.synchronize();
    nExecutions++;
    std::cout << "UpdaterReadCrosslinkConnections::execute " << nExecutions << " times.\n\n";
    //close the filestream and return false if the file has ended and thus the updater has nothing more to do
    if (stream.eof()) {
        stream.close();
        return false;
    }else{
        stream.close();
        return true;
    }
}
