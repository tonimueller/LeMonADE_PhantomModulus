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

#ifndef LEMONADE_PM_UPDATER_MOVES_MOVENONLINEARFORCEEQUILIBRIUM_H
#define LEMONADE_PM_UPDATER_MOVES_MOVENONLINEARFORCEEQUILIBRIUM_H
#include <limits>
#include <fstream>
#include <LeMonADE_PM/updater/moves/MoveNonLinearForceEquilibriumBase.h>
#include <LeMonADE/utility/DistanceCalculation.h>
#include <LeMonADE_PM/utility/neighborX.h>

/*****************************************************************************/
/**
 * @file
 *
 * @class MoveNonLinearForceEquilibrium
 *
 * @brief Standard local bfm-move on simple cubic lattice for the scBFM.
 *
 * @details The class is a specialization of MoveLocalBase using the (CRTP) to avoid virtual functions.
 * Here implemented for networks made of monodisperse chains.
 **/
/*****************************************************************************/

class MoveNonLinearForceEquilibrium:public MoveForceEquilibriumBase<MoveNonLinearForceEquilibrium>{
public:
    MoveNonLinearForceEquilibrium(std::filename input_, double relaxationChain_):input(input_),relaxationChain(relaxationChain_),bondlength(2.68){};

    // overload initialise function to be able to set the moves index and direction if neccessary
    template <class IngredientsType> void init(const IngredientsType& ing);
    template <class IngredientsType> void init(const IngredientsType& ing, uint32_t index);
    template <class IngredientsType> void init(const IngredientsType& ing, uint32_t index, VectorDouble3 dir );

    template <class IngredientsType> bool check(IngredientsType& ing);
    template< class IngredientsType> void apply(IngredientsType& ing);

private:
    //average square bond length 
    const double bondlength;
    //input filename for the force extension curve 
    std::string input;
    
    //minimum force in the file 
    double min_force:
    //maximum force in the file 
    double max_force;
    //force steps 
    double force_steps;
    
    //min_extension in the file 
    double min_extension;
    //maximum extensions in the file 
    double max_extension;
    //extension steps 
    double extension_steps;
    
    //force extension mapping
    std::map<double, double> extension_force
    //extension force mapping 
    std::vector<double> force_extension;

    //an equivalent chain which relaxes the cross link
    double relaxationChain;



    //read file to create a lookup table of the force extension curve 
    void createTable();

    //Gaussina force extension relation 
    //f=R*3/(N*b^2)
    VectorDouble3 FE(VectorDouble3 extensionVector, double nSegs){
        return extensionVector*3./(std::sqrt(nSegs)*bondlength*bondlength);
    }
    //Gaussian extension force relation 
    //R=-f/3*N*b^2
    VectorDouble3 EF(VectorDouble3 force, double nSegs){
        return force/(-3.)*std::sqrt(nSegs)*bondlength*bondlength;
    }

    //calculate the shift for the cross link
    template< class IngredientsType >
    VectorDouble3 CalculateShift(IngredientsType& ing ){
        std::vector<neighborX> Neighbors(ing.getCrossLinkNeighborIDs(this->getIndex()) );
        VectorDouble3 force(0.,0.,0.);
        VectorDouble3 shift(0.,0.,0.);
        double avNSegments(0.);
        if (Neighbors.size() > 0) {
            VectorDouble3 Position(ing.getMolecules()[this->getIndex()].getVector3D());      
                // std::cout << "CorssLinkPos=" << Position << std::endl;
            for (size_t i = 0; i < Neighbors.size(); i++){
                VectorDouble3 vec(Position-ing.getMolecules()[Neighbors[i].ID].getVector3D()-Neighbors[i].jump);
                // avNSegments+=1./Neighbors[i].segDistance;
                // force+=FE(vec,Neighbors[i].segDistance);
                force+=force_extension[vec.getLength()]*vec.normalize();
            }
            // force/=(1.*Neighbors.size());  
            shift=EF(force,relaxationChain);
        }
        return shift;
    };


};
/////////////////////////////////////////////////////////////////////////////
/////////// implementation of the members ///////////////////////////////////
void MoveNonLinearForceEquilibrium::createTable(){
    std::ifstream in(input);
    uint32_t counter(0);
    while(in.good()){
        std::string line;
        getline(stream, line);
        //ignore comments and blank lines 
        if (line.empty())
            findStartofData = true;
        else if (line.at(0) == '#') //go to next line
            continue;
        //read data 
        double force, extension;
        std::stringstream ss << line;
        ss>>force >> extension; 
        if(min_force > force ) min_force=force; 
        if(max_force < force ) max_force=force; 
        if(min_extension > extension ) min_extension=extension; 
        if(min_extension < extension ) max_extension=extension; 
        extension_force.insert(extension_force.end(),std::pair<double,double>(force, extension));
        if(counter==1){
            force_steps=max_force-min_force;
        }else if(counter>1)
        counter++;
    }
    in.close();
     //make lookup for the extension force relation 
     //make a entry from 0 to max_extension in steps of 1 
    for ( auto r=0;r<static_cast<uint32_t>(max_extension); r++   ){
        if(r==0)
            extension_force.push_back(0.);
        else{
            auto it_last=extension_force.begin()
            for (auto it=extension_force.begin(); it !=extension_force.end();it++){
                if(it->second > r){
                    //at the force linear interpolated in between the two current forces
                    auto deltaForce(it->first-it_last->first);
                    auto deltaExtension(it->second-it_last->second);
                    auto factor( (static_cast<double>(r)-it_last->second)/deltaExtension );
                    //at interpolated force
                    force_extension.push_back(it_last->first+deltaForce*factor);
                    break;
                }
                it_last=it;
            }

        }
    }
}


/*****************************************************************************/
/**
 * @brief Initialize the move.
 *
 * @details Resets the move probability to unity. Dice a new random direction and
 * Vertex (monomer) index inside the graph.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 **/
template <class IngredientsType>
void MoveNonLinearForceEquilibrium::init(const IngredientsType& ing)
{
    this->resetProbability();

    //draw index
    this->setIndex( (this->randomNumbers.r250_rand32()) %(ing.getMolecules().size()) );

    //calculate the shift of the cross link
    this->setShiftVector(CalculateShift(ing));
}

/*****************************************************************************/
/**
 * @brief Initialize the move with a given monomer index.
 *
 * @details Resets the move probability to unity. Dice a new random direction.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 * @param index index of the monomer to be connected
 **/
template <class IngredientsType>
void MoveNonLinearForceEquilibrium::init(const IngredientsType& ing, uint32_t index)
{
  this->resetProbability();

  //set index
  if( (index >= 0) && (index <= (ing.getMolecules().size()-1)) )
    this->setIndex( index );
  else
    throw std::runtime_error("MoveNonLinearForceEquilibrium::init(ing, index): index out of range!");

  //calculate the shift of the cross link
  this->setShiftVector(CalculateShift(ing));
  
}

/*****************************************************************************/
/**
 * @brief Initialize the move with a given monomer index.
 *
 * @details Resets the move probability to unity. Dice a new random direction.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 * @param index index of the monomer to be connected
 * @param bondpartner index of the monomer to connect to 
 **/
template <class IngredientsType>
void MoveNonLinearForceEquilibrium::init(const IngredientsType& ing, uint32_t index, VectorDouble3 dir )
{
  this->resetProbability();

  //set index
  if( (index >= 0) && (index <= (ing.getMolecules().size()-1)) )
    this->setIndex( index );
  else
    throw std::runtime_error("MoveNonLinearForceEquilibrium::init(ing, index, bondpartner): index out of range!");

  //calculate the shift of the cross link
  this->setShiftVector(dir);
}

/*****************************************************************************/
/**
 * @brief Check if the move is accepted by the system.
 *
 * @details This function delegates the checking to the Feature.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 * @return True if move is valid. False, otherwise.
 **/
template <class IngredientsType>
bool MoveNonLinearForceEquilibrium::check(IngredientsType& ing)
{
  //send the move to the Features to be checked
  return ing.checkMove(ing,*this);
}

/*****************************************************************************/
/**
 * @brief Apply the move to the system , e.g. add the displacement to Vertex (monomer) position.
 *
 * @details As first step: all Feature should apply the move using applyMove().\n
 * Second: Modify the positions etc. of the Vertex etc.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 **/
template< class IngredientsType>
void MoveNonLinearForceEquilibrium::apply(IngredientsType& ing)
{

	//move must FIRST be applied to the features
	ing.applyMove(ing,*this);
	//THEN the position can be modified
	ing.modifyMolecules()[this->getIndex()]+=this->getShiftVector();
}

#endif /*MOVEFORCEEQUILIBRIUM_H*/
