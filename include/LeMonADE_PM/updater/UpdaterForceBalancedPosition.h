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
#ifndef LEMONADE_PM_UPDATER_UPDATERFORCEBALANCEPOSITION_H
#define LEMONADE_PM_UPDATER_UPDATERFORCEBALANCEPOSITION_H


#include <LeMonADE/updater/AbstractUpdater.h>
#include <LeMonADE_PM/updater/moves/MoveForceEquilibrium.h>
#include <vector>
 /**
 * @class UpdaterForceBalancedPosition
 * @tparam IngredientsType
 */

 template <class IngredientsType>
class UpdaterForceBalancedPosition:public AbstractUpdater
{
public:
    
    UpdaterForceBalancedPosition(IngredientsType& ing_, double threshold_):ing(ing_),threshold(threshold_){};
    void initialize();
    bool execute();
    virtual void cleanup(){};
    uint32_t getNCrossLinks() const {return NCrossLinks;}
  
private:
    
    //!copy of the main container for the system informations 
    IngredientsType& ing;
    
    //! threshold for the certainty 
    double threshold;
    
    //! look up table for the cross link ids to monomer ids
    std::vector<uint32_t> CrossLinkIDs;
    
    //! number of cross links 
    uint32_t NCrossLinks;
    
    //! move to equilibrate the cross links by force equilibrium
    MoveForceEquilibrium move;
    
    //! random number generator 
    RandomNumberGenerators rng;
    
    //! Hold the value of lattice size in X
    uint32_t _boxX;

    //! Hold the value of lattice size in Y
    uint32_t _boxY;

    //! Hold the value of lattice size in Z
    uint32_t _boxZ;
    
    //! Functions for folding absolute coordinates into the lattice in X
    double foldBackX(double value) const;

    //! Functions for folding absolute coordinates into the lattice in Y
    double foldBackY(double value) const;

    //! Functions for folding absolute coordinates into the lattice in Z
    double foldBackZ(double value) const;
};
/**
 * Fold back the absolute coordinate into the relative coordinate in X by modulo operation.
 *
 * @param value absolute coordinate in X
 * @return \a double relative coordinate in X
 */
template <class IngredientsType>
inline double UpdaterForceBalancedPosition<IngredientsType>::foldBackX(double value) const{
    while ( value > (_boxX) ) value-=_boxX;
    while ( value < 0       ) value+=_boxX;
    return value;
}

/**
 * Fold back the absolute coordinate into the relative coordinate in Y by modulo operation.
 *
 * @param value absolute coordinate in Y
 * @return \a double relative coordinate in Y
 */
template <class IngredientsType>
inline double UpdaterForceBalancedPosition<IngredientsType>::foldBackY(double value) const{
    while ( value > (_boxZ) ) value-=_boxY;
    while ( value < 0       ) value+=_boxY;
    return value;
}

/**
 * Fold back the absolute coordinate into the relative coordinate in Z by modulo operation.
 *
 * @param value absolute coordinate in Z
 * @return \a double relative coordinate in Z
 */
template <class IngredientsType>
inline double UpdaterForceBalancedPosition<IngredientsType>::foldBackZ(double value) const{
    while ( value > (_boxZ) ) value-=_boxZ;
    while ( value < 0       ) value+=_boxZ;
    return value;
}

template <class IngredientsType>
void UpdaterForceBalancedPosition<IngredientsType>::initialize(){
    uint32_t k(0);
    for (size_t i = 0 ; i < ing.getMolecules().size();i++){
        if( ing.getMolecules()[i].isReactive()  && ing.getMolecules()[i].getNumMaxLinks() > 2 ){
          CrossLinkIDs.push_back(i);
          if (k < 10 )
            std::cout << "CrosslinksID: " << CrossLinkIDs.back() <<std::endl;
          k++;
        }	 
    }
    NCrossLinks=CrossLinkIDs.size();
    _boxX=ing.getBoxX(); 
    _boxY=ing.getBoxY();
    _boxZ=ing.getBoxZ();
}

template <class IngredientsType>
bool UpdaterForceBalancedPosition<IngredientsType>::execute(){
    std::cout << "UpdaterForceBalancedPosition::execute(): Start equilibration" <<std::endl;
    double avShift(threshold*1.1);
    uint32_t StartMCS(ing.getMolecules().getAge());
    while (avShift > threshold  ){
        double NSuccessfulMoves(0.);
        avShift=0.0;
        for (uint32_t i =0 ; i<NCrossLinks ; i++){
            uint32_t RandomMonomer(CrossLinkIDs[ rng.r250_rand32() % NCrossLinks]);
            move.init(ing, RandomMonomer);
            if(move.check(ing)){
              move.apply(ing);
              avShift+=move.getShiftVector().getLength();
              NSuccessfulMoves++;
            }
        }
        avShift/=(NSuccessfulMoves);
        ing.modifyMolecules().setAge(ing.getMolecules().getAge()+1);
        std::cout << "MCS: " << ing.getMolecules().getAge() << "  and average shift: " << avShift<< std::endl;
    }
    for(size_t i = 0; i < NCrossLinks;i++){
        uint32_t ID(CrossLinkIDs[i]);
        ing.modifyMolecules()[ID].modifyVector3D().setAllCoordinates(
            foldBackX(ing.modifyMolecules()[ID].getX()), 
            foldBackY(ing.modifyMolecules()[ID].getY()), 
            foldBackZ(ing.modifyMolecules()[ID].getZ()) );
    }
    std::cout << "Finish equilibration with average shift per cross link < " << avShift << " after " << ing.getMolecules().getAge()-StartMCS <<std::endl;
    ing.modifyMolecules().setAge(StartMCS);
    return false;
}
#endif /* LEMONADE_PM_UPDATER_UPDATERFORCEBALANCEPOSITION_H*/