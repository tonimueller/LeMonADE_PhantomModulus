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
#ifndef LEMONADE_PM_UPDATER_UPDATERFORCEBALANCEPOSITIONTENDOMER_H
#define LEMONADE_PM_UPDATER_UPDATERFORCEBALANCEPOSITIONTENDOMER_H


#include <LeMonADE/updater/AbstractUpdater.h>
#include <vector>
 /**
 * @class UpdaterForceBalancedPositionTendomer
 * @tparam IngredientsType
 */

 template <class IngredientsType, class moveType >
class UpdaterForceBalancedPositionTendomer:public AbstractUpdater
{
public:
    //! constructor for UpdaterForceBalancedPositionTendomer
    UpdaterForceBalancedPositionTendomer(IngredientsType& ing_, double threshold_, double factor_=0.995 ):
    ing(ing_),threshold(threshold_),factor(factor_){};
    
    virtual void initialize(){relaxationParameter=move.getRelaxationParameter();}// init(move);};
    bool execute();
    virtual void cleanup(){};  

    void setFilename(const std::string filename) {move.setFilename(filename); }
    void setRelaxationParameter( const double relax ) {move.setRelaxationParameter(relax);}
private:
    //!copy of the main container for the system informations 
    IngredientsType& ing;
    
    //! threshold for the certainty 
    double threshold;

    double relaxationParameter;

    double factor;
    
    //! move to equilibrate the cross links by force equilibrium
    moveType move;
    
    //! random number generator 
    RandomNumberGenerators rng;
    
};
template <class IngredientsType, class moveType>
bool UpdaterForceBalancedPositionTendomer<IngredientsType,moveType>::execute(){
    std::cout << "UpdaterForceBalancedPositionTendomer::execute(): Start equilibration" <<std::endl;
    auto counter=0;
    double avShift(threshold*1.1);
    uint32_t StartMCS(ing.getMolecules().getAge());
    setRelaxationParameter(relaxationParameter);
    //! get look up table for the cross link ids to monomer ids
    auto CrossLinkIDs(ing.getCrosslinkIDs());
    //! number of cross links 
    auto NCrossLinks(CrossLinkIDs.size());
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
        if( NSuccessfulMoves>0 ){
            avShift/=(NSuccessfulMoves);
            counter=0;
        }else {
            avShift=threshold*1.1;
            counter++;
        }
        if (counter > 10 ){break;}
        ing.modifyMolecules().setAge(ing.getMolecules().getAge()+1);
        if (ing.getMolecules().getAge() %1000 == 0 ){
            std::cout << "MCS: " << ing.getMolecules().getAge() << "  and average shift: " << avShift << std::endl;
            setRelaxationParameter(move.getRelaxationParameter()*factor );
        }
    }
    std::cout << "Finish equilibration with average shift per cross link < " << avShift << " after " << ing.getMolecules().getAge()-StartMCS <<std::endl;
    ing.modifyMolecules().setAge(StartMCS);
    return false;
}
#endif /* LEMONADE_PM_UPDATER_UPDATERFORCEBALANCEPOSITION_H*/