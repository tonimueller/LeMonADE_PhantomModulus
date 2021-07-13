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
#ifndef LEMONADE_PM_UPDATER_UPDATERAFFINEDEFORMATION_H
#define LEMONADE_PM_UPDATER_UPDATERAFFINEDEFORMATION_H


#include <LeMonADE/updater/AbstractUpdater.h>
#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE_PM/utility/neighborX.h>
#include <vector>
#include <cmath>
 /**
 * @class UpdaterAffineDeformation
 * @brief This class performs a affine deformation accordign to the stretching factor (lambda)
 * for the positions and the jump vector. 
 * @tparam IngredientsType
 */

 template <class IngredientsType >
class UpdaterAffineDeformation:public AbstractUpdater
{
public:
    //! constructor for UpdaterAffineDeformation
    UpdaterAffineDeformation(IngredientsType& ing_, double stretching_factor_ ):
    ing(ing_),stretching_factor(stretching_factor_),stretching_factor_XY(std::sqrt(stretching_factor)){};
    
    virtual void initialize();
    virtual bool execute(){return false;};
    virtual void cleanup(){};  

private:
    //!copy of the main container for the system informations 
    IngredientsType& ing;
    
    //! threshold for the certainty 
    double stretching_factor;
    double stretching_factor_XY;
    VectorDouble3 deform( VectorDouble3 vec ){
        return VectorDouble3(vec.getX()*stretching_factor,vec.getY()*stretching_factor_XY,vec.getZ()*stretching_factor_XY);
    }
};
template <class IngredientsType>
void UpdaterAffineDeformation<IngredientsType>::initialize(){
    std::cout << "UpdaterAffineDeformation<IngredientsType>::initialize():"<< std::endl;
    //adjusting the box size is not neccessary, because it is used only once in the FeatureCrosslinkConnections*
    //there the jump vectors are calculated
    
    //The jump vector needs to be adjusted according the deformation labmda! 
    auto CrossLinkIDs(ing.getCrosslinkIDs());
    for (uint32_t i =0 ; i<CrossLinkIDs.size() ; i++){
        uint32_t ID(CrossLinkIDs[i]);
        if (i < 20 ) 
            std::cout << "ID="<< i << "initial jump " << ing.getMolecules()[ID].getVector3D() << " ";  
        std::vector<neighborX> Neighbors(ing.getCrossLinkNeighborIDs(ID) );
        int32_t number_of_neighbors(Neighbors.size());
        if (number_of_neighbors > 0) {
            Neighbors[i].jump= deform(Neighbors[i].jump);
        }
        if (i < 20 ) 
            std::cout << "deformed jump " << ing.getMolecules()[ID].getVector3D() << "\n";  
    }

    //The initial posistions of all monomers needs to be adjusted according the deformation labmda! 
    for(size_t i = 0; i< ing.getMolecules().size();i++){
        if (i < 20 ) 
            std::cout << "ID="<< i << "initial position " << ing.getMolecules()[i].getVector3D() << " ";  
        ing.modifyMolecules()[i].modifyVector3D()=deform(ing.modifyMolecules()[i].getVector3D());
        if (i < 20 ) 
            std::cout << "deformed position " << ing.getMolecules()[i].getVector3D() << "\n";  
    }
    std::cout << "UpdaterAffineDeformation<IngredientsType>::initialize():done.\n";
}
#endif /*LEMONADE_PM_UPDATER_UPDATERAFFINEDEFORMATION_H*/