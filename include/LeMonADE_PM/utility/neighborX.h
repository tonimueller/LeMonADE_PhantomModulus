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

#ifndef LENONADE_PM_UTILITY_NEIGHBORX_H
#define LENONADE_PM_UTILITY_NEIGHBORX_H


//! crosslink neighbor with information about the ID the distance in segments 
// and the jumps in multiples of the box sizes for connections across 
//periodic boundaries 
struct neighborX{
    // neighborX(){};
    // neighborX(uint32_t ID_, uint32_t segDistance_, VectorInt3 jump_):ID(ID_), segDistance(segDistance_), jump(jump_){};
    //id of the neighbor 
    uint32_t ID;
    //segmental distance between the neighbors 
    uint32_t segDistance;
    //jump vector across periodic boundaries  
    VectorInt3 jump;
};

#endif /*LENONADE_PM_UTILITY_NEIGHBORX_H*/