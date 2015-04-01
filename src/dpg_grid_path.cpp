/*
Implementação do plugin Distributed Independent Random Fields SGSIM.

(c) 2014, LPM/UFRGS, Luiz Gustavo Rasera, Péricles Lopes Machado

BASED ON AR2GEMS SGSIM IMPLEMENTATION
*/



/* -----------------------------------------------------------------------------
** Copyright (c) 2013 Advanced Resources and Risk Technology, LLC
** All rights reserved.
**
** This file is part of Advanced Resources and Risk Technology, LLC (AR2TECH) 
** version of the open source software sgems.  It is a derivative work by 
** AR2TECH (THE LICENSOR) based on the x-free license granted in the original 
** version of the software (see notice below) and now sublicensed such that it 
** cannot be distributed or modified without the explicit and written permission 
** of AR2TECH.
**
** Only AR2TECH can modify, alter or revoke the licensing terms for this 
** file/software.
**
** This file cannot be modified or distributed without the explicit and written 
** consent of AR2TECH.
**
** Contact Dr. Alex Boucher (aboucher@ar2tech.com) for any questions regarding
** the licensing of this file/software
**
** The open-source version of sgems can be downloaded at 
** sourceforge.net/projects/sgems.
** ----------------------------------------------------------------------------*/

//#include <math/random_numbers.h>

#include "dpg_grid_path.h"

#include <QMutex>

QMutex randomize_lock;

DPG_Grid_path::iterator DPG_Grid_path::at(GsTLInt _path_index)
{
    return DPG_Grid_path::iterator(grid_, prop_,
        _path_index, DPG_Grid_path_.size(),
        TabularMapIndex(&DPG_Grid_path_));
}

DPG_Grid_path::const_iterator DPG_Grid_path::at(GsTLInt _path_index) const
{
    return DPG_Grid_path::const_iterator(grid_, prop_,
        _path_index, DPG_Grid_path_.size(),
        TabularMapIndex(&DPG_Grid_path_));

}

DPG_Grid_path::iterator
DPG_Grid_path::begin()
{
  return DPG_Grid_path::iterator( grid_, prop_,
                                   0, DPG_Grid_path_.size(),
                                   TabularMapIndex(&DPG_Grid_path_) );
}

DPG_Grid_path::const_iterator
DPG_Grid_path::begin() const
{
  return DPG_Grid_path::const_iterator( grid_, prop_,
                                   0, DPG_Grid_path_.size(),
                                   TabularMapIndex(&DPG_Grid_path_) );
}

DPG_Grid_path::iterator
DPG_Grid_path::end()
{
  return DPG_Grid_path::iterator( grid_, prop_,
                                  DPG_Grid_path_.size(), DPG_Grid_path_.size(),
                                  TabularMapIndex(&DPG_Grid_path_) );
}

DPG_Grid_path::const_iterator
DPG_Grid_path::end() const
{
  return DPG_Grid_path::const_iterator( grid_, prop_,
                                  DPG_Grid_path_.size(), DPG_Grid_path_.size(),
                                  TabularMapIndex(&DPG_Grid_path_) );
}

GsTLInt 
DPG_Grid_path::size() const
{
  return DPG_Grid_path_.size();
}

GsTLInt 
DPG_Grid_path::node_id( GsTLInt _path_index ) const
{
  return grid_->node_id( DPG_Grid_path_[_path_index] );
}

Geovalue 
DPG_Grid_path::geovalue( GsTLInt _path_index ) const
{
  GsTLInt node_id = grid_->node_id( DPG_Grid_path_[_path_index] );
  return Geovalue(grid_, prop_, node_id);
}

DPG_Grid_path::DPG_Grid_path() : grid_(0),prop_(0)  {}

DPG_Grid_path::DPG_Grid_path(Geostat_grid * _grid, Grid_continuous_property * _prop, Grid_region * _region, bool ini)
{
  //if ( !_grid || !_prop ) return;
  this->grid_ = _grid;
  this->prop_ = _prop;
  this->region_ = _region;
  
  if (!ini) return;

  // init random path
  if(_region == 0) {
    DPG_Grid_path_.resize( grid_->size() );
    for( int i = 0; i < grid_->size(); i++ ) {
      DPG_Grid_path_[i] = i;
    }
  }
  else {
    DPG_Grid_path_.clear();
    for( int i = 0; i < grid_->size(); i++ ) {
      if(_region->is_inside_region(i)) DPG_Grid_path_.push_back( i );
    }
  }
}

DPG_Grid_path::DPG_Grid_path(const DPG_Nodes& _nodes, Geostat_grid * _grid, Grid_continuous_property * _prop, Grid_region * _region, bool ini)
{
  //if ( !_grid || !_prop ) return;
  this->grid_ = _grid;
  this->prop_ = _prop;
  this->region_ = _region;

  if (!ini) return;

  // init random path
  if(_region == 0) {
    DPG_Grid_path_.resize( _nodes.size() );
    for( int i = 0; i < _nodes.size(); i++ ) {
      DPG_Grid_path_[i] = _nodes[i];
    }
  }
  else {
    DPG_Grid_path_.clear();
    for( int i = 0; i < _nodes.size(); i++ ) {
      if(_region->is_inside_region(_nodes[i])) DPG_Grid_path_.push_back( _nodes[i] );
    }
  }
}

/*
DIRF_Grid_path::DIRF_Grid_path(Geostat_grid * _grid,  Grid_region * _region)
{
  if ( !_grid  ) return;
  this->grid_ = _grid;
  this->prop_ = 0;
  
  // init random path
  if(_region == 0) {
    DIRF_Grid_path_.resize( grid_->size() );
    for( int i = 0; i < grid_->size() ; i++ ) {
      DIRF_Grid_path_[i] = i;
    }
  }
  else {
    DIRF_Grid_path_.reserve( grid_->size() );
    for( int i = 0; i <  grid_->size() ; i++ ) {
      if(_region->is_inside_region(i)) DIRF_Grid_path_.push_back( i );
    }
  }
}
*/

void DPG_Grid_path::set_multi_grid_level(const DPG_Nodes& _nodes, int level, double dx, double dy, double dz)
{
    DPG_Grid_path_.clear();

	if (region_ == 0) {
		for (int i = 0; i < _nodes.size(); ++i) {
			int ix = static_cast<int>(grid_->location(_nodes[i]).x() / dx);
			int iy = static_cast<int>(grid_->location(_nodes[i]).y() / dy);
			int iz = static_cast<int>(grid_->location(_nodes[i]).z() / dz);
			if (ix % level == 0 &&
				iy % level == 0 &&
				iz % level == 0) {
                DPG_Grid_path_.push_back(_nodes[i]);
			}
		}
	}
	else {
		for (int i = 0; i < _nodes.size(); ++i) {
			int ix = static_cast<int>(grid_->location(_nodes[i]).x() / dx);
			int iy = static_cast<int>(grid_->location(_nodes[i]).y() / dy);
			int iz = static_cast<int>(grid_->location(_nodes[i]).z() / dz);
			if (ix % level == 0 &&
				iy % level == 0 &&
				iz % level == 0) {
                if (region_->is_inside_region(_nodes[i])) DPG_Grid_path_.push_back(_nodes[i]);
			}
		}

	}
}

void DPG_Grid_path::set_multi_grid_level(int level, double dx, double dy, double dz)
{
	if (region_ == 0) {
        DPG_Grid_path_.clear();
		for (int i = 0; i < grid_->size(); ++i) {

			int ix = static_cast<int>(grid_->location(i).x() / dx);
			int iy = static_cast<int>(grid_->location(i).y() / dy);
			int iz = static_cast<int>(grid_->location(i).z() / dz);

			if (ix % level == 0 &&
				iy % level == 0 &&
				iz % level == 0) {
                DPG_Grid_path_.push_back(i);
			}
		}
	}
	else {
        DPG_Grid_path_.clear();
		for (int i = 0; i < grid_->size(); ++i) {
			int ix = static_cast<int>(grid_->location(i).x() / dx);
			int iy = static_cast<int>(grid_->location(i).y() / dy);
			int iz = static_cast<int>(grid_->location(i).z() / dz);

			if (ix % level == 0 &&
				iy % level == 0 &&
				iz % level == 0) {
                if (region_->is_inside_region(i)) DPG_Grid_path_.push_back(i);
			}
		}
	}
}

DPG_Grid_path::~DPG_Grid_path(void)
{
  this->grid_ = NULL;
  this->prop_ = NULL;
}

bool DPG_Grid_path::set_property(std::string prop_name){
  this->prop_ = grid_->property(prop_name);
  return prop_ != 0;
}

void DPG_Grid_path::randomize(){
	//STL_generator gen;
	//STL_generator gen;
  std::random_shuffle( DPG_Grid_path_.begin(), DPG_Grid_path_.end(), tgen );
}


