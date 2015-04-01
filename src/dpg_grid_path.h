/*
Implementação do plugin Distributed Independent Random Fields SGSIM.

(c) 2014, LPM/UFRGS, Luiz Gustavo Rasera, Péricles Lopes Machado

BASED ON AR2GEMS SGSIM IMPLEMENTATION
*/


/* -----------------------------------------------------------------------------
** Copyright (c) 2012 Advanced Resources and Risk Technology, LLC
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


#ifndef __PLUGINS_LPM_UFRGS_DPG_GRID_ORDERED_PATH_H__
#define __PLUGINS_LPM_UFRGS_DPG_GRID_ORDERED_PATH_H__
 
#include <grid/common.h>
#include <math/gstlvector.h> 
#include <grid/geovalue.h>  
#include <grid/geostat_grid.h>
#include <grid/gval_iterator.h>
#include <grid/grid_property.h>
#include <grid/grid_region.h>

#include <vector>  
#include <qmessagebox.h>
#include <qthread.h>
#include <geostat/common.h>
#include <geostat/geostat_algo.h> 
#include <geostat/utilities.h> 
#include <grid/geostat_grid.h> 
#include <grid/property_copier.h> 

#include <grid/two_point_statistics.h>
#include <geostat/kriging_system.h>
#include <geostat/kriging_constraint.h>
#include <geostat/gaussian_cdf_kestimator.h>

#include <qtimer.h>
#include <QTime>
#include <qtextstream.h>
#include <qmessagebox.h>

#include <geostat/sgsim.h>
#include <geostat/parameters_handler.h>
#include <geostat/utilities.h>
#include <utils/gstl_messages.h>
#include <utils/error_messages_handler.h>
#include <utils/string_manipulation.h>
#include <grid/geostat_grid.h>
#include <grid/combined_neighborhood.h>
#include <grid/gval_iterator.h>
#include <grid/cartesian_grid.h>
#include <grid/point_set.h>
#include <grid/grid_path.h>
#include <utils/manager_repository.h>

#include <qmutex.h>

#include <iterator>
#include <vector>
#include <fstream>
#include <random>
#include <algorithm>

#include <grid/reduced_grid.h>


#include "common.h"
 
typedef std::vector<GsTLInt> DPG_Nodes;

/** 
 * Grid_ordered_path class represents a traversal path on a Geostat_grid 
 * in the order from the beginning of the grid to the end of the grid.
 * The path can be retricted to a region.
 */
class PLUGINS_LPM_UFRGS_DECL DPG_Grid_path
{
public:
  typedef Gval_iterator<TabularMapIndex> iterator;
  typedef Gval_const_iterator<TabularMapIndex> const_iterator;


public:
  // get the reference to the begin of the path
  iterator begin();
  const_iterator begin() const;

  // get the end of the path
  iterator end();
  const_iterator end() const;

  // get the size of the path
  GsTLInt size() const;

  // get the node_id (used in the grid) of a node on the path. The argument _path_index is the index on the path
  GsTLInt node_id( GsTLInt _path_index ) const;

  iterator at(GsTLInt _path_index);
  const_iterator at(GsTLInt _path_index) const;

  // get the Geovalue associated with a node on the path. The argument _path_index is the index on the path
  Geovalue geovalue( GsTLInt _path_index ) const;

  // randomize the path.  A new ordering is created each time this function is called
  // The constructor create a ordered path
  virtual void randomize();

  void set_multi_grid_level(int level, double dx, double dy, double dz);
  void set_multi_grid_level(const DPG_Nodes& _nodes, int level, double dx, double dy, double dz);

  virtual bool set_property(std::string prop_name);

  void setSeed(int seed_) {
	  seed = seed_;
	  //tgen.seed(seed);
	 tgen.seed(seed);
  }

protected:
  std::vector<GsTLInt> DPG_Grid_path_;
  Geostat_grid * grid_;
  Grid_continuous_property * prop_;
  Grid_region * region_;

  class MyGen {
  public:
	  void seed(int seed_) {
		  r.seed(seed_);
	  }
	  size_t operator()(size_t N) {
		  return r() % N;
	  }
  private:
	  std::default_random_engine r;
  };

  // engine;
  MyGen tgen;
  //std::mt19937 tgen;
  //Rand48_generator tgen;
  //STL_generator tgen;
  int seed;

public:
  DPG_Grid_path();
//  DIRF_Grid_path(Geostat_grid * _grid, Grid_region * _region = 0);
  DPG_Grid_path(Geostat_grid * _grid, Grid_continuous_property * _prop, Grid_region * _region =0, bool ini = true);
  DPG_Grid_path(const DPG_Nodes& _nodes, Geostat_grid * _grid, Grid_continuous_property * _prop, Grid_region * _region = 0, bool ini = true);
  virtual ~DPG_Grid_path(void);
};


template<class TargetCdfType, class DataCdfType>
void cdf_transform( Grid_continuous_property* prop, size_t io, size_t inc,
	      const DataCdfType& from, const TargetCdfType& to) {
  
  for(size_t i = io; i < prop->size(); i += inc) {
	  double P = from.prob(prop->get_value(i));
      prop->set_value(to.inverse(P), i);
    }
}





#endif
