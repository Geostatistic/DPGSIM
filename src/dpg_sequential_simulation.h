/*
Implementação do plugin Distributed Independent Random Fields SGSIM.

(c) 2014, LPM/UFRGS, Luiz Gustavo Rasera, Péricles Lopes Machado

BASED ON AR2GEMS SGSIM IMPLEMENTATION
*/

#ifndef __DPG_SEQUENTIAL_SIMULATION_H__
#define __DPG_SEQUENTIAL_SIMULATION_H__


#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <Eigen/Core>

#include <utils/gstl_messages.h>
#include <math/continuous_distribution.h>
#include <grid/geostat_grid.h>
#include <geostat/cokriging_system.h>
#include <geostat/kriging_system.h>
#include <math/continuous_distribution.h>

#include <math/continuous_distribution.h>

#include "dpg_grid_path.h"

namespace LPM_UFRGS {




/**
 * Functions in this file behave exactly functions in the same file in GsTL library.
 * The difference is that the functions here use thread-safe method find_neighbors(center, neighbors)
 * of Neighborhood classes in order to support parallelism.
 */






/** This function implements the sequential simulation algorithm
 * @return 0 if no problem were encountered during the simulation,
 * -1 if the execution was aborted, n (n>0) if the ccdf could not
 * be estimated n times.
 *
 * @param [begin,end) is the range of geovalues to be simulated
 * @param neighbors is the neighborhood that will be used to retrieve the 
 * conditioning data at each node to be simulated
 * @param ccdf is the conditional cdf that will be estimated at each 
 * node. <pre>ccdf</pre> is changed at each node.
 * @param estim is a functor the estimates the ccdf.
 * @param marginal is the marginal cdf.
 * @param samp is the functor that draws a realization from the ccdf.
 */
template
<
  class Cdf,
  class CdfEstimator,
  class MarginalCdf,
  class Sampler
>
inline int 
sequential_simulation(
              DPG_Grid_path* path,
			  GsTLInt ini, GsTLInt fin,
		      Neighborhood& neighborhood,
		      Cdf& ccdf,
		      CdfEstimator& estim,
		      const MarginalCdf& marginal,
		      Sampler& samp,
			  bool is_cross_validation = false
		      ) {
  int ok = 0;

  for (GsTLInt p = ini; p < fin; ++p) {
	Geovalue node = path->geovalue(p);

	if (node.is_informed()  && !is_cross_validation) continue;
	if (!node.is_informed() && is_cross_validation) continue;

    Neighbors neighbors;
    neighborhood.find_neighbors( node , neighbors );
    
    DEBUG_PRINT_LOCATION( "center", node.location() );

    if( neighbors.is_empty() ){
      //if we don't have any conditioning data, we simply draw from the
      // marginal
      WRITE_TO_DEBUG_STREAM( "drawing from marginal" << std::endl );
      samp(node, marginal);
    }
    else {
      DEBUG_PRINT_NEIGHBORHOOD( "neighbors", &neighbors );

      int status = estim( node, neighbors, ccdf);

      if(status == 0) {
      	samp(node, ccdf);
      }
      else {
      	// the ccdf could not be estimated. Draw from the marginal
	      samp(node, marginal);
	      ok++;
      }
    }
  }

  return ok;
}







/** This function implements the sequential simulation algorithm.
 * @return 0 if no problem were encountered during the simulation,
 * -1 if the execution was aborted, n (n>0) if the ccdf could not
 * be estimated n times.
 *
 * @param [begin,end) is the range of geovalues to be simulated
 * @param neighbors is the neighborhood that will be used to retrieve the 
 * conditioning data at each node to be simulated
 * @param ccdf is the conditional cdf that will be estimated at each 
 * node. <pre>ccdf</pre> is changed at each node.
 * @param estim is a functor the estimates the ccdf.
 * @param marginal is the marginal cdf.
 * @param samp is the functor that draws a realization from the ccdf.
 * @param notifier is an object that has a <pre>notify()</pre> function.
 * That function is executed after a node has been simulated. It can be 
 * used to give information on the progress of the simulation.
 */
template
<
  class Cdf,
  class CdfEstimator,
  class MarginalCdf,
  class Sampler,
  class Notifier
>
inline int 
sequential_simulation(
            DPG_Grid_path* path,
			GsTLInt ini, GsTLInt fin,
		    Neighborhood& neighborhood,
		    Cdf& ccdf,
		    CdfEstimator& estim,
		    const MarginalCdf& marginal,
		    Sampler& samp,
		    Notifier* notifier,
			bool is_cross_validation = false
		      ) {

  int bad = 0;

  for (GsTLInt p = ini; p < fin; ++p) {
	 Geovalue node = path->geovalue(p);
	 if (node.is_informed() && !is_cross_validation) continue;
	 if (!node.is_informed() && is_cross_validation) continue;


    Neighbors neighbors;
    neighborhood.find_neighbors(node, neighbors);

    DEBUG_PRINT_LOCATION( "center", begin->location() );

    if( neighbors.is_empty() ){
      //if we don't have any conditioning data, we simply draw from the
      // marginal
      WRITE_TO_DEBUG_STREAM( "drawing from marginal" << std::endl );
      samp(node, marginal);
    }
    else {
      DEBUG_PRINT_NEIGHBORHOOD( "neighbors", &neighbors );

      int status = estim( node, neighbors, ccdf);

      if(status == 0) {
          samp(node, ccdf);
      }
      else {
      	// the ccdf could not be estimated. Draw from the marginal
        WRITE_TO_DEBUG_STREAM( "Can't estimate ccdf. drawing from marginal\n" );
        samp(node, marginal);
      	bad++;
      }
    }

    if(notifier) {
      if( !notifier->notify() ) return -1;
    }
  }

  return bad;
}





template
<
  class Random_number_generator_t
>
inline int
sequential_cokriging_simulation(
                    DPG_Grid_path* path,
					GsTLInt ini, GsTLInt fin,
                    const std::vector<Neighborhood *> & neighborhoods, 
                    const std::vector<std::pair<int, Grid_continuous_property*> > & simul_props,
                    Cokriging_system & kriging_system,
                    Random_number_generator_t& gen,
					bool is_cross_validation = false
                        ) {
  int bad = 0;
  const int number_variables = neighborhoods.size();
  const int number_target_variables = simul_props.size();

  Gaussian_distribution gauss(0, 1);

  // Build unknown-to-unknown matrix (same matrix for all nodes)
//  Cokriging_system::Kriging_matrix unknown2unknown_matrix;
//  kriging_system.build_colocated_unknown2unknown_matrix(number_target_variables, unknown2unknown_matrix);

  for (GsTLInt p = ini; p < fin; ++p) {
	  Geovalue node = path->geovalue(p);


	  
    // 1) Find neighbors
    kriging_system.localize(node);
    std::vector<Neighbors> neighbors_set(number_variables);
    int nb_conditioning_data = 0;
    for(int i = 0; i < neighborhoods.size(); i++) {
      neighborhoods[i]->find_neighbors( node , neighbors_set[i] );
      nb_conditioning_data += neighbors_set[i].size();
    }

    //Chekc if all targetvar are informed
    bool all_informed = true;
    for(int i=0; i<number_target_variables && all_informed; ++i ) {
      if( !simul_props[i].second->is_informed( node.node_id() ) ) {
        all_informed  =false;
      }
    }

    if (all_informed && !is_cross_validation) continue;
	if (!all_informed && is_cross_validation) continue;

 
    // Build unknown-to-unknown matrix (same matrix for all nodes)
    Cokriging_system::Kriging_matrix unknown2unknown_matrix;
    kriging_system.build_colocated_unknown2unknown_matrix(number_target_variables, unknown2unknown_matrix);


    // 2) Build the LHS matrix
    Cokriging_system::Kriging_matrix LHS;
    kriging_system.build_data2data_matrix(neighbors_set, LHS); 
    
    // 3) Build RHS matrix
    Cokriging_system::Kriging_matrix RHS = Cokriging_system::Kriging_matrix::Zero(LHS.rows(), number_target_variables);
    for (int var = 0; var < number_target_variables; var++) {
      Cokriging_system::Kriging_vector b;
      kriging_system.build_data2unknown_vector(node, neighbors_set, simul_props[var].first, b);
      RHS.col(var) = b;
    }

    // 4) Build combined matrix (LU matrix)
    const int combine_matrix_dimension = LHS.rows() + number_target_variables;
    Cokriging_system::Kriging_matrix combined_matrix = Cokriging_system::Kriging_matrix::Zero(combine_matrix_dimension, combine_matrix_dimension);
    combined_matrix.block( 0, 0, LHS.rows(), LHS.cols() ) = LHS;
    combined_matrix.block( 0, LHS.cols(), LHS.rows(), number_target_variables ) = RHS;
    combined_matrix.block( LHS.rows(), 0, number_target_variables, LHS.rows() ) = RHS.transpose();
    combined_matrix.block( LHS.rows(), LHS.cols(), number_target_variables, number_target_variables ) = unknown2unknown_matrix;

    Eigen::LLT<Cokriging_system::Kriging_matrix> llt(combined_matrix); // compute the Cholesky decomposition of combined_matrix
    Cokriging_system::Kriging_matrix L = llt.matrixL(); // retrieve factor L  in the decomposition
    
    // 5) Solve the linear equation system for L without unknown variables
    Cokriging_system::Kriging_matrix L_data = L.topLeftCorner(L.rows() - number_target_variables, L.cols() - number_target_variables);
    Cokriging_system::Kriging_vector Z_data(L_data.rows());
    Cokriging_system::Kriging_vector W(L.rows());
    // build Z-data
    int current_row = 0;
    for (int var = 0; var < number_variables; var++) {
      // will put a function here for parallel
      for (int n = 0; n < neighbors_set[var].size(); n++, current_row++) {
        Z_data(current_row) = neighbors_set[var][n].property_value();
      }
    }
    // solve
    Cokriging_system::Kriging_vector W_data;
    W.head(Z_data.rows()) = L_data.triangularView<Eigen::Lower>().solve(Z_data);
    
    // 6) Solve for unknown variables to get simulated values
    // generate random values
    for (int var = 0; var < number_target_variables; var++) {
      double rand = gauss.quantile(gen.operator()());
      W(nb_conditioning_data + var) = rand;
    }
    Cokriging_system::Kriging_vector Z_simulated = L.bottomRows(number_target_variables) * W;
    
    // 7) Finally, set to each simulation property
    int center_node_id = node.node_id();
    for (int var = 0; var < number_target_variables; var++) {
      simul_props[var].second->set_value(Z_simulated(var), center_node_id);
    }
  }
  
  return bad;
}




template
<
  class Sampler,
  class Notifier
>
inline int
sequential_gaussian_cosimulation(
            DPG_Grid_path* path,
			GsTLInt ini, GsTLInt fin,
			const Neighborhood* primary_neighborhood, 
			const Neighborhood* secondary_neighborhood, 
			Grid_continuous_property* simul_prop,
			Cokriging_system & kriging_system,
			Sampler& samp,
			Notifier* notifier,
			bool is_cross_validation = false
                        ) {
  int bad = 0;

  Gaussian_distribution marginal(0, 1);

 
  for (GsTLInt p = ini; p < fin; ++p) {
      DPG_Grid_path::iterator begin = path->at(p);
    if(notifier) notifier->notify();
    kriging_system.localize(*begin);

    // 1) Find neighbors
    std::vector<Neighbors> neighbors_set;
    neighbors_set.push_back(Neighbors());
    neighbors_set.push_back(Neighbors());
    primary_neighborhood->find_neighbors( *begin , neighbors_set[0] );
    secondary_neighborhood->find_neighbors( *begin , neighbors_set[1] );

    int nb_conditioning_data = neighbors_set[0].size() + neighbors_set[1].size();

    if(nb_conditioning_data == 0) {
      samp(*begin, marginal);
      continue;
    }

    // check if all target center are informed



    Cokriging_system::Kriging_matrix LHS;
    kriging_system.build_data2data_matrix(neighbors_set, LHS); 

    Cokriging_system::Kriging_solve_result result;
    kriging_system.solve(*begin,neighbors_set, 0,LHS,result);

    //Get the data point
    Cokriging_system::Kriging_vector z_data(LHS.rows());
    int n_primary_data = neighbors_set[0].size();
    for (int n = 0; n < n_primary_data; n++) {
      z_data(n) = neighbors_set[0][n].property_value();
    }
    for (int n = 0; n < neighbors_set[1].size(); n++) {
      z_data(n + n_primary_data) = neighbors_set[1][n].property_value();
    }

    double mean = z_data.dot(result.weights);
    double var = kriging_system.compute_kriging_variance(0,result);

    if(var <= 0) {
      begin->set_property_value(mean);
    }
    else {
      Gaussian_distribution cond_distribution(mean, std::sqrt(var));
      samp(*begin, cond_distribution);
    }
  }

  return 0;

}



} // namespace ar2gems_geostat

#endif
