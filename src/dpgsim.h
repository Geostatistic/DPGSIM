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



#ifndef __PLUGINS_LPM_UFRGS_DPGSGIM_H___
#define __PLUGINS_LPM_UFRGS_DPGSGIM_H___

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

//#include <GsTL/geometry/covariance.h> 
#include <GsTL/cdf/non_param_cdf.h> 
//#include <GsTL/kriging/kriging_constraints.h> 
//#include <GsTL/kriging/kriging_combiner.h> 
#include <GsTL/utils/smartptr.h> 
#include <string> 
#include <geostat/parameters_handler_impl.h>
#include <QDomDocument>

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
//#include <math/random_numbers.h>
#include <GsTL/math/random_number_generators.h>
#include <appli/utilities.h>

#include <geostat/sequential_simulation.h>
#include <geostat/kriging_constraint.h>
#include <math/continuous_distribution.h>
#include <math/non_parametric_distribution.h>


#include <GsTL/cdf/gaussian_cdf.h>
#include <GsTL/sampler/monte_carlo_sampler.h>
#include <GsTL/univariate_stats/cdf_transform.h>
#include <GsTL/univariate_stats/build_cdf.h>

#include <grid/point_set_neighborhood.h>

#include <qmutex.h>

#include <iterator>
#include <vector>
#include <algorithm>
#include <fstream>

#include <grid/reduced_grid.h>

#include <omp.h>


#include "common.h"
#include "dpg_grid_path.h"
#include "dpg_sequential_simulation.h"

#include "report.h"

class Neighborhood; 
class Parameters_handler; 
class Error_messages_handler; 
class DPGSIMProcess;

template<class TargetCdfType, class DataCdfType>
class  DPGSIMCDFTransform : public QThread {
public slots:
    DPGSIMCDFTransform(size_t io, size_t inc,
							const DataCdfType* from, const TargetCdfType* to, Grid_continuous_property* prop = 0)
							: io(io), inc(inc), from(from), to(to), prop(prop) {}

	void setProp(Grid_continuous_property* prop)
	{
		this->prop = prop;
	}

	void run()
	{
		cdf_transform(prop, io, inc, *from, *to);
	}

private:
	size_t io; 
	size_t inc;
	const DataCdfType* from;
	const TargetCdfType* to;
	Grid_continuous_property* prop;
};

class PLUGINS_LPM_UFRGS_DECL DPGSgsim : public Geostat_algo {
public:
    DPGSgsim();
    ~DPGSgsim();

    virtual bool initialize( const Parameters_handler* parameters,
                             Error_messages_handler* errors,
                             Progress_notifier* notifier = 0);
    virtual int execute( GsTL_project* proj=0, Progress_notifier* notifier = 0 );

	int execute_sgsim(int id_thread = 0, int n_threads = 1, Progress_notifier* notifier = 0);

    virtual std::string name() const { return "lpm_ufrgs_dpgsim"; }

	void set_report_is_open(bool open) { report_is_open_ = open; }
	bool report_is_open() { return report_is_open_; }

public:
    static Named_interface* create_new_interface( std::string& );
	
protected:
	 int total_steps;
	 int frequency;
	 int nIRFX, nIRFY, nIRFZ;
	 int n_Revisitations;
	 double dx, dy, dz;

   typedef Geostat_grid::location_type Location;


    std::vector<Geostat_grid*> simul_grid_;
	std::vector<Multi_property_generator*> multireal_property_;

	std::vector<Multi_property_generator*> multireal_stage_property_;

    std::vector<Geostat_grid*> harddata_grid_;
    std::vector<Grid_continuous_property*> harddata_property_;

    std::string harddata_property_name_;
	std::string property_name_;
	std::string simul_grid_name_;
    //  Grid_initializer* initializer_;

    std::vector<SmartPtr<Property_copier> > property_copier_;

    std::vector<SmartPtr<Neighborhood> > neighborhood_;
	//vector<SmartPtr<Neighborhood> > my_neigh;

    long int seed_;
    int nb_of_realizations_;
	int nb_of_parallel_realizations_;
	bool execute_cross_validation_;

    std::vector<Two_point_nested_structure> covar_;
    std::vector<Kriging_system> kriging_system_;
    std::vector<Kriging_estimator*> kriging_estimator_;


    bool use_target_hist_;
    //geostat_utils::NonParametricCdfType target_cdf_;
    std::vector<SmartPtr<Continuous_distribution> > target_cdf_;
    bool clear_temp_properties_;

    // Set up the regions
    std::vector<Grid_region*> target_grid_region_;
    std::vector<Grid_region*> hd_grid_region_;

	int n_threads_;


protected:

    void clean(int id_thread, Grid_continuous_property* prop = 0 );
    void clean( std::vector<DPG_Grid_path* >& path );
    void clean( std::vector<DPGSIMProcess*>& threads);
    void clean( std::vector<DPGSIMCDFTransform<Continuous_distribution, Gaussian_cdf>*>* threads);

    DPGSGSIM_report* report_;
	QMutex mutex_;
	bool report_is_open_;
};

class PLUGINS_LPM_UFRGS_DECL DPGExecuteSGSIMProcess : public QThread {
public:
    DPGExecuteSGSIMProcess(DPGSgsim* sgsim = 0, int id_thread = 0, int n_threads = 1, Progress_notifier* notifier = 0)
	: sgsim_(sgsim), id_thread_(id_thread), n_threads_(n_threads), notifier_(notifier)  {

	}

	void run()
	{
		sgsim_->execute_sgsim(id_thread_, n_threads_, notifier_);
	}

private:
    DPGSgsim* sgsim_;
	int id_thread_, n_threads_;
	Progress_notifier* notifier_;
};

class PLUGINS_LPM_UFRGS_DECL DPGSIMProcess : public QThread {
public slots:
    DPGSIMProcess(
					 int size_rev = 0,
					 int n_rev = 0,
					 size_t n_threads = 0, 
				     size_t id_thread = 0, 
                     std::vector<DPG_Grid_path* >* fields_path = 0,
					 std::vector<int>* fields = 0, 
					 Grid_continuous_property* prop = 0, 
					 SmartPtr<Neighborhood>* neighborhood_ = 0,
					 Gaussian_cdf* ccdf = 0,
					 Gaussian_cdf* marginal = 0,
                     Monte_carlo_sampler_t< Mersenne_generator >* sampler = 0,
					 Gaussian_cdf_kestimator* cdf_estimator = 0,
					 Progress_notifier* notifier = 0,
					 bool execute_cross_validation = false) 
	: size_rev(size_rev), n_rev(n_rev), n_threads(n_threads), id_thread(id_thread), fields_path(fields_path), 
	  fields(fields), prop(prop), neighborhood_(neighborhood_),
	  ccdf(ccdf), marginal(marginal), 
	  sampler(sampler), 
	  cdf_estimator(cdf_estimator), notifier(notifier),
	  execute_cross_validation(execute_cross_validation) {
	}

    ~DPGSIMProcess() {
		if (sampler) delete sampler;
		if (ccdf) delete ccdf;
		if (marginal) delete marginal;
		if (cdf_estimator) delete cdf_estimator;
	}

	void setValues( int size_rev = 0,
					int n_rev = 0,
		            size_t n_threads = 0, 
				    size_t id_thread = 0, 
                    std::vector<DPG_Grid_path* >* fields_path = 0,
					std::vector<int>* fields = 0, 
					Grid_continuous_property* prop = 0, 
					SmartPtr<Neighborhood>* neighborhood_ = 0,
					Gaussian_cdf* ccdf = 0,
					Gaussian_cdf* marginal = 0,
                    Monte_carlo_sampler_t< Mersenne_generator >* sampler = 0,
					Gaussian_cdf_kestimator* cdf_estimator = 0,
					Progress_notifier* notifier = 0,
					bool execute_cross_validation = false)
	{
		this->size_rev = size_rev;
		this->n_rev = n_rev;
		this->n_threads = n_threads;
		this->id_thread = id_thread;
		this->fields_path = fields_path; 
		this->fields = fields; 
		this->prop = prop; 
		this->neighborhood_ = neighborhood_;
		this->ccdf = ccdf; 
		this->marginal = marginal; 
		this->sampler = sampler; 
		this->cdf_estimator = cdf_estimator; 
		this->notifier = notifier;
		this->execute_cross_validation = execute_cross_validation;
	}

    void setSampler(Monte_carlo_sampler_t< Mersenne_generator >* sampler_)
	{
		this->sampler = sampler_;
	}

    Monte_carlo_sampler_t< Mersenne_generator >* getSampler()
	{
		return this->sampler;
	}

	void run();

	int getStatus() { return status; }

private:
	int size_rev;
	int n_rev;
	int n_threads;
	int id_thread;
    std::vector<DPG_Grid_path* >* fields_path;
	std::vector<int>* fields;
	Grid_continuous_property* prop;
	SmartPtr<Neighborhood>* neighborhood_;
	Gaussian_cdf* ccdf;
	Gaussian_cdf* marginal;
	Gaussian_cdf_kestimator* cdf_estimator;
    Monte_carlo_sampler_t< Mersenne_generator >* sampler;
	Progress_notifier* notifier;
	int status;
	bool execute_cross_validation;
};


#endif 
