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


#include "dpgsim.h"

int  DPGSgsim::execute_sgsim(int id_thread_real, int n_threads_real, Progress_notifier* notifier)
{
	QTime clock_;
	QString out_buf;
	QTextStream out(&out_buf);
	clock_.start();

	// In sequential gaussian simulation, the marginal is a Gaussian cdf, with mean 0 and variance 1.
	Gaussian_cdf marginal(0.0, 1.0);
	Gaussian_cdf ccdf;

	// work on the fine grid
	if (simul_grid_[id_thread_real]) {
		if (dynamic_cast<RGrid*>(simul_grid_[id_thread_real])) {
			RGrid* rgrid = dynamic_cast<RGrid*>(simul_grid_[id_thread_real]);
			rgrid->set_level(1);
		}
	}

	// set up the cdf-estimator
	Gaussian_cdf_kestimator cdf_estimator(&covar_[id_thread_real], &kriging_system_[id_thread_real], kriging_estimator_[id_thread_real]);

	// set up the sampler
	std::default_random_engine gen_ini;

	gen_ini.seed(this->seed_);

	if (id_thread_real > 0) {
		int r = 0;

		for (int i = 0; i < id_thread_real; ++i) {
			r = gen_ini() % 26666667 + 2954951;
		}
		
		if (r < 0) r = -r;
		if (r % 2 == 0) ++r;

		gen_ini.seed(r);
	}

	/*
	// Random_number_generator gen;
	// Monte_carlo_sampler_t< Random_number_generator > sampler( gen );
	*/

	bool from_scratch = true;

	// loop on all realizations

	mutex_.lock();
	Grid_continuous_property* prop_stages = 0;

	prop_stages = multireal_stage_property_[id_thread_real]->new_realization();
	
	double Lx = 0, Ly = 0, Lz = 0;
	double Xmin = 0, Xmax = 0, Ymin = 0, Ymax = 0, Zmin = 0, Zmax = 0;

	for (int i = 0; i < simul_grid_[id_thread_real]->size(); ++i) {
		double x = simul_grid_[id_thread_real]->location(i).x();
		double y = simul_grid_[id_thread_real]->location(i).y();
		double z = simul_grid_[id_thread_real]->location(i).z();

		if (i == 0) {
			Xmin = Xmax = x;
			Ymin = Ymax = y;
			Zmin = Zmax = z;
		}
		else {
			if (x < Xmin) Xmin = x;
			if (x > Xmax) Xmax = x;

			if (y < Ymin) Ymin = y;
			if (y > Ymax) Ymax = y;

			if (z < Zmin) Zmin = z;
			if (z > Zmax) Zmax = z;
		}
	}

	Lx = Xmax - Xmin;
	Ly = Ymax - Ymin;
	Lz = Zmax - Zmin;

	if (fabs(Zmax - Zmin) < 1e-20) Lz = 1;
	if (fabs(Ymax - Ymin) < 1e-20) Ly = 1;
	if (fabs(Xmax - Xmin) < 1e-20) Lx = 1;

	double Dx = Lx / this->nIRFX, Dy = Ly / this->nIRFY, Dz = Lz / this->nIRFZ;

	if (fabs(Dx) < 1e-20) Dx = 1;
	if (fabs(Dy) < 1e-20) Dy = 1;
	if (fabs(Dz) < 1e-20) Dz = 1;

	std::vector<std::vector<int> > fields(this->nIRFX * this->nIRFY * this->nIRFZ + 1);
	std::vector<std::vector<int> > stages(8);
	std::vector<int> time_stages(8, 0);

	std::map<int, bool> mark;

	for (int i = 0; i < simul_grid_[id_thread_real]->size(); ++i) {
		double x = simul_grid_[id_thread_real]->location(i).x();
		double y = simul_grid_[id_thread_real]->location(i).y();
		double z = simul_grid_[id_thread_real]->location(i).z();

		int ix = static_cast<int>(floor(fabs(x - Xmin) / Dx));
		int iy = static_cast<int>(floor(fabs(y - Ymin) / Dy));
		int iz = static_cast<int>(floor(fabs(z - Zmin) / Dz));

		if (ix >= this->nIRFX) ix = this->nIRFX - 1;
		if (iy >= this->nIRFY) iy = this->nIRFY - 1;
		if (iz >= this->nIRFZ) iz = this->nIRFZ - 1;

		int id_field = ix * this->nIRFY * this->nIRFZ + iy * this->nIRFZ + iz;
		ix %= 2;
		iy %= 2;
		iz %= 2;
		int id_stage = ix * 4 + iy * 2 + iz;

		fields[id_field].push_back(i);
		if (!mark[id_field]) {
			mark[id_field] = true;
			stages[id_stage].push_back(id_field);
		}

		prop_stages->set_value(id_stage, i);	
	}
	mutex_.unlock();

	//Grid_path path(simul_grid_, simul_grid_->selected_property(), target_grid_region_);
    std::vector<DPG_Grid_path*> fields_path(this->nIRFX * this->nIRFY * this->nIRFZ + 1);


	for (size_t i = 0; i < fields.size(); ++i) {
        fields_path[i] = new DPG_Grid_path(fields[i], simul_grid_[id_thread_real], simul_grid_[id_thread_real]->selected_property(), target_grid_region_[id_thread_real]);

		int r = gen_ini() % 26666667 + 2954951;
		if (r < 0) r = -r;
		if (r % 2 == 0) ++r;

		fields_path[i]->setSeed(r);
	}

	size_t n_threads_max = std::min(this->n_threads_, QThread::idealThreadCount());
    std::vector<DPGSIMProcess*> threads(n_threads_max);
    std::vector<DPGSIMCDFTransform< Continuous_distribution, Gaussian_cdf>*>* threads_cdf;
    if (use_target_hist_) threads_cdf = new std::vector<DPGSIMCDFTransform< Continuous_distribution, Gaussian_cdf>*>(n_threads_max);

	for (size_t i = 0; i < n_threads_max; ++i) {
		int r = gen_ini() % 26666667 + 2954951;
		if (r < 0) r = -r;
		if (r % 2 == 0) ++r;

        Mersenne_generator tgen;
		tgen.seed(r);
        Monte_carlo_sampler_t< Mersenne_generator >* tsampler = new Monte_carlo_sampler_t< Mersenne_generator >(tgen);

        threads[i] = new DPGSIMProcess;
		threads[i]->setSampler(tsampler);

        if (use_target_hist_) (*threads_cdf)[i] = new DPGSIMCDFTransform< Continuous_distribution, Gaussian_cdf>(i, n_threads_max, &marginal, target_cdf_[id_thread_real].raw_ptr());
	}

	for (int nreal = id_thread_real; nreal < nb_of_realizations_; nreal += n_threads_real) {

		if (notifier) {
			notifier->write(("Working on realization " + QString::number(nreal + 1)).toLatin1().data(), 0);
			if (!notifier->notify()) {
				clean(fields_path);
				if (use_target_hist_) clean(threads_cdf);
				return 1;
			}
			while (notifier->is_paused()) {
				notifier->write("Paused", 0);
			}
		}

		Grid_continuous_property* prop = multireal_property_[id_thread_real]->new_realization();
		
		prop->swap_to_memory();
		simul_grid_[id_thread_real]->select_property(prop->name());
		neighborhood_[id_thread_real]->select_property(prop->name());

		// initialize the new realization with the hard data, if that was requested
		if (property_copier_[id_thread_real]) {
			property_copier_[id_thread_real]->copy(harddata_grid_[id_thread_real], harddata_property_[id_thread_real],
				simul_grid_[id_thread_real], prop);
		}

		prop->set_parameters(parameters_);
		
		int status = 0;

		for (int nrev = 0; nrev <= n_Revisitations; ++nrev) {
			for (size_t s = 0; s < 8; ++s) {
				QTime clock_s;
				if (stages[s].size() == 0) continue;

				clock_s.start();

				//Create SGSIM threads

				size_t n_threads = std::min(n_threads_max, stages[s].size());

				for (size_t id_thread = 0; id_thread < n_threads; ++id_thread) {

					threads[id_thread]->setValues(
						n_Revisitations,
						nrev,
						n_threads,
						id_thread,
						&fields_path,
						&stages[s],
						prop,
						&neighborhood_[id_thread_real],
						//&this->my_neigh[id_thread],
						new Gaussian_cdf,//&ccdf,
						new Gaussian_cdf(0.0, 1.0),//&marginal,
						threads[id_thread]->getSampler(),
						new Gaussian_cdf_kestimator(&covar_[id_thread_real], &kriging_system_[id_thread_real], kriging_estimator_[id_thread_real]),//&cdf_estimator,
						notifier,
						this->execute_cross_validation_);
					threads[id_thread]->start();
					//threads[id_thread]->wait();
				}

				for (int id_thread = 0; id_thread < n_threads; ++id_thread) {
					threads[id_thread]->wait();

					if (status != -1) {
						status = threads[id_thread]->getStatus();
					}
				}

				time_stages[s] += clock_s.elapsed();

			}

			if (status == -1) {
				clean(id_thread_real, prop);
				clean(fields_path);
				clean(threads);
				if (use_target_hist_) clean(threads_cdf);
				return 1;
			}
		}



		// back-transform if needed
		if (use_target_hist_) {
			for (size_t i = 0; i < n_threads_max; ++i) {
				(*threads_cdf)[i]->setProp(prop);
				(*threads_cdf)[i]->start();
			}
			for (size_t i = 0; i < n_threads_max; ++i) {
				(*threads_cdf)[i]->wait();
			}
		}
		prop->swap_to_disk();
	}

	mutex_.lock();
	clean(id_thread_real);

	if (notifier) notifier->write("Simulation complete ", 0);
	
	clean(fields_path);
	clean(threads);

	if (use_target_hist_) clean(threads_cdf);

	int t = clock_.elapsed();


	out << "Thread " << id_thread_real << " report\n";
	
	for (size_t i = 0; i < 8; ++i) {
		out << "\t Stage " << (i + 1) << " took " << time_stages[i] << " ms to run\n";
	}
	out << "Number of threads : " << n_threads_max << "\n";
	out << "Process took " << t << "ms to run\n";

	out << "Number of IRF X : " << this->nIRFX << "\n";
	out << "Number of IRF Y : " << this->nIRFY << "\n";
	out << "Number of IRF Z : " << this->nIRFZ << "\n";
	out << "----------------------------------------------------\n";

	report_->append(out_buf);
	
	mutex_.unlock();

	return 0;
}

int DPGSgsim::execute(GsTL_project*, Progress_notifier* notifier)
{
	QTime clock_;
	QString out_buf;
	QTextStream out(&out_buf);
	clock_.start();

    std::vector<DPGExecuteSGSIMProcess*> process;

	for (int id_thread = 0; id_thread < nb_of_parallel_realizations_; ++id_thread) {
        process.push_back(new DPGExecuteSGSIMProcess(this, id_thread, nb_of_parallel_realizations_, notifier));
	}

	for (int id_thread = 0; id_thread < nb_of_parallel_realizations_; ++id_thread) {
		process[id_thread]->start();
	}

	for (int id_thread = 0; id_thread < nb_of_parallel_realizations_; ++id_thread) {
		process[id_thread]->wait();
	}

	for (int id_thread = 0; id_thread < nb_of_parallel_realizations_; ++id_thread) {
		delete process[id_thread];
	}
	
	int t = clock_.elapsed();

	
	out << "Total time : " << t << "ms to run\n";
	report_->append(out_buf);
	set_report_is_open(true);

	return 0;
}

void DPGSIMProcess::run()
{
	for (size_t f = this->id_thread; f < this->fields->size(); f += this->n_threads) {
		if (notifier) { if (!notifier->notify()) break; }
		size_t id_field = fields->at(f);

		if (n_rev == 0) {
			fields_path->at(id_field)->set_property(prop->name());
			fields_path->at(id_field)->randomize();
		}


		// do the simulation
		GsTLInt size_visit = fields_path->at(id_field)->size() / (this->size_rev + 1);
		GsTLInt ini = n_rev * size_visit;
		GsTLInt fin = std::min(fields_path->at(id_field)->size(), (n_rev + 1) * size_visit);

		if (n_rev == this->size_rev) {
			fin = fields_path->at(id_field)->size();
		}

		int status =
			LPM_UFRGS::sequential_simulation(fields_path->at(id_field), ini,
			fin,
			*(neighborhood_->raw_ptr()),
			*ccdf,
			*cdf_estimator,
			*marginal,
			*sampler,
			notifier,
			this->execute_cross_validation
			);
		this->status = status;
		if (notifier) { if (!notifier->notify()) break; }
		if (status == -1) break;
	}
}


bool DPGSgsim::initialize(const Parameters_handler* parameters,
	Error_messages_handler* errors, Progress_notifier* notifier)
{
	OPEN_DEBUG_STREAM("sgsim.dbg");

	nb_of_realizations_ =
		String_Op::to_number<int>(parameters->value("Nb_Realizations.value"));
	nb_of_parallel_realizations_ =
		String_Op::to_number<int>(parameters->value("nb_parallel_real.value"));

	seed_ = String_Op::to_number<int>(parameters->value("Seed.value"));
	
	this->n_Revisitations = String_Op::to_number<int>(parameters->value("nRevisitations.value")) - 1;

	this->dx = String_Op::to_number<double>(parameters->value("dx.value"));
	this->dy = String_Op::to_number<double>(parameters->value("dy.value"));
	this->dz = String_Op::to_number<double>(parameters->value("dz.value"));

	if (String_Op::to_number<int>(parameters->value("enableCustomizeIRF.value"))) {
		this->nIRFX = String_Op::to_number<int>(parameters->value("nIRFX.value"));
		this->nIRFY = String_Op::to_number<int>(parameters->value("nIRFY.value"));
		this->nIRFZ = String_Op::to_number<int>(parameters->value("nIRFZ.value"));
	}
	else {
		this->nIRFX = 1;
		this->nIRFY = 1;
		this->nIRFZ = 1;
	}

	this->execute_cross_validation_ = String_Op::to_number<int>(parameters->value("executeCrossValidation.value"));

	this->n_threads_ = String_Op::to_number<int>(parameters->value("nThreads.value"));

	std::string harddata_grid_name = parameters->value("Hard_Data.grid");

	// Extract the parameters input by the user from the parameter handler

	//-------------
	// The "simulation" grid parameters

	// Extract the parameters as a XML
	//  const Parameters_handler_xml *paraCopy =
	//          static_cast<const Parameters_handler_xml*> (parameters);
	//  this->currentParametersXML = paraCopy->getDoc();
	std::string simul_grid_name = parameters->value("Grid_Name.value");
	this->simul_grid_name_ = simul_grid_name;
	errors->report(simul_grid_name.empty(),
		"Grid_Name", "No grid selected");
	std::string property_name = parameters->value("Property_Name.value");
	this->property_name_ = property_name;
	errors->report(property_name.empty(),
		"Property_Name", "No property name specified");

	// Get the simulation grid from the grid manager
	if (simul_grid_name.empty()) return false;

	simul_grid_.resize(nb_of_parallel_realizations_);
	harddata_grid_.resize(nb_of_parallel_realizations_);
	target_cdf_.resize(nb_of_parallel_realizations_);
	covar_.resize(nb_of_parallel_realizations_);
	kriging_system_.resize(nb_of_parallel_realizations_);
	property_copier_.resize(nb_of_parallel_realizations_);

	for (int id_thread = 0; id_thread < nb_of_parallel_realizations_; ++id_thread) {
		
		bool ok = geostat_utils::create(simul_grid_[id_thread], simul_grid_name,
			"Grid_Name", errors);
		if (!ok) return false;

		// create  a multi-realization property
		std::string suffix = QString::number(id_thread).toStdString();
		multireal_property_.push_back(simul_grid_[id_thread]->add_multi_realization_property(property_name + "_t" + suffix));

		multireal_stage_property_.push_back(simul_grid_[id_thread]->add_multi_realization_property("stages_t" + suffix));

		target_grid_region_.push_back(simul_grid_[id_thread]->region(parameters->value("Grid_Name.region")));


		//-------------
		// The hard data parameters
		if (!harddata_grid_name.empty()) {
			std::string hdata_prop_name = parameters->value("Hard_Data.property");
			errors->report(hdata_prop_name.empty(),
				"Hard_Data", "No property name specified");

			// Get the harddata grid from the grid manager
			bool ok = geostat_utils::create(harddata_grid_[id_thread], harddata_grid_name,
				"Hard_Data", errors);
			if (!ok) return false;

			harddata_property_.push_back(harddata_grid_[id_thread]->property(hdata_prop_name));

			if (!harddata_property_[id_thread]) {
				std::ostringstream error_stream;
				error_stream << harddata_grid_name
					<< " does not have a property called "
					<< hdata_prop_name;
				errors->report("Hard_Data", error_stream.str());
				return false;
			}

		}
		if (dynamic_cast<Point_set*>(harddata_grid_[id_thread])) {
			harddata_grid_[id_thread]->set_coordinate_mapper(simul_grid_[id_thread]->coordinate_mapper());
		}

		std::string hd_region_name = parameters->value("Hard_Data.region");
		if (harddata_grid_[id_thread]) {
			hd_grid_region_.push_back(harddata_grid_[id_thread]->region(hd_region_name));
		}


		// hard data assignement and transform is only needed if we have a valid
		// hard data grid and property.  We always assigne the data if it belongs
		// the same grid
		bool  assign_harddata =
			String_Op::to_number<bool>(parameters->value("Assign_Hard_Data.value"));
		if (harddata_grid_[id_thread] == NULL) assign_harddata = false;
		else if (harddata_grid_[id_thread] == simul_grid_[id_thread]) assign_harddata = true;

		if (assign_harddata) {
			property_copier_[id_thread] = 
				Property_copier_factory::get_copier(harddata_grid_[id_thread], simul_grid_[id_thread]);
			if (!property_copier_[id_thread]) {
				std::ostringstream message;
				message << "It is currently not possible to copy a property from a "
					<< harddata_grid_[id_thread]->classname() << " to a "
					<< simul_grid_[id_thread]->classname();
				errors->report(!property_copier_[id_thread], "Assign_Hard_Data", message.str());
				return false;
			}
		}
		else {
			property_copier_[id_thread] = 0;
		}


		//-------------
		// Target histogram

		use_target_hist_ =
			String_Op::to_number<bool>(parameters->value("Use_Target_Histogram.value"));

		if (use_target_hist_) {
			bool ok = distribution_utils::get_continuous_cdf(target_cdf_[id_thread], parameters, errors, "nonParamCdf");
			if (!ok) return false;

			if (harddata_property_[id_thread]) {
				harddata_property_[id_thread] =
					distribution_utils::gaussian_transform_property(harddata_property_[id_thread], target_cdf_[id_thread].raw_ptr(), harddata_grid_[id_thread], hd_grid_region_[id_thread]);
				if (!harddata_property_[id_thread]) return false;

				clear_temp_properties_ = true;
				harddata_grid_[id_thread]->select_property(harddata_property_[id_thread]->name());
			}
		}

		//-------------
		// Variogram (covariance) initialization

		bool init_cov_ok =
			geostat_utils::initialize_two_point_nested_structure(&covar_[id_thread], "covariance_input",
			parameters, errors, simul_grid_[id_thread]);
		if (!init_cov_ok) return false;

		//-------------
		// Set up the search neighborhood

		int max_neigh =
			String_Op::to_number<int>(parameters->value("Max_Conditioning_Data.value"));
		int max_neigh_simul =
			String_Op::to_number<int>(parameters->value("Max_Conditioning_Simul_Data.value"));

		GsTLTriplet ranges;
		GsTLTriplet angles;
		bool extract_ok =
			geostat_utils::extract_ellipsoid_definition(ranges, angles,
			"Search_Ellipsoid.value",
			parameters, errors);
		if (!extract_ok) return false;


		// If the hard data are not "relocated" on the simulation grid,
		// use a "combined neighborhood", otherwise use a single
		// neighborhood
		// The octant search is only used on the hard data
		// The max size is set for each neighborhood not for the combined


		//this->my_neigh.clear();

		if (!harddata_grid_[id_thread] || assign_harddata) {
			/*
			neighborhood_ = SmartPtr<Neighborhood>(
			simul_grid_->neighborhood( ranges, angles, &covar_ ) );
			*/
			//for (int i = 0; i < this->n_threads_; ++i) {
			neighborhood_.push_back(SmartPtr<Neighborhood>(
				simul_grid_[id_thread]->neighborhood(ranges, angles)));
			neighborhood_[id_thread]->max_size(max_neigh_simul);
			//this->my_neigh.push_back(neighborhood_);
			//}
			//   geostat_utils::set_advanced_search(neighborhood_,
			//                      "AdvancedSearch", parameters, errors);


		}
		else {
			//Neighborhood* simul_neigh  = simul_grid_->neighborhood( ranges, angles, &covar_ );
			Neighborhood* simul_neigh = simul_grid_[id_thread]->neighborhood(ranges, angles);

			simul_neigh->max_size(max_neigh_simul);
			//   geostat_utils::set_advanced_search(simul_neigh,
			//                      "AdvancedSearch", parameters, errors);


			harddata_grid_[id_thread]->select_property(harddata_property_[id_thread]->name());

			Neighborhood* harddata_neigh;
			if (dynamic_cast<Point_set*>(harddata_grid_[id_thread])) {
				harddata_neigh =
					harddata_grid_[id_thread]->neighborhood(ranges, angles, 0, true, hd_grid_region_[id_thread]);
			}
			else {
				harddata_neigh =
					harddata_grid_[id_thread]->neighborhood(ranges, angles, 0, false, hd_grid_region_[id_thread]);
			}
			/*
			if( dynamic_cast<Point_set*>(harddata_grid_) ) {
			harddata_neigh =
			harddata_grid_->neighborhood( ranges, angles, &covar_, true, hd_grid_region_ );
			}
			else {
			harddata_neigh =
			harddata_grid_->neighborhood( ranges, angles, &covar_, false, hd_grid_region_ );
			}
			*/
			harddata_neigh->max_size(max_neigh);
			geostat_utils::set_advanced_search(harddata_neigh,
				"AdvancedSearch", parameters, errors);
			//  harddata_neigh->select_property( harddata_property_->name() );

			//for (int i = 0; i < this->n_threads_; ++i) {
			neighborhood_.push_back(
				SmartPtr<Neighborhood>(new Combined_neighborhood(harddata_neigh,
				simul_neigh, 0)));

			//this->my_neigh.push_back(neighborhood_);
			//}
			/*
			neighborhood_ =
			SmartPtr<Neighborhood>( new Combined_neighborhood( harddata_neigh,
			simul_neigh, &covar_) );
			*/

		}

		if (execute_cross_validation_) neighborhood_[id_thread]->includes_center(false);

		// neighborhood_->max_size( max_neigh ); // The constraint is on the individual neighborhood
		//geostat_utils::set_advanced_search(neighborhood_,
		//                    "AdvancedSearch", parameters, errors);


		kriging_system_[id_thread].kriging_constraint_is(SK_constraint());
		kriging_system_[id_thread].lhs_covariance_is(covar_[id_thread]);
		kriging_system_[id_thread].rhs_covariance_is(covar_[id_thread]);
		kriging_estimator_.push_back(new SK_estimator());
		((SK_estimator *)kriging_estimator_[id_thread])->kriging_mean_is(0.0);


		//-----------------
		// The kriging constraints and combiner
		/*
		geostat_utils::KrigTagMap tags_map;
		tags_map[ geostat_utils::KT  ] = "Trend.value";
		tags_map[ geostat_utils::LVM ] = "Local_Mean_Property.value";

		geostat_utils::KrigDefaultsMap defaults;
		defaults[ geostat_utils::SK ] = "0.0";

		std::string kriging_base_type = "Kriging_Type";
		geostat_utils::initialize_kriging_system_( kriging_base_type, kriging_system_, kriging_estimator_,
		parameters, errors,
		simul_grid_, NULL, defaults );
		*/
	}

	//----------------
	// Report errors if any

	if (!errors->empty()) {
		
		for (int i = 0; i < nb_of_parallel_realizations_; ++i) {
			clean(i);
		}

		return false;
	}

	this->extract_parameters(parameters);


	// total_steps = simul_grid_->size() * (nb_of_realizations_);
	//frequency = std::max( total_steps / 20, 1 );
	//progress_notifier =
	// utils::create_notifier( "Running Sgsim",
	//		    total_steps, frequency );

	if (notifier) {
		int active_cells = target_grid_region_[0] ? target_grid_region_[0]->active_size() : simul_grid_[0]->size();
		total_steps = nb_of_realizations_*active_cells;
		notifier->total_steps(total_steps);
		notifier->frequency(20);
	}

	return true;
}



void DPGSgsim::clean(int id_thread, Grid_continuous_property* prop)
{
	if (prop)
		if (simul_grid_.size() > 0)
		simul_grid_[id_thread]->remove_property(prop->name());

	if (harddata_grid_.size() > 0 && harddata_property_.size() > 0)
	if (clear_temp_properties_ && harddata_property_[id_thread] && harddata_grid_[id_thread]) {
		harddata_grid_[id_thread]->remove_property(harddata_property_[id_thread]->name());
		harddata_property_[id_thread] = 0;
	}
}

Named_interface* DPGSgsim::create_new_interface(std::string&)
{
    return new DPGSgsim;
}


DPGSgsim::DPGSgsim()
{
	use_target_hist_ = false;
	clear_temp_properties_ = false;
    report_ = new DPGSGSIM_report(0);
	report_->setSelf(report_);
	report_is_open_ = false;
}


DPGSgsim::~DPGSgsim()
{

	for (int i = 0; i < nb_of_parallel_realizations_; ++i) {
		if (kriging_estimator_.size() > 0) if (kriging_estimator_[i] != 0) delete kriging_estimator_[i];

		if (harddata_grid_.size() > 0)
			if (harddata_grid_[i] != 0 && dynamic_cast<Point_set*>(harddata_grid_[i])) {
				harddata_grid_[i]->set_coordinate_mapper(0);
			}

		clean(i);
	}

	if (!report_is_open_) delete report_;
}

void DPGSgsim::clean(std::vector<DPG_Grid_path* >& path)
{
	for (size_t i = 0; i < path.size(); ++i) delete path[i];
}

void DPGSgsim::clean(std::vector<DPGSIMProcess*>& threads)
{
	for (int id_thread = 0; id_thread < threads.size(); ++id_thread) {
		delete threads[id_thread];
	}
}

void DPGSgsim::clean(std::vector<DPGSIMCDFTransform<Continuous_distribution, Gaussian_cdf>*>* threads)
{
	for (int id_thread = 0; id_thread < threads->size(); ++id_thread) {
		delete threads->at(id_thread);
	}
	delete threads;
}
