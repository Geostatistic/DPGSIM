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

#include <geostat/library_geostat_init.h>
#include <geostat/parameters_handler_impl.h>

#include <utils/manager_repository.h>
#include <utils/gstl_messages.h>
#include <geostat/two_point_statistics_initialization.h>
#include <geostat/two_point_statistics_local_initialization.h>
#include <geostat/two_point_lmc_initialization.h>
#include <math/continuous_distribution.h>
#include <math/non_parametric_distribution.h>

#include <GsTL/utils/smartptr.h>
#include <utils/gstl_messages.h>
#include <utils/manager_repository.h>


#include "common.h"

#include "dpgsim.h"


extern "C" PLUGINS_LPM_UFRGS_DECL
int plugin_init() 
{
	// ALGORITHMS
    GsTLlog << "\n\n registering DPGSGSIM Algorithms" << "\n";
	SmartPtr<Named_interface> ni = Root::instance()->interface(geostatAlgo_manager);
	Manager* dir = dynamic_cast<Manager*>(ni.raw_ptr());
	if( !dir ) {
		GsTLlog << "Directory " << geostatAlgo_manager << " does not exist \n";
		return 1;
	}
    dir->factory(DPGSgsim().name(), DPGSgsim::create_new_interface);
	

	return 0;
}
