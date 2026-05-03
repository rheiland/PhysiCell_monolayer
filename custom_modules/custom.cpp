/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include <algorithm>    // for std::remove

#include "./custom.h"

double time_90pct;  // time at which the outer cells reach 90% max width
double xpos_90pct;
bool cells11_flag = true;

double rest_length_factor = 1.0;

double tmin=0;
double tmax= 1;

void create_cell_types( void )
{
    if (parameters.doubles.find_index("rest_length_factor") != -1)
	{
        rest_length_factor = parameters.doubles("rest_length_factor");
	}
    else
    {
        std::cout << "Error: rest_length_factor needs to be defined in user params\n" << std::endl;
        std::exit(-1);
    }
	
	//   Put any modifications to default cell definition here if you 
	//   want to have "inherited" by other cell types. 
	   
	//   This is a good place to set default functions. 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	// cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
    cell_defaults.functions.cell_division_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	//   This parses the cell definitions in the XML config file. 
	initialize_cell_definitions_from_pugixml(); 

	//   This builds the map of cell definitions and summarizes the setup. 
	build_cell_definitions_maps(); 

	//   This intializes cell signal and response dictionaries 
	setup_signal_behavior_dictionaries(); 	

    //   Cell rule definitions 
	setup_cell_rules(); 

	//   Put any modifications to individual cell definitions here. 
	//   This is a good place to set custom functions. 
	
	// cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_cell_rule; 
    cell_defaults.functions.update_velocity = custom_update_cell_velocity;
	// cell_defaults.functions.update_velocity = standard_update_cell_velocity;

    xpos_90pct = (0.9 * cell_defaults.phenotype.geometry.radius * 2 * 10.0) / 2.0;
    std::cout << "-------- xpos_90pct= " << xpos_90pct << std::endl;
	
	//   This builds the map of cell definitions and summarizes the setup. 
	display_cell_definitions( std::cout ); 
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	initialize_microenvironment(); 	
	return; 
}

void setup_tissue( void )
{
    // First, decide if we're doing the 11 or 21 cell relaxation model
    cells11_flag = parameters.bools("cells11");

	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 
	

	// create some of each type of cell 
	
	Cell* pC;
	
	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
		{
			std::vector<double> position = {0,0,0}; 
			position[0] = Xmin + UniformRandom()*Xrange; 
			position[1] = Ymin + UniformRandom()*Yrange; 
			position[2] = Zmin + UniformRandom()*Zrange; 
			
			pC = create_cell( *pCD ); 
			pC->assign_position( position );
		}
	}
	std::cout << std::endl; 
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml();
	set_parameters_from_distributions();
	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }


// do every mechanics dt
void custom_update_cell_velocity( Cell* pCell, Phenotype& phenotype, double dt )
{
    pCell->custom_data["cell_ID"] = pCell->ID;
    static double effective_repulsion = phenotype.mechanics.cell_cell_repulsion_strength;  // e.g., 10

    // if (pCell->ID == 1 && PhysiCell_globals.current_time > tmin && PhysiCell_globals.current_time < tmax )
    // {
    //     std::cout << "------ nearby_cells\n";
    //     for (const auto& nbr : pCell->nearby_cells ())   // nearby_cells or nearby_interacting_cells
    //         std::cout << __FUNCTION__ << ": t=" <<PhysiCell_globals.current_time <<",         ID= " << pCell->ID << " has nbr->ID = " << nbr->ID << std::endl;

    //     std::cout << "------ nearby_interacting_cells\n";
    //     for (const auto& nbr : pCell->nearby_interacting_cells())   
    //         std::cout << __FUNCTION__ << ": t=" <<PhysiCell_globals.current_time <<",         ID= " << pCell->ID << " has nbr->ID = " << nbr->ID << std::endl;
    //     // std::cout << " ~~nearby_interacting_cells are " << pCell->nearby_interacting_cells() << std::endl;
    // }

    // ----------  Step 1) determine who are my nbrs  ------------
	pCell->state.neighbors.clear();
    for (const auto& nbr : pCell->nearby_interacting_cells())  // choice: nearby_cells or nearby_interacting_cells
    {
        if (pCell->ID == nbr->ID)  // let's NOT count me as a nbr of me (as is done in core)
            continue; 

        double rest_length = (pCell->phenotype.geometry.radius + nbr->phenotype.geometry.radius) * rest_length_factor;
        double dx = nbr->position[0] - pCell->position[0];
        double dy = nbr->position[1] - pCell->position[1];  // should be 0.0
        double distance = std::sqrt(dx*dx + dy*dy);

        // if (pCell->ID == 1 && PhysiCell_globals.current_time > tmin && PhysiCell_globals.current_time < tmax )
        // {
        //     std::cout << "   t=" <<PhysiCell_globals.current_time <<" :  " << pCell->ID << " ->"<< nbr->ID << ": dist=" << distance <<" , rest_length="<<rest_length << std::endl;
        // }

        if (distance <= rest_length)
        {
            pCell->state.neighbors.push_back(nbr);
        }
    }

    // if (pCell->ID == 1 && PhysiCell_globals.current_time > tmin && PhysiCell_globals.current_time < tmax )
    // {
    //     std::cout << "          -- state.neighbors  after constraining" << std::endl;
    //     for (const auto& nbr : pCell->state.neighbors)
    //         std::cout << ",         ID= " << pCell->ID << " has nbr->ID = " << nbr->ID << std::endl;
    // }

    // if (pCell->ID == 1 && PhysiCell_globals.current_time > tmin && PhysiCell_globals.current_time < tmax )
    // {
    //     std::cout << "    after 3) ~~~~~~~~~~~~~   final cull:  pCell->state.neighbors:" << std::endl;
    //     for (const auto& nbr : pCell->state.neighbors)
    //         std::cout << "             ID= " << pCell->ID << " has nbr->ID = " << nbr->ID << std::endl;
    // }


    pCell->custom_data["num_nbrs"] = 0;

    // ----------  Step 2) compute my velocity using my (corrected) nbrs ------------
    pCell->velocity = {0.0, 0.0, 0.0};
    for (const auto& pNeighbor : pCell->state.neighbors)
    {
        // double rest_length = (pCell->phenotype.geometry.radius + pNeighbor->phenotype.geometry.radius) * rest_length_factor;

        // NOTE: we do NOT use rest_length_factor here
        double rest_length = (pCell->phenotype.geometry.radius + pNeighbor->phenotype.geometry.radius);

        // if (pCell->ID <= 1)
        // {
        //     std::cout << "---   nbr ID= " << pNeighbor->ID << ", rest_length= " << rest_length << std::endl;
        // }

        // Vector from A (pCell) toward B (pNeighbor)
        double dx = pNeighbor->position[0] - pCell->position[0];

        // if (PhysiCell_globals.current_time < 0.05 && (pCell->ID == 0) && (pNeighbor->ID==1))
        // {
        //     std::cout << "------ " << PhysiCell_globals.current_time << ", rest_length= " << rest_length << ", dx= "<<dx<< std::endl;   
        //     // ------ 0, rest_length= 10, dx= 5
        //     // ------ 0.01, rest_length= 10, dx= 4.99625
        //     // ------ 0.02, rest_length= 10, dx= 4.99374
        // }
        double dy = pNeighbor->position[1] - pCell->position[1];  // should be 0.0
        // double dy = 0.;
        double dz = 0.; // pNeighbor->position[2] - pCell->position[2];
        //jdouble distance = std::sqrt( dx*dx + dy*dy + dz*dz );
        // double distance = std::sqrt( dx*dx + dy*dy);
        double distance = std::abs(dx);   // keep it simple

        if( distance < 1e-12 ) continue;  // avoid division by zero

        double overlap = distance - rest_length;  // negative when cells overlap

        if( overlap >= 0.0 )
        {
            continue;
        }
        else
        {
            pCell->custom_data["num_nbrs"] += 1;   // now performed in custom_cell_rule[_slow], after update_pos

            // if ((pCell->custom_data["num_nbrs"] == 1) && pCell->ID == 1)
            // {
            //     // std::cout << "---   cell ID= " << pCell->ID << " has just 1 nbr: ID = " << pNeighbor->ID << ", t="<<PhysiCell_globals.current_time << std::endl;
            //     nbr_ID = pNeighbor->ID;
            //     t_saved = PhysiCell_globals.current_time;
            //     // true_nbrs.push_back();

            //     // std::cout << "---   cell ID= " << pCell->ID << " had nbrs = " << pCell->state.neighbors <<", latest has ID = " << pNeighbor->ID << std::endl;
            //     // if (pNeighbor->ID != 0)
            //         // std::cout << "---   cell ID= " << pCell->ID << " has nbr ID = " << pNeighbor->ID << std::endl;
            // }

            //------- PhysiCell standard velocity calculation
            // double temp_r = -distance; // -d

            // double temp_r = -distance; // -d
            // temp_r /= rest_length; // -d/R
            // temp_r += 1.0; // 1-d/R
            // temp_r *= temp_r; // (1-d/R)^2 

            // double effective_repulsion = sqrt( phenotype.mechanics.cell_cell_repulsion_strength * other_agent->phenotype.mechanics.cell_cell_repulsion_strength ); 
            // static double effective_repulsion = phenotype.mechanics.cell_cell_repulsion_strength * other_agent->phenotype.mechanics.cell_cell_repulsion_strength ); 

            // temp_r *= effective_repulsion; 

            double temp = 1.0 - (distance / rest_length);   // normalized overlap fraction
            // double magnitude = stiffness * temp * temp;
            double magnitude = effective_repulsion * temp * temp;


            // --- cf. PhysiCell /core: void Cell::add_potentials(Cell* other_agent)
            //	velocity[i] += displacement[i] * temp_r; 
    	    // axpy( &velocity , temp_r , displacement );    # core/PhysiCell_cell.cpp

            // quadratic repulsion
            // pCell->velocity[0] += magnitude * dx / distance;
            // pCell->velocity[1] += magnitude * dy / distance;
            // pCell->velocity[2] += magnitude * dz / distance;
            // pCell->velocity[2] = 0.0;

            // PhysiCell
            // temp_r /= distance;
            // for( int i = 0 ; i < 3 ; i++ ) 
            // {
            //	velocity[i] += displacement[i] * temp_r; 
            // }
            // axpy( &velocity , temp_r , displacement ); 
            // pCell->velocity[0] += temp_r * dx / distance;
            // pCell->velocity[1] += temp_r * dy / distance;
            // pCell->velocity[0] += temp_r * dx;
            // pCell->velocity[0] += magnitude * dx / distance;
            pCell->velocity[0] -= magnitude * dx / distance;     // rwh: yipeee: negate for repulsion

            // pCell->velocity[1] += temp_r * dy;
            // pCell->velocity[1] = 0.0;
        }
    }

    // if ((pCell->custom_data["num_nbrs"] == 1) && pCell->ID == 1)
    // {
    //     std::cout << "--- FINAL:   cell ID= " << pCell->ID << " has 1 nbr: ID = " << nbr_ID << ", t_saved="<< t_saved << std::endl;
    //     // std::cout << "---   cell ID= " << pCell->ID << " had nbrs = " << pCell->state.neighbors <<", latest has ID = " << pNeighbor->ID << std::endl;
    //     // if (pNeighbor->ID != 0)
    //         // std::cout << "---   cell ID= " << pCell->ID << " has nbr ID = " << pNeighbor->ID << std::endl;
    // }
}

// called every dt_mech; pC->functions.custom_cell_rule( pC,pC->phenotype,time_since_last_mechanics );
void custom_cell_rule( Cell* pCell, Phenotype& phenotype , double dt )
{ 
    static bool reached_90 = false;
    std::stringstream ss;

    pCell->custom_data["cell_ID"] = pCell->ID;
    pCell->custom_data["num_nbrs"] = pCell->state.neighbors.size();
    pCell->custom_data["vel_mag"] = std::sqrt( pCell->previous_velocity[0]*pCell->previous_velocity[0] + pCell->previous_velocity[1]*pCell->previous_velocity[1] );


    // if (pCell->ID == 10 && pCell->position[0] >= 90.0)  // 90%, 90pct
    if (pCell->ID == 10)
    {
        // std::cout << "------- " << __FUNCTION__ << ":  ID= " << pCell->ID <<":  x= " << pCell->position[0] << std::endl;
        // full width: -5*d to 5*d, e.g. d=10 -> 50-(-50)=100; 90% is 90; divide by 2 = 45 (x-position)
        // full width: -5*d to 5*d, e.g. d=10 -> 50-(-50)=100; 90% is 90; divide by 2 = 45 (x-position)
        // if (pCell->position[0] >= 45.0)  // for symmetric test (hard-coded for diam=10)
        if (pCell->position[0] >= xpos_90pct)  // for symmetric test (hard-coded for diam=10)
        {
            if (!reached_90)
            {
                std::cout <<"---- "<< __FUNCTION__ << ": cell radius = " << phenotype.geometry.radius << std::endl;
                time_90pct = PhysiCell_globals.current_time;
                std::cout <<"---- "<< __FUNCTION__ << ": Width reached 90% , t= " << time_90pct << std::endl;

                std::ofstream outFile;
                // sprintf( filename , "%s/output%08u" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index );
                ss << PhysiCell_settings.folder.c_str()  << "/time_90pct.txt" ;
                std::string out_file = ss.str();
                std::cout <<"---- "<< __FUNCTION__ << ": time_pct out_file=" << out_file << std::endl;
                outFile.open(out_file);   // read by  ../analysis/plot_11cells_sweep.py  for plotting results
                outFile << time_90pct;
                outFile.close();

                reached_90 = true;
                // std::exit();
            }
            else
            {
                // if (PhysiCell_globals.current_time / time_90pct > 10.0)
                if (cells11_flag && PhysiCell_globals.current_time / time_90pct > 10.0)
                {
                    std::cout <<"---- "<< __FUNCTION__ << "  calibrated T > 10.  Exit simulation!" << std::endl;
                    std::exit(-1);
                }

            }
        }
    }
} 