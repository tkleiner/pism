/* Copyright (C) 2019, 2020 Thomas Kleiner, Johannes Sutter
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */


#include "ISMIP6.hh"
#include "ISMIP6Physics.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/Vars.hh"
#include <gsl/gsl_math.h>

#include "pism/coupler/util/options.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/util/Mask.hh"
#include "pism/util/Time.hh"
#include "pism/util/iceModelVec.hh"
#include "pism/util/pism_options.hh"

#include <algorithm> // std::max
#include <cmath>     // std::abs
#include <string>
#include <vector>

namespace pism {
namespace ocean {

//
// The Example in doc/sphinx/technical/code/Example.cc|hh is a very good starting point
//
ISMIP6::ISMIP6(IceGrid::ConstPtr g) : CompleteOceanModel(g, nullptr), m_physics(new ISMIP6Physics(g)) {

  // option call for activating 3d ocean forcing for ismip6 -ocean ismip6,...

  ForcingOptions opt(*m_grid->ctx(), "ocean.ismip6");

  if (opt.period > 0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "ISMIP6::ISMIP6(): period = %d not implemented yet.",
                                  opt.period);
  }

  if (opt.reference_time > 0.0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "ISMIP6::ISMIP6(): reference_year = %f not implemented yet.",
                                  opt.reference_time);
  }

  // Thermal forcing is a 4d field here TF(x,y,z,t), but
  // grid_info contains only the coordinates std::vector<double> x, y, z;
  // PIO file(m_grid->com, "guess_mode", opt.filename, PISM_READONLY);
  File file(m_grid->com, opt.filename, PISM_GUESS, PISM_READONLY);


  grid_info info(file, "thermal_forcing", m_sys, m_grid->registration());
  file.close();
  m_log->message(2,
                 "  - thermal_forcing found in '%s'\n"
                 "    with dimensions x: %d, y: %d, z: %d, t: %d \n",
                 opt.filename.c_str(), info.x_len, info.y_len, info.z_len, info.t_len);

  if (info.t_len < 2) {
    m_time_dependent = false;
  } else {
    m_time_dependent = true; // this should be the case for most applications
  }

  // PISM expects all coordinate variables in input
  // files (and internally in PISM) to be strictly increasing. Thus we start with the most negative values
  // at the bottom of the ocean (e.g. -5500m ) and go up to 0m (or slightly below depending on the data).
  // z_min -> bottom, z_max -> ocean surface, less or equal 0 m
  m_Moz = info.z_len;
  m_Loz = -info.z_min + info.z_max;
  m_log->message(2, "  - ocean z range: %f m to %f m  (Loz = %f,  Moz = %d)\n", info.z_min, info.z_max, m_Loz, m_Moz);


  // todo: think about z_bounds
  // ...
  //  -180, -120,
  //  -120, -60,
  //  -60, -0 ;   -> thus we have the surface at 0 m instead of the reported range of -30m from the z values

  // use the z levels from input
  m_z_ocean_surface_layer    = info.z_max; // should be 0 ore something slightly below the surface (-6m in COSMOS data)
  m_z_ocean_bottom_layer     = info.z_min; // total water column
  std::vector<double> zocean = info.z;


  // create instance of IceModelVec3Custom (todo: stable1.1 version not working in dev)
  {
    std::map<std::string, std::string> attrs;
    attrs["units"]     = "m";
    attrs["long_name"] = "Z-coordinate in ocean";
    attrs["axis"]      = "Z";
    attrs["positive"]  = "up";
    m_thermal_forcing.create(m_grid, "thermal_forcing", "zOcean", zocean, attrs);
  }

  // I prefer degC as in NCL and CDO examples
  m_thermal_forcing.set_attrs("climate_forcing",
                              "Thermal Forcing (in situ Temperature minus freezing temperature)",
                              "degC", "", "", 0);

  // this is not really needed and could be derived from  shelf_base_mass_flux in the post processing
  m_forcing.create(m_grid, "thermal_forcing_at_depth", WITHOUT_GHOSTS);
  m_forcing.set_attrs("diagnostic",
                      "ISMIP6 thermal forcing at ice shelf base",
                      "degC", "", "", 0);

  // this output is not really needed and could be derived from
  // shelf_base_mass_flux in the post-processing
  m_basal_melt_rate.create(m_grid, "ismip6_basal_melt_rate", WITHOUT_GHOSTS);
  m_basal_melt_rate.set_attrs("model_state",
                              "ISMIP6 sub-shelf melt rate", "m s-1",
                              "m year-1", "", 0);

  // dummy time series data to get excess to the times stored in the forcing file
  if (m_time_dependent) {
    m_tser     = nullptr;
    auto tname = m_config->get_string("time.dimension_name");
    m_tser     = new Timeseries(*m_grid, tname, tname);
    m_tser->dimension().set_string("units", m_grid->ctx()->time()->units_string());
    // IceGrid.cc grid_info::report() time units are expected to be in seconds
    m_tser->variable().set_string("units", "seconds");            // read will not work without proper units
    m_tser->variable().set_string("glaciological_units", "days"); // read will not work without proper units

    // this is not working, since this uses udunits for the conversion instead of the 365_day calendar
    //  m_tser->variable().set_string("units", "years"); // read will not work without proper units
    //  m_tser->variable().set_string("glaciological_units", "years"); // read will not work without proper units

    m_log->message(2, "  - forcing from '%s' is time dependent\n", opt.filename.c_str());
  } else {
    m_log->message(2, "  - forcing from '%s' is time independent\n", opt.filename.c_str());
  }

  // for debug dumps
  m_first_update = true;

  // as in BedDef.cc
  m_t_last = GSL_NAN;
}

ISMIP6::~ISMIP6() {
  // empty
}

void ISMIP6::init_impl(const Geometry &geometry) {

  (void)geometry;
  ForcingOptions opt(*m_grid->ctx(), "ocean.ismip6");

  m_log->message(2, "* Initializing the ISMIP6 model for the ocean ...\n");

  // as in BedDef.cc init with the first year of simulation
  m_t_last = m_grid->ctx()->time()->start();

  // read basin file and deltaT file
  m_physics->init_impl();

  // init with size and value
  m_basin_count = m_physics->get_basin_count();
  m_basin_area.resize(m_basin_count, 0.0);
  m_basin_average_forcing.resize(m_basin_count, 0.0);

  m_thermal_forcing.set(0.0);

  // READ TIMES FROM TF(x,y,z,t) file
  if (m_time_dependent) {
    m_log->message(2, "ISMIP6::init_impl(...): reading times only from file '%s'\n", opt.filename.c_str());
    //PIO file(m_grid->com, "guess_mode", opt.filename, PISM_READONLY);
    File file(m_grid->com, opt.filename, PISM_GUESS, PISM_READONLY);
    m_tser->read(file, *m_grid->ctx()->time(), *m_grid->ctx()->log());
    file.close();
#if (PISM_DEBUG == 1)
    // report all time steps from TF file
    std::vector<double> times = m_tser->times();
    std::vector<double>::iterator itr;
    const Time &time = *m_grid->ctx()->time();
    for (itr = times.begin(); itr < times.end(); ++itr) {
      m_log->message(2, "ISMIP6::init_impl(...): time '%.0f s' -> day '%.0f' -> end of year '%f' \n", *itr,
                     time.convert_time_interval(*itr, "days"), time.convert_time_interval(*itr, "years"));
    }
#endif
  } else {
    // READING in the data directly
    m_log->message(2, "ISMIP6::init_impl(...): regrid thermal forcing from file '%s'\n", opt.filename.c_str());
    m_thermal_forcing.regrid(opt.filename, CRITICAL); //TODO: check if CRITICAL_FILL_MISSING is beneficial
  }
}

// Derived from pico.cc|hh
void ISMIP6::define_model_state_impl(const File &output) const {

  //m_basin_numbers.define(output);
  OceanModel::define_model_state_impl(output);
}

// Derived from pico.cc|hh
void ISMIP6::write_model_state_impl(const File &output) const {

  //m_basin_numbers.write(output);
  OceanModel::write_model_state_impl(output);
}


MaxTimestep ISMIP6::max_timestep_impl(double t) const {
  (void)t;
  // Use this to disable the time step restriction
  return MaxTimestep("ocean ISMIP6");
}

// According to example_melt_NON_LOCAL_param.f90 this formula results in
// a basal melt rate of pure water in m per year!!!
// ... 'Total Antarctic ice-shelf basal mass balance = ', sum(total), ' Gt/yr   [ Ctrl value = 1533.82 ]'
void ISMIP6::update_impl(const Geometry &geometry, double t, double dt) {


  m_log->message(5, "ISMIP6::update_impl(...): t = %.0f, dt = %.0f seconds\n", t, dt);

  const IceModelVec2S &ice_thickness    = geometry.ice_thickness;
  const IceModelVec2S &ice_surface      = geometry.ice_surface_elevation; // debug only
  const IceModelVec2CellType &cell_type = geometry.cell_type;

  // Geometric part of ISMIP6. Updates ice_draft and ice_shelf_mask
  m_physics->update_geometry(ice_thickness, ice_surface, cell_type);


  // Override ice shelf mask from bedmap2 input for testing
  // This option is only available for debugging and has therefore no proper equivalent in the config!
  //#if (PISM_DEBUG == 1)
  {
    std::string option = "-ocean_geometry_override_file";
    options::String filename(option, "Specifies a file containing bedmap2 'draft' and 'ice_shelf_mask'");
    if (filename.is_set()) {
      m_physics->update_geometry_from_file(filename);
    }
  }
  //#endif

  if (m_time_dependent) {
    // get the year of this time step
    const Time &time = *m_grid->ctx()->time();

    // ISMIP6 requires forcing updates every year
    // pism_config:time_stepping.hit_multiples = 1.0;
    double curr_year = time.convert_time_interval(t, "year");
    double last_year = time.convert_time_interval(m_t_last, "year");

    if (t == m_grid->ctx()->time()->start() || curr_year >= last_year + 1.0) {

      m_log->message(3, "ISMIP6::update_impl(): needs an update for this year %f with last year  %f \n", curr_year,
                     last_year);

      // TODO: For now we read in forcing data year by year.
      //       This should be changed later to behave more natural within PISM,
      //       thus use time bounds if present (interpreting data as piecewise-constant)
      //       and use linear interpolation otherwise.
      unsigned int tidx = get_time_index(t, m_tser->times());
      ForcingOptions opt(*m_grid->ctx(), "ocean.ismip6"); // opt required for forcing file name
      m_thermal_forcing.update(opt.filename, tidx);

      // done update
      m_t_last = t; // Note, bedef uses t_final =  t + dt here.
    }               // end update forcing data from netcdf
  }                 // end time dependent

  // lets continue with the current thermal forcing ...

  // for all i,j get ice draft and vertical interpolate TF to the ice draft position
  m_log->message(5, "ISMIP6::update_impl(...): m_forcing start\n");
  m_physics->ocean_forcing_at_draft(m_thermal_forcing, m_forcing);
  m_log->message(5, "ISMIP6::update_impl(...): m_forcing end\n");

  // compute basin averages for non local methods or diagnostics
  m_log->message(5, "ISMIP6::update_impl(...): basin averages start\n");
  m_physics->compute_ocean_forcing_per_basin(m_forcing, m_basin_area, m_basin_average_forcing);
  m_log->message(5, "ISMIP6::update_impl(...): basin averages end\n");

  // basal melt
  m_log->message(5, "ISMIP6::update_impl(...): ice shelf basal melt start\n");
  {
    std::string method = m_config->get_string("ocean.ismip6.thermal_forcing_method");
    m_log->message(5, "ISMIP6::update_impl(...): with method <%s>\n", method.c_str());

    if (method == "quadratic_local") {

      m_physics->local_quadratic_melt(m_forcing, m_basal_melt_rate);

    } else if (method == "quadratic") {

      m_physics->non_local_quadratic_melt(m_basin_average_forcing, m_forcing, m_basal_melt_rate);

    } else {
      // should not be reached
      throw RuntimeError(PISM_ERROR_LOCATION, "ISMIP6 invalid method");
    }
  }
  m_log->message(5, "ISMIP6::update_impl(...): ice shelf basal melt end\n");


  // print basin summary
  if (m_first_update) {
    m_log->message(2, "ISMIP6::update_impl(...): Report basins if first update\n");
    m_physics->report_basin_summary(m_basin_average_forcing, m_basin_area, m_basal_melt_rate);
  }

  // RESULTS
  // shelf_base_mass_flux, melting_point_temperature and melange_back_pressure_fraction
  //
  // convert m fresh water m/s (with gamma0 in m/s) into kg.m-2.s-1
  m_shelf_base_mass_flux->copy_from(m_basal_melt_rate); // m.s-1
  m_shelf_base_mass_flux->scale(1000.0);          // kg.m-2.s-1, rho_fw = 1000 kg.m-3
  melting_point_temperature(ice_thickness, *m_shelf_base_temperature);
  m_melange_back_pressure_fraction->set(0.0);


#if (PISM_DEBUG == 1)
  // dump fields for debug
  if (m_first_update) {
    std::string output_file = m_config->get_string("output.file_name");
    std::string base_name   = output_file.substr(0, output_file.find_last_of('.'));
    std::string dump_file;
    dump_file = base_name + "_-_dump_ice_draft.nc";
    m_physics->ice_draft().dump(dump_file.c_str());
    dump_file = base_name + "_-_dump_ice_shelf_mask.nc";
    m_physics->ice_shelf_mask().dump(dump_file.c_str());
    dump_file = base_name + "_-_dump_thermal_forcing.nc";
    m_forcing.dump(dump_file.c_str());
    dump_file = base_name + "_-_dump_ice_shelf_melt.nc";
    m_basal_melt_rate.dump(dump_file.c_str());
    dump_file = base_name + "_-_dump_shelf_base_mass_flux.nc";
    m_shelf_base_mass_flux->dump(dump_file.c_str());
  }
#endif

  m_first_update = false;
}

// adapted from double Timeseries::operator()(double t) const ...
// piecewise-linear case
unsigned int ISMIP6::get_time_index(double t, std::vector<double> times) const {

  int result = 0;
  auto end   = times.end();
  auto j     = lower_bound(times.begin(), end, t); // binary search
  if (j == end) {
    result = times.size() - 1; // out of range (on the right)
  } else {
    result = (int)(j - times.begin()) - 1;
    if (result < 0) {
      result = 0; // out of range (on the left)
    }
  }

  {
    const Time &time = *m_grid->ctx()->time(); // needed for time conversion
    m_log->message(3, "ISMIP6::get_time_index: idx = %d for time = %0.f s -> year %f \n", result, times[result],
                   time.convert_time_interval(times[result], "year"));
  }


  return result;
}

// Write diagnostic variables to extra files if requested
DiagnosticList ISMIP6::diagnostics_impl() const {

  // TODO:
  //  add diagnostic "total_basal_mass_flux_floating"
  //  similar to basal_mass_flux_floating (average basal mass flux over the reporting interval)
  //  but for each basin and in units Gt/a as in surface_melt_rate
  //  This would be the same information as given std::vector<double> total(m_basin_count, 0.0);

  // todo: drop prefix "ismip6_"
  DiagnosticList result = {
    { "ismip6_basin_number", Diagnostic::wrap(m_physics->basin_numbers()) },
    { "ismip6_ice_shelf_mask", Diagnostic::wrap(m_physics->ice_shelf_mask()) },
    { "ismip6_thermal_forcing_at_depth", Diagnostic::wrap(m_forcing) },
    { "ismip6_basal_melt_rate", Diagnostic::wrap(m_basal_melt_rate) },
    { "ismip6_basal_temperature", Diagnostic::wrap(*m_shelf_base_temperature) },
  };

  return combine(result, OceanModel::diagnostics_impl());
}

//! Copy from PIK::melting_point_temperature(...)
void ISMIP6::melting_point_temperature(const IceModelVec2S &ice_thickness, IceModelVec2S &result) const {
  const double T0          = m_config->get_number("constants.fresh_water.melting_point_temperature"), // K
      beta_CC              = m_config->get_number("constants.ice.beta_Clausius_Clapeyron"),
               g           = m_config->get_number("constants.standard_gravity"),
               ice_density = m_config->get_number("constants.ice.density");


  IceModelVec::AccessList list{ &ice_thickness, &result };

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    const double pressure = ice_density * g * ice_thickness(i, j); // FIXME task #7297
    // result is set to melting point at depth
    result(i, j) = T0 - beta_CC * pressure;
  }
}


} // end of namespace ocean
} // end of namespace pism
