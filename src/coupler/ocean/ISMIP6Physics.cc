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


#include <algorithm> // max_element

#include <algorithm> // std::max
#include <cctype>
#include <cmath> // std::abs
#include <iostream>
#include <numeric> // for std::accumulate
#include <string>
#include <vector>


#include "ISMIP6Physics.hh"
//#include "pism/calving/connected_components.hh"
#include "pism/util/IceModelVec2CellType.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {
namespace ocean {

ISMIP6Physics::ISMIP6Physics(IceGrid::ConstPtr grid) : Component(grid) {

  m_ice_shelf_mask.create(m_grid, "ice_shelf_mask", WITHOUT_GHOSTS);
  m_ice_shelf_mask.set_attrs("diagnostic", "ice shelf mask (floating ice connected to the open ocean)", "", "");
  m_ice_shelf_mask.metadata().set_numbers("flag_values", { 0.0, 1.0 });
  m_ice_shelf_mask.metadata().set_string("flag_meanings", "non_floating ice, floating_ice");

  m_ice_draft.create(m_grid, "ice_draft", WITHOUT_GHOSTS);
  m_ice_draft.set_attrs("diagnostic",
                        "ice_draft",
                        "m", "", "", 0);

  m_basin_numbers.create(m_grid, "basinNumber", WITHOUT_GHOSTS); // check why we have GHOSTS
  m_basin_numbers.set_attrs("climate_forcing",
                            "mask determines basins for ISMIP6 ocean forcing",
                            "", "", "", 0);

  m_deltaT_basin.create(m_grid, "deltaT_basin", WITH_GHOSTS); // check why we have GHOSTS
  m_deltaT_basin.set_attrs("climate_forcing",
                           "basin delta T",
                           "K", "", "", 0);

  m_gamma0      = 9999.0; // should crash
  m_basin_count = 0;
}

ISMIP6Physics::~ISMIP6Physics() {
  // empty
}

const IceModelVec2Int &ISMIP6Physics::basin_numbers() const {
  return m_basin_numbers;
}

const IceModelVec2Int &ISMIP6Physics::ice_shelf_mask() const {
  return m_ice_shelf_mask;
}

const IceModelVec2S &ISMIP6Physics::ice_draft() const {
  return m_ice_draft;
}

/*!
 * Compute masks needed by the ISMIP6 Ocean.
 *
 * After this call ice_shelf_mask() will be up to date.
 */
void ISMIP6Physics::update_geometry(const IceModelVec2S &ice_thickness, const IceModelVec2S &ice_surface,
                                    const IceModelVec2CellType &cell_type) {

  IceModelVec::AccessList list{ &cell_type, &m_ice_shelf_mask, &m_ice_draft, &ice_thickness, &ice_surface };

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    if (cell_type.floating_ice(i, j)) {
      m_ice_shelf_mask(i, j) = 1.0;
    } else {
      m_ice_shelf_mask(i, j) = 0.0;
    }

    // according to the example_melt_NON_LOCAL_param.f90,
    // the draft is computed everywhere from surface and thickness
    m_ice_draft(i, j) = ice_surface(i, j) - ice_thickness(i, j);
    // PISM usually uses
    // double ice_shelf_draft = -(ice_density / sea_water_density) * ice_thickness(i, j);
  }
}

void ISMIP6Physics::update_geometry_from_file(const std::string &filename) {

  m_log->message(2, "  ### DEBUG ################################################################\n");
  m_log->message(2, "  - Reading 'ice_draft' and 'ice_shelf_mask' from '%s'...\n", filename.c_str());
  m_log->message(2, "  ##########################################################################\n");
  m_ice_shelf_mask.regrid(filename, CRITICAL);
  m_ice_draft.regrid(filename, CRITICAL);
}

void ISMIP6Physics::ocean_forcing_at_draft(const IceModelVec3Custom &ocean_forcing3, IceModelVec2S &result) const {

  IceModelVec::AccessList list{ &m_ice_draft, &m_ice_shelf_mask, &ocean_forcing3, &result };

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    if (m_ice_shelf_mask(i, j) > 0.0) {
      //double ice_shelf_draft = -(ice_density / sea_water_density) * ice_thickness(i, j);
      const double ice_shelf_draft = m_ice_draft(i, j);
      // Thermal forcing at ice shelf base
      result(i, j) = ocean_forcing3.getValZocean(i, j, ice_shelf_draft); // getValZocean is very similar to getValZ
    } else {
      result(i, j) = 0.0; // TODO: set to missing value
    }
  }
}

void ISMIP6Physics::init_impl() {

  //
  //  read ocean basins
  //
  auto basin_N_file = m_config->get_string("ocean.ismip6.basin_file");
  if (not basin_N_file.empty()) {
    m_log->message(2, "  - Reading basins from '%s'...\n", basin_N_file.c_str());
    // Note, for IceModelVec2Int the m_interpolation_type = NEAREST;
    m_basin_numbers.regrid(basin_N_file, CRITICAL);
  } else {
    m_log->message(2, "  - Option ocean.ismip6.basin_file is not set. Using default '%d'...\n", 0);
    m_basin_numbers.set(0); // Note, must be zero, as basin count starts at zero
  }

  // It is assumed that the numbers are from 0...lastBasin without gaps!
  m_basin_count = m_basin_numbers.max() + 1;
  m_log->message(2, "  - We have '%d' basin(s)...\n", m_basin_count);
  m_log->message(2, "    basin min=%f, max=%f\n", m_basin_numbers.min(), m_basin_numbers.max());

  //
  //  Read deltaT field and gamma0
  //

  auto deltaT_file = m_config->get_string("ocean.ismip6.deltaT_file");

  if (deltaT_file.empty()) {

    m_deltaT_basin.set(0.0); // todo: read from config
    m_gamma0 = m_config->get_double("ocean.ismip6.gamma0", "meter seconds-1");
    m_log->message(2, "  - WARNING: Option ocean.ismip6.deltaT_file is not set.\n");
    m_log->message(2, "             Using default deltaT '%f'...\n", 0.0);
    m_log->message(2, "             Using default gamma0 '%f'...\n", m_gamma0);

  } else {

    m_log->message(2, "  - Reading deltaT field and gamma0 from '%s'...\n", deltaT_file.c_str());
    PIO file(m_grid->com, "guess_mode", deltaT_file, PISM_READONLY);
    m_deltaT_basin.regrid(file, CRITICAL);
    if (file.inq_var("gamma0")) {
      std::string units = file.get_att_text("gamma0", "units");
      double gamma0     = 0.0;
      file.get_var_double("gamma0", &gamma0);
      m_gamma0 = units::convert(m_sys, gamma0, units, "meter seconds-1"); // convert m/a or m/s -> m/s
    } else {
      m_log->message(2, "  - WARNING: gamma0 not found in '%s'...\n", deltaT_file.c_str());
    }
    file.close();
  }

  // Note, (https://github.com/Unidata/UDUNITS-2/blob/master/lib/udunits2-common.xml)
  //      Julian_year = 365.25 days
  //   Gregorian_year = 365.2425 days
  //    tropical_year = 365.242198781 days
  m_log->message(2, "  - We use gamma0 = %.16e m/s (%16.10f m/a, udunits2) ...\n", m_gamma0,
                 units::convert(m_sys, m_gamma0, "meter seconds-1", "meter years-1"));
}

// could be used for any quantity e.g. temperature, salinity or thermal_forcing,
// thus, do not use m_basin_area or m_basin_average from ISMIP6Physics
void ISMIP6Physics::compute_ocean_forcing_per_basin(const IceModelVec2S &ocean_forcing2,
                                                    std::vector<double> &basin_area,
                                                    std::vector<double> &basin_average) {
#if (PISM_DEBUG == 1)
  m_log->message(2, "ISMIP6Physics: size basin_area=%d basin_average=%d, should be %d\n",
                 static_cast<int>(basin_area.size()), static_cast<int>(basin_average.size()), m_basin_count);
#endif

  // reset
  std::fill(basin_area.begin(), basin_area.end(), 0.0);
  std::fill(basin_average.begin(), basin_average.end(), 0.0);

  // Note, cell_area in stable1.1 is always dx*dy and is not available via -extra_vars.
  // Therefore cell_area is outside the loop in pism1.1.
  // This is different to prev. pism versions!
  auto cell_area = m_grid->cell_area();
  std::vector<int> count(m_basin_count, 0);
  IceModelVec::AccessList list{ &m_basin_numbers, &m_ice_shelf_mask, &ocean_forcing2 };
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    if (m_ice_shelf_mask.as_int(i, j) > 0) {
      const int id = m_basin_numbers.as_int(i, j);
      count[id] += 1; // may not required
      basin_area[id] += cell_area;
      basin_average[id] += (ocean_forcing2(i, j) * cell_area);
    }
  }

  // compute global sums in serial loop, required after parallel loop
  for (int id = 0; id < m_basin_count; id++) {
    count[id]         = GlobalSum(m_grid->com, count[id]);
    basin_area[id]    = GlobalSum(m_grid->com, basin_area[id]);
    basin_average[id] = GlobalSum(m_grid->com, basin_average[id]);

    if (count[id] == 0 || basin_area[id] == 0.0) {
      m_log->message(2, "ISMIP6 ocean WARNING: basin '%d' contains no cells with ice shelf.\n"
                        "This should not cause any problems as the basin average is only used\n"
                        "for melting at ice shelf location (absent in this case).\n",
                     id);
      basin_average[id] = 0.0;
    } else {
      basin_average[id] /= basin_area[id];
    }

#if (PISM_DEBUG == 1)
    m_log->message(2, "ISMIP6Physics: basin = %02d, mean(TF) = %5.3f degC, area = %9.2f km2, counts = %d\n", id,
                   basin_average[id], basin_area[id] * 1.0e-6, count[id]);
#endif
  }
}

int ISMIP6Physics::get_basin_count() const {
  return m_basin_count;
}


void ISMIP6Physics::local_quadratic_melt(const IceModelVec2S &forcing, IceModelVec2S &result) const {

  IceModelVec::AccessList list{ &m_basin_numbers, &m_ice_shelf_mask, &forcing, &m_deltaT_basin, &result };

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    if (m_ice_shelf_mask(i, j) > 0.0) {
      const double deltaT = m_deltaT_basin(i, j);
      const double TF     = forcing(i, j);
      const double tmp    = std::max(TF + deltaT, 0.0);
      // m/s fresh water with gamma0 in m/s (updated Dec 2018 at ISMIP6 pre-AGU meeting)
      result(i, j) = m_gamma0 * m_coeff2 * tmp * tmp;
    } else {
      result(i, j) = 0.0; // TODO: set to missing value
    }
  }
}


void ISMIP6Physics::non_local_quadratic_melt(const std::vector<double> &basin_average_forcing,
                                             const IceModelVec2S &forcing, IceModelVec2S &result) const {


  IceModelVec::AccessList list{ &m_basin_numbers, &m_ice_shelf_mask, &forcing, &m_deltaT_basin, &result };

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    if (m_ice_shelf_mask(i, j) > 0.0) {
      const int id           = m_basin_numbers.as_int(i, j);
      const double deltaT    = m_deltaT_basin(i, j);
      const double averageTF = basin_average_forcing[id];
      const double TF        = forcing(i, j);
      // m/s fresh water with gamma0 in m/s (updated Dec 2018 at ISMIP6 pre-AGU meeting)
      // todo: test modifications proposed in Lipscomb et al., 2020 (doi: 10.5194/tc-2019-334)
      result(i, j) = m_gamma0 * m_coeff2 * (TF + deltaT) * std::abs(averageTF + deltaT);
    } else {
      result(i, j) = 0.0; // TODO: set to missing value
    }
  }

}

void ISMIP6Physics::report_basin_summary(const std::vector<double> &basin_average_forcing,
                                         const std::vector<double> &basin_area,
                                         const IceModelVec2S &basal_melt_rate) const {

  // -------------------------------------------------------------------
  // Testing total melt rates per basin
  // (so that you can check while you adapt this program to your ice sheet model):
  std::vector<double> total(m_basin_count, 0.0);
  auto cell_area = m_grid->cell_area();
  const double sec_per_year = units::convert(m_sys, 1.0, "year", "second"); // use udunits
  m_log->message(2, "ISMIP6: seconds per year used to convert %13.11e\n", sec_per_year);

  IceModelVec::AccessList list{ &m_ice_shelf_mask, &basal_melt_rate, &m_basin_numbers };
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    if (m_ice_shelf_mask(i, j) > 0.0) {
      const int id = m_basin_numbers.as_int(i, j);
      total[id] += (basal_melt_rate(i, j) * cell_area) * 1.0e-9 * sec_per_year;
    }
  }

  // compute global sum, required after parallel loop
  for (int id = 0; id < m_basin_count; id++) {
    total[id] = GlobalSum(m_grid->com, total[id]);
    m_log->message(2, "ISMIP6: melt in basin %02d = %8.2f Gt/a, mean(TF) = %5.3f degC, area = %9.2f km2\n", id,
                   total[id], basin_average_forcing[id], basin_area[id] * 1.0e-6);
  }

  // totals
  m_log->message(2, "ISMIP6: total ice shelf basal melt = %8.2f Gt/a, total ice shelf area = %9.2f km2\n",
                 std::accumulate(total.begin(), total.end(), 0.0),
                 std::accumulate(basin_area.begin(), basin_area.end(), 0.0) * 1.0e-6);
}

} // end of namespace ocean
} // end of namespace pism
