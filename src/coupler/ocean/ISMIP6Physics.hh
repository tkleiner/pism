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


#ifndef ISMIP6PHYSICS_H
#define ISMIP6PHYSICS_H

#include "pism/util/Component.hh"
#include "pism/util/iceModelVec.hh"
#include "pism/util/iceModelVec3Custom.hh"

namespace pism {

class IceModelVec2S;
class IceModelVec2CellType;

namespace ocean {


/*!
 * This class isolates geometric and other computations performed by the ISMIP6 ocean model.
 */
class ISMIP6Physics : public Component {
public:
  ISMIP6Physics(IceGrid::ConstPtr grid);
  virtual ~ISMIP6Physics();

  void init_impl();

  void update_geometry(const IceModelVec2S &ice_thickness, const IceModelVec2S &ice_surface,
                       const IceModelVec2CellType &cell_type);


  void update_geometry_from_file(const std::string &filename);

  // get thermal forcing at ice shelf draft elevation
  void ocean_forcing_at_draft(const IceModelVec3Custom &ocean_forcing3, IceModelVec2S &result) const;

  void compute_ocean_forcing_per_basin(const IceModelVec2S &ocean_frocing2, std::vector<double> &basin_area,
                                       std::vector<double> &basin_average);

  // This should be part of some helper class
  void local_quadratic_melt(const IceModelVec2S &forcing, IceModelVec2S &result) const;
  void non_local_quadratic_melt(const std::vector<double> &basin_average_forcing, const IceModelVec2S &forcing,
                                IceModelVec2S &result) const;

  void report_basin_summary(const std::vector<double> &basin_average_forcing, const std::vector<double> &basin_area,
                            const IceModelVec2S &basal_melt_rate) const;


  const IceModelVec2Int &ice_shelf_mask() const;
  const IceModelVec2Int &basin_numbers() const;
  const IceModelVec2S &ice_draft() const;
  int get_basin_count() const;


private:
  // no interpolation is allowed; implemented and stored as floating point
  IceModelVec2Int m_basin_numbers;
  IceModelVec2S m_deltaT_basin; //! Temperature offset per basin
  double m_gamma0;              //! uniform melt coefficient to be read from input file


  // Those fields are only used to allow for bypassing the ice draft and mask from an external file.
  // Consider to remove this in later versions and use pism mask with a local ice draft
  IceModelVec2Int m_ice_shelf_mask;
  IceModelVec2S m_ice_draft;

  // number of ocean basins
  int m_basin_count;

  // fixme: parameter differ between ISMIP6 and PISM defaults
  const double rho_i = 918.0,  // Ice density (kg/m^3)
      rho_sw         = 1028.0, // Sea water density (kg/m^3)
      L_f            = 3.34e5, // Latent heat of fusion (J/kg)
      c_pw           = 3974.0; // Specific heat of sea water (J/kg/K)

  // coefficient to non-dimensionalize TF
  const double coeff1   = rho_sw * c_pw / (rho_i * L_f); // units (K^-1) linear
  const double m_coeff2 = coeff1 * coeff1;               // units (K^-2) quadratic
};

} // end of namespace ocean
} // end of namespace pism

#endif /* ISMIP6PHYSICS_H */
