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


#ifndef _ISMIP6_H_
#define _ISMIP6_H_

#include "CompleteOceanModel.hh"
#include "pism/util/IceModelVec2CellType.hh"
#include "pism/util/iceModelVec2T.hh"
#include "pism/util/iceModelVec3Custom.hh"


namespace pism {
namespace ocean {

class ISMIP6Physics;

class ISMIP6 : public CompleteOceanModel {

public:
  ISMIP6(IceGrid::ConstPtr g);
  virtual ~ISMIP6();

protected:
  void init_impl(const Geometry &geometry);
  void update_impl(const Geometry &geometry, double t, double dt);
  MaxTimestep max_timestep_impl(double t) const;

  void define_model_state_impl(const File &output) const;
  void write_model_state_impl(const File &output) const;

  std::map<std::string, Diagnostic::Ptr> diagnostics_impl() const;

  // copy from PIK::melting_point_temperature()
  void melting_point_temperature(const IceModelVec2S &ice_thickness, IceModelVec2S &result) const;

  unsigned int m_Moz; //!< number of vertical levels within the ocean
  double m_Loz;       //!< thickness of the ocean layer, in meters

  // split 4d TF(x,y,z,t) into 3d TF(x,y,z) and time
  IceModelVec3Custom::Ptr m_thermal_forcing; //! spatial data in TF(x,y,z) for given time
  Timeseries *m_tser;                        //! time dimension data from TF
  //todo: use also tidxPrev, tidxCurr for interpolation in time also
  unsigned int tidxNext;

protected:
  IceModelVec2S m_forcing;         //! Thermal forcing at ice shelf base
  IceModelVec2S m_basal_melt_rate; // in units m/s or m/a scaled from

private:
  std::unique_ptr<ISMIP6Physics> m_physics;

  // store for sanity check during simulation
  double m_z_ocean_surface_layer, m_z_ocean_bottom_layer;

  bool m_first_update;

  int m_basin_count;
  std::vector<double> m_basin_area;            // only known after update
  std::vector<double> m_basin_average_forcing; // only known after update

  //! flag indicating if we have TF(t,z,y,x) or TF(z,y,x)
  bool m_time_dependent;

  //! time of the last forcing update (similar to BedDef.hh m_t_beddef_last)
  double m_t_last;

  unsigned int get_time_index(double t, std::vector<double> times) const;
};

} // end of namespace ocean
} // end of namespace pism
#endif /* _ISMIP6_H_ */
