/**
 * tlibs2
 * magnon dynamics
 * @author Tobias Weber <tweber@ill.fr>
 * @date january-2022
 * @license GPLv3, see 'LICENSE' file
 *
 * References:
 *   - (Toth 2015) S. Toth and B. Lake, J. Phys.: Condens. Matter 27 166002 (2015):
 *                 https://doi.org/10.1088/0953-8984/27/16/166002
 *                 https://arxiv.org/abs/1402.6069
 *   - (Heinsdorf 2021) N. Heinsdorf, manual example calculation for a simple
 *                      ferromagnetic case, personal communications, 2021/2022.
 *
 * @desc This file implements the formalism given by (Toth 2015).
 * @desc For further references, see the 'LITERATURE' file.
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2023  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ----------------------------------------------------------------------------
 */

#ifndef __TLIBS2_MAGDYN_H__
#define __TLIBS2_MAGDYN_H__

#include <vector>
#include <array>
#include <tuple>
#include <unordered_map>
#include <string>

#include <algorithm>
#include <numeric>

#include <iostream>
#include <fstream>
#include <iomanip>

#include <boost/container_hash/hash.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#ifndef USE_LAPACK
	#define USE_LAPACK 1
#endif
#include "maths.h"
#include "units.h"
#include "phys.h"
#include "algos.h"
#include "expr.h"

// enables debug output
//#define __TLIBS2_MAGDYN_DEBUG_OUTPUT__
//#define __TLIBS2_MAGDYN_DEBUG_PY_OUTPUT__


namespace tl2_mag {

// ----------------------------------------------------------------------------
// helper functions
// ----------------------------------------------------------------------------

/**
 * rotate spin vector for incommensurate structures, i.e. helices
 */
template<class t_mat, class t_vec, class t_real = typename t_mat::value_type>
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_mat<t_mat> && tl2::is_vec<t_vec>
#endif
void rotate_spin_incommensurate(t_vec& spin_vec,
	const t_vec& sc_vec, const t_vec& ordering, const t_vec& rotaxis,
	t_real eps = std::numeric_limits<t_real>::epsilon())
{
	t_real sc_angle = t_real(2) * tl2::pi<t_real> *
		tl2::inner<t_vec>(ordering, sc_vec);

	if(!tl2::equals_0<t_real>(sc_angle, t_real(eps)))
	{
		t_mat sc_rot = tl2::rotation<t_mat, t_vec>(rotaxis, sc_angle);
		spin_vec = sc_rot * spin_vec;
	}
}


/**
 * polarisation matrix
 * @see https://doi.org/10.1088/1361-6463/aa7573
 */
template<class t_mat, class t_cplx = typename t_mat::value_type>
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_mat<t_mat>
#endif
t_mat get_polarisation(int channel = 0, bool in_chiral_basis = true)
{
	if(in_chiral_basis)
	{
		t_mat pol = tl2::zero<t_mat>(3);

		// just pick the selected component on the diagonal
		if(channel >=0 && channel < 3)
			pol(channel, channel) = t_cplx(1);

		return pol;
	}
	else
	{
		constexpr const t_cplx halfi = t_cplx(0, 0.5);
		constexpr const t_cplx half = t_cplx(0.5, 0);

		// TODO: check coordinate system
		switch(channel)
		{
			case 0: return tl2::create<t_mat>({
				{   half, +halfi,  0 },
				{ -halfi,   half,  0 },
				{      0,      0,  0 } });
			case 1: return tl2::create<t_mat>({
				{   half, -halfi,  0 },
				{ +halfi,   half,  0 },
				{      0,      0,  0 } });
			case 2: return tl2::create<t_mat>({
				{      0,      0,  0 },
				{      0,      0,  0 },
				{      0,      0,  1 } });
		}
	}

	return tl2::zero<t_mat>(3);
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// input- and output struct templates
// ----------------------------------------------------------------------------
/**
 * magnetic sites
 */
template<class t_vec_real, class t_mat, class t_real, class t_size>
struct t_MagneticSite
{
	std::string name{};          // identifier
	t_size index{};              // index

	t_vec_real pos{};            // magnetic site position

	std::string spin_dir[3];     // expression for spin direction
	std::string spin_ortho[3];   // spin orthogonal plane

	t_real spin_mag{};           // spin magnitude
	t_mat g{};                   // g factor
};


/**
 * temporary per-site calculation results
 */
template<class t_vec>
struct t_MagneticSiteCalc
{
	t_vec spin_dir{};        // spin direction

	t_vec u{};               // spin orthogonal plane vector 1
	t_vec v{};               // spin orthogonal plane vector 2
	t_vec u_conj{};
};


/**
 * couplings between magnetic sites
 */
template<class t_vec_real, class t_size>
struct t_ExchangeTerm
{
	std::string name{};      // identifier
	t_size index{};          // index

	t_size site1{};          // index of first magnetic site
	t_size site2{};          // index of second magnetic site
	t_vec_real dist{};       // distance between unit cells

	std::string J{};         // parsable expression for Heisenberg interaction
	std::string dmi[3];      // parsable expression for Dzyaloshinskij-Moriya interaction
	std::string Jgen[3][3];  // parsable expression for a general exchange interaction
};


/**
 * temporary per-term calculation results
 */
template<class t_mat, class t_vec, class t_cplx>
struct t_ExchangeTermCalc
{
	t_cplx J{};              // Heisenberg interaction
	t_vec dmi{};             // Dzyaloshinskij-Moriya interaction
	t_mat Jgen{};            // general exchange interaction
};


/**
 * terms related to an external magnetic field
 */
template<class t_vec_real, class t_real>
struct t_ExternalField
{
	bool align_spins{};      // align spins along external field
	t_vec_real dir{};        // field direction
	t_real mag{};            // field magnitude
};


/**
 * eigenenergies and spin-spin correlation matrix
 */
template<class t_mat, class t_real>
struct t_EnergyAndWeight
{
	t_real E{};

	// full dynamical structure factor
	t_mat S{};
	t_real weight_full{};
	t_real weight_channel_full[3] = {0., 0., 0.};

	// projected dynamical structure factor for neutron scattering
	t_mat S_perp{};
	t_real weight{};
	t_real weight_channel[3] = {0., 0., 0.};
};


/**
 * variables for the expression parser
 */
template<class t_cplx>
struct t_Variable
{
	std::string name{};
	t_cplx value{};
};
// ----------------------------------------------------------------------------



/**
 * calculates magnon dynamics,
 * implementing the formalism given in (Toth 2015)
 */
template<
	class t_mat, class t_vec,
	class t_mat_real, class t_vec_real,
	class t_cplx = typename t_mat::value_type,
	class t_real = typename t_mat_real::value_type,
	class t_size = std::size_t>
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_mat<t_mat> && tl2::is_vec<t_vec> &&
	tl2::is_mat<t_mat_real> && tl2::is_vec<t_vec_real>
#endif
class MagDyn
{
public:
	// --------------------------------------------------------------------
	// structs and types
	// --------------------------------------------------------------------
	using MagneticSite = t_MagneticSite<t_vec_real, t_mat, t_real, t_size>;
	using MagneticSiteCalc = t_MagneticSiteCalc<t_vec>;
	using ExchangeTerm = t_ExchangeTerm<t_vec_real, t_size>;
	using ExchangeTermCalc = t_ExchangeTermCalc<t_mat, t_vec, t_cplx>;
	using ExternalField = t_ExternalField<t_vec_real, t_real>;
	using EnergyAndWeight = t_EnergyAndWeight<t_mat, t_real>;
	using Variable = t_Variable<t_cplx>;

	using t_indices = std::pair<t_size, t_size>;
	using t_Jmap = std::unordered_map<t_indices, t_mat, boost::hash<t_indices>>;
	// --------------------------------------------------------------------


public:
	MagDyn() = default;
	~MagDyn() = default;

	MagDyn(const MagDyn&) = default;
	MagDyn& operator=(const MagDyn&) = default;


	// --------------------------------------------------------------------
	// cleanup functions
	// --------------------------------------------------------------------
	/**
	 * clear all
	 */
	void Clear()
	{
		ClearVariables();
		ClearMagneticSites();
		ClearExchangeTerms();
		ClearExternalField();

		// clear temperature, -1: don't use
		m_temperature = -1.;

		// clear ordering wave vector
		m_ordering = tl2::zero<t_vec_real>(3);

		// reset rotation axis
		m_rotaxis = tl2::create<t_vec_real>({1., 0., 0.});
	}


	/**
	 * clear all parser variables
	 */
	void ClearVariables()
	{
		m_variables.clear();
	}


	/**
	 * clear all magnetic sites
	 */
	void ClearMagneticSites()
	{
		m_sites.clear();
		m_sites_calc.clear();
	}


	/**
	 * clear all couplings
	 */
	void ClearExchangeTerms()
	{
		m_exchange_terms.clear();
		m_exchange_terms_calc.clear();
	}


	/**
	 * clear the external field settings
	 */
	void ClearExternalField()
	{
		m_field.dir.clear();
		m_field.mag = 0.;
		m_field.align_spins = false;
	}
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	// getter
	// --------------------------------------------------------------------
	const std::vector<Variable>& GetVariables() const { return m_variables; }
	const std::vector<MagneticSite>& GetMagneticSites() const { return m_sites; }
	const std::vector<MagneticSiteCalc>& GetMagneticSitesCalc() const { return m_sites_calc; }
	const std::vector<ExchangeTerm>& GetExchangeTerms() const { return m_exchange_terms; }
	const std::vector<ExchangeTermCalc>& GetExchangeTermsCalc() const { return m_exchange_terms_calc; }

	const ExternalField& GetExternalField() const { return m_field; }
	const t_vec_real& GetRotationAxis() const { return m_rotaxis; }
	const t_vec_real& GetOrderingWavevector() const { return m_ordering; }

	t_real GetTemperature() const { return m_temperature; }
	t_real GetBoseCutoffEnergy() const { return m_bose_cutoff; }


	bool IsIncommensurate() const
	{
		return m_is_incommensurate || m_force_incommensurate;
	}


	/**
	 * get number of magnetic sites with the given name (to check if the name is unique)
	 */
	std::vector<const MagneticSite*> FindMagneticSites(const std::string& name) const
	{
		std::vector<const MagneticSite*> sites;

		for(const auto& site : m_sites)
		{
			if(site.name == name)
				sites.push_back(&site);
		}
		return sites;
	}
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	// setter
	// --------------------------------------------------------------------
	void SetEpsilon(t_real eps) { m_eps = eps; }
	void SetPrecision(int prec) { m_prec = prec; }

	void SetTemperature(t_real T) { m_temperature = T; }
	void SetBoseCutoffEnergy(t_real E) { m_bose_cutoff = E; }
	void SetUniteDegenerateEnergies(bool b) { m_unite_degenerate_energies = b; }
	void SetForceIncommensurate(bool b) { m_force_incommensurate = b; }

	void SetPhaseSign(t_real sign) { m_phase_sign = sign; }
	void SetCholeskyMaxTries(t_size max_tries) { m_tries_chol = max_tries; }
	void SetCholeskyInc(t_real delta) { m_delta_chol = delta; }


	void SetExternalField(const ExternalField& field)
	{
		m_field = field;
		t_real len = tl2::norm<t_vec_real>(m_field.dir);
		if(!tl2::equals_0<t_real>(len, m_eps))
			m_field.dir /= len;
	}


	void RotateExternalField(const t_vec_real& axis, t_real angle)
	{
		t_mat_real rot = tl2::rotation<t_mat_real, t_vec_real>(
			axis, angle, false);
		m_field.dir = rot * m_field.dir;
	}


	void RotateExternalField(t_real x, t_real y, t_real z, t_real angle)
	{
		RotateExternalField(tl2::create<t_vec_real>({x, y, z}), angle);
	}


	void SetOrderingWavevector(const t_vec_real& ordering)
	{
		m_ordering = ordering;
		m_is_incommensurate = !tl2::equals_0<t_vec_real>(m_ordering, m_eps);
	}


	void SetCalcHamiltonian(bool H, bool Hp, bool Hm)
	{
		m_calc_H = H;
		m_calc_Hp = Hp;
		m_calc_Hm = Hm;
	}


	void SetRotationAxis(const t_vec_real& axis)
	{
		m_rotaxis = axis;
		t_real len = tl2::norm<t_vec_real>(m_rotaxis);
		if(!tl2::equals_0<t_real>(len, m_eps))
			m_rotaxis /= len;
	}


	void AddVariable(Variable&& var)
	{
		m_variables.emplace_back(std::forward<Variable&&>(var));
	}


	void SetVariable(Variable&& var)
	{
		// is a variable with the same name already registered?
		auto iter = std::find_if(m_variables.begin(), m_variables.end(),
			[&var](const auto& thevar)
		{
			return thevar.name == var.name;
		});

		if(iter == m_variables.end())
		{
			// add a new variable
			AddVariable(std::forward<Variable&&>(var));
		}
		else
		{
			// replace the value of an existing variable
			iter->value = var.value;
		}
	}


	void AddMagneticSite(MagneticSite&& site, bool set_index = true)
	{
		if(set_index)
			site.index = GetMagneticSites().size();
		m_sites.emplace_back(std::forward<MagneticSite&&>(site));
	}


	void AddExchangeTerm(ExchangeTerm&& term, bool set_index = true)
	{
		if(set_index)
			term.index = GetExchangeTerms().size();
		m_exchange_terms.emplace_back(std::forward<ExchangeTerm&&>(term));
	}
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	/**
	 * get an expression parser object with registered variables
	 */
	tl2::ExprParser<t_cplx> GetExprParser() const
	{
		tl2::ExprParser<t_cplx> parser;

		// register all variables
		parser.SetAutoregisterVariables(false);
		for(const Variable& var : m_variables)
			parser.register_var(var.name, var.value);

		return parser;
	}
	// --------------------------------------------------------------------


	/**
	 * get the needed supercell ranges from the exchange terms
	 */
	std::tuple<t_vec_real, t_vec_real> GetSupercellMinMax() const
	{
		t_vec_real min = tl2::zero<t_vec_real>(3);
		t_vec_real max = tl2::zero<t_vec_real>(3);

		for(const ExchangeTerm& term : m_exchange_terms)
		{
			for(t_size i=0; i<3; ++i)
			{
				min[i] = std::min(min[i], term.dist[i]);
				max[i] = std::max(max[i], term.dist[i]);
			}
		}

		return std::make_tuple(min, max);
	}


	// --------------------------------------------------------------------
	// calculation functions
	// --------------------------------------------------------------------
	/**
	 * calculate the spin rotation trafo for the magnetic sites
	 */
	void CalcMagneticSites()
	{
		const t_size num_sites = m_sites.size();
		if(num_sites == 0)
			return;

		m_sites_calc.clear();
		m_sites_calc.reserve(num_sites);

		bool use_field =
			(!tl2::equals_0<t_real>(m_field.mag, m_eps) || m_field.align_spins)
			&& m_field.dir.size() == 3;

		if(use_field)
		{
			// rotate field to [001] direction
			m_rot_field = tl2::convert<t_mat>(
				tl2::trans<t_mat_real>(
					tl2::rotation<t_mat_real, t_vec_real>(
						-m_field.dir, m_zdir, &m_rotaxis, m_eps)));

#ifdef __TLIBS2_MAGDYN_DEBUG_OUTPUT__
			std::cout << "Field rotation from:\n";
			tl2::niceprint(std::cout, -m_field.dir, 1e-4, 4);
			std::cout << "\nto:\n";
			tl2::niceprint(std::cout, m_zdir, 1e-4, 4);
			std::cout << "\nmatrix:\n";
			tl2::niceprint(std::cout, m_rot_field, 1e-4, 4);
			std::cout << std::endl;
#endif
		}

		tl2::ExprParser parser = GetExprParser();

		for(const MagneticSite& site : m_sites)
		{
			try
			{
				MagneticSiteCalc site_calc{};
				bool has_explicit_uv = true;

				site_calc.spin_dir = tl2::zero<t_vec>(3);
				site_calc.v = tl2::zero<t_vec>(3);
				site_calc.u = tl2::zero<t_vec>(3);
				site_calc.u_conj = tl2::zero<t_vec>(3);

				for(t_size dir_idx=0; dir_idx<3; ++dir_idx)
				{
					// non-empty spin direction string?
					if(site.spin_dir[dir_idx].size())
					{
						if(bool dir_ok = parser.parse(site.spin_dir[dir_idx]); dir_ok)
						{
							site_calc.spin_dir[dir_idx] = parser.eval();
							site_calc.v[dir_idx] = site_calc.spin_dir[dir_idx];
						}
						else
						{
							std::cerr << "Error parsing spin direction \""
								<< site.spin_dir[dir_idx] << "\""
								<< " for site " << site.index
								<< " and component " << dir_idx
								<< "." << std::endl;
						}
					}

					// non-empty spin-plane string?
					if(site.spin_ortho[dir_idx].size())
					{
						if(bool dir_ok = parser.parse(site.spin_ortho[dir_idx]); dir_ok)
						{
							site_calc.u[dir_idx] = parser.eval();
							site_calc.u_conj[dir_idx] = std::conj(site_calc.u[dir_idx]);
						}
						else
						{
							has_explicit_uv = false;

							std::cerr << "Error parsing spin orthogonal plane \""
								<< site.spin_ortho[dir_idx] << "\""
								<< " for site " << site.index
								<< " and component " << dir_idx
								<< "." << std::endl;
						}
					}
					else
					{
						has_explicit_uv = false;
					}
				}

				// spin rotation of equation (9) from (Toth 2015)
				if(m_field.align_spins)
				{
					std::tie(site_calc.u, site_calc.v) = R_to_uv(m_rot_field);
				}
				else
				{
					if(!has_explicit_uv)
					{
						// calculate u and v from the spin rotation
						std::tie(site_calc.u, site_calc.v) =
							spin_to_uv(site_calc.spin_dir);
					}

					// TODO: normalise the v vector as well as the real and imaginary u vectors
					// in case they are explicitly given

#ifdef __TLIBS2_MAGDYN_DEBUG_OUTPUT__
					std::cout << "Site " << site.index << " u = "
						<< site_calc.u[0] << " " << site_calc.u[1] << " " << site_calc.u[2]
						<< std::endl;
					std::cout << "Site " << site.index << " v = "
						<< site_calc.v[0] << " " << site_calc.v[1] << " " << site_calc.v[2]
						<< std::endl;
#endif
				}

				site_calc.u_conj = tl2::conj(site_calc.u);
				m_sites_calc.emplace_back(std::move(site_calc));
			}
			catch(const std::exception& ex)
			{
				std::cerr << ex.what() << std::endl;
			}
		}
	}


	/**
	 * parse the exchange term expressions
	 */
	void CalcExchangeTerms()
	{
		if(m_exchange_terms.size() == 0)
			return;

		tl2::ExprParser parser = GetExprParser();

		m_exchange_terms_calc.clear();
		m_exchange_terms_calc.reserve(m_exchange_terms.size());

		for(const ExchangeTerm& term : m_exchange_terms)
		{
			ExchangeTermCalc calc;

			try
			{
				// symmetric interaction
				if(parser.parse(term.J))
				{
					t_cplx J = parser.eval();
					calc.J = J;
				}
				else
				{
					std::cerr << "Error parsing J term \""
						<< term.J << "\"."
						<< std::endl;
				}


				// dmi interaction
				calc.dmi = tl2::zero<t_vec>(3);

				for(t_size dmi_idx=0; dmi_idx<3; ++dmi_idx)
				{
					// empty string?
					if(!term.dmi[dmi_idx].size())
						continue;

					if(parser.parse(term.dmi[dmi_idx]))
					{
						calc.dmi[dmi_idx] = parser.eval();
					}
					else
					{
						std::cerr << "Error parsing DMI term \""
							<< term.dmi[dmi_idx]
							<< "\" (index " << dmi_idx << ")"
							<< "." << std::endl;
					}
				}


				// general exchange interaction
				calc.Jgen = tl2::zero<t_mat>(3, 3);

				for(t_size i=0; i<calc.Jgen.size1(); ++i)
				for(t_size j=0; j<calc.Jgen.size2(); ++j)
				{
					// empty string?
					if(!term.Jgen[i][j].size())
						continue;

					if(parser.parse(term.Jgen[i][j]))
					{
						calc.Jgen(i, j) = parser.eval();
					}
					else
					{
						std::cerr << "Error parsing general term \""
							<< term.Jgen[i][j]
							<< "\" (indices " << i << ", " << j << ")"
							<< "." << std::endl;
					}
				}
			}
			catch(const std::exception& ex)
			{
				std::cerr << ex.what() << std::endl;
			}

			m_exchange_terms_calc.emplace_back(std::move(calc));
		}
	}


	/**
	 * update the site and term indices
	 */
	void CalcIndices()
	{
		for(t_size site_idx=0; site_idx<m_sites.size(); ++site_idx)
		{
			MagneticSite& site = m_sites[site_idx];
			site.index = site_idx;
		}

		for(t_size term_idx=0; term_idx<m_exchange_terms.size(); ++term_idx)
		{
			ExchangeTerm& term = m_exchange_terms[term_idx];
			term.index = term_idx;
		}
	}


	/**
	 * calculate the real-space interaction matrix J of
	 * equations (10) - (13) from (Toth 2015)
	 */
	t_mat CalcRealJ(const ExchangeTerm& term) const
	{
		if(term.index >= m_exchange_terms_calc.size())
		{
			std::cerr << "Error: Coupling terms not yet calculated." << std::endl;
			return t_mat{};
		}

		const ExchangeTermCalc& term_calc = m_exchange_terms_calc[term.index];

		// symmetric part of the exchange interaction matrix, see (Toth 2015) p. 2
		t_mat J = tl2::diag<t_mat>(
			tl2::create<t_vec>({ term_calc.J, term_calc.J, term_calc.J }));

		// dmi as anti-symmetric part of interaction matrix
		// using a cross product matrix, see (Toth 2015) p. 2
		if(term_calc.dmi.size() == 3)
			J += tl2::skewsymmetric<t_mat, t_vec>(-term_calc.dmi);

		// general J matrix
		if(term_calc.Jgen.size1() == 3 && term_calc.Jgen.size2() == 3)
			J += term_calc.Jgen;

		// incommensurate case: rotation wrt magnetic unit cell
		// equations (21), (6), (2) as well as section 10 from (Toth 2015)
		if(IsIncommensurate())
		{
			t_real rot_UC_angle = s_twopi * tl2::inner<t_vec_real>(m_ordering, term.dist);
			if(!tl2::equals_0<t_real>(rot_UC_angle, m_eps))
			{
				t_mat rot_UC = tl2::convert<t_mat>(
					tl2::rotation<t_mat_real, t_vec_real>(
						m_rotaxis, rot_UC_angle));
				J = J * rot_UC;

#ifdef __TLIBS2_MAGDYN_DEBUG_OUTPUT__
				std::cout << "Coupling rot_UC = " << term.index << ":\n";
				tl2::niceprint(std::cout, rot_UC, 1e-4, 4);
#endif
			}
		}

		return J;
	}


	/**
	 * calculate the reciprocal interaction matrices J(Q) and J(-Q) of
	 * equations (12) and (14) from (Toth 2015)
	 */
	std::tuple<t_Jmap, t_Jmap> CalcReciprocalJs(const t_vec_real& Qvec) const
	{
		t_Jmap J_Q, J_Q0;

		// no exchange terms given
		if(m_exchange_terms.size() == 0)
			return std::make_tuple(J_Q, J_Q0);

		if(m_exchange_terms.size() != m_exchange_terms_calc.size())
		{
			std::cerr << "Error: Coupling terms not yet calculated." << std::endl;
			return std::make_tuple(J_Q, J_Q0);
		}

		// iterate couplings to pre-calculate corresponding J matrices
		for(const ExchangeTerm& term : m_exchange_terms)
		{
			// insert or add an exchange matrix at the given site indices
			auto insert_or_add = [](t_Jmap& J, const t_indices& indices, const t_mat& J33)
			{
				if(auto iter = J.find(indices); iter != J.end())
					iter->second += J33;
				else
					J.emplace(std::move(std::make_pair(indices, J33)));
			};

			if(term.site1 >= m_sites.size() || term.site2 >= m_sites.size())
			{
				std::cerr << "Error: Site index out of bounds for coupling term "
					<< term.index << "." << std::endl;
				continue;
			}

			const t_indices indices = std::make_pair(term.site1, term.site2);
			const t_indices indices_t = std::make_pair(term.site2, term.site1);

			t_mat J = CalcRealJ(term);
			if(J.size1() == 0 || J.size2() == 0)
				continue;

			// get J in reciprocal space by fourier trafo
			// equations (14), (12), (11), and (52) from (Toth 2015)
			insert_or_add(J_Q, indices, J *
				std::exp(m_phase_sign * s_imag * s_twopi *
					tl2::inner<t_vec_real>(term.dist, Qvec)));

			t_mat J_T = tl2::trans(J);
			insert_or_add(J_Q, indices_t, J_T *
				std::exp(m_phase_sign * s_imag * s_twopi *
					tl2::inner<t_vec_real>(term.dist, -Qvec)));

			insert_or_add(J_Q0, indices, J);
			insert_or_add(J_Q0, indices_t, J_T);
		}  // end of iteration over couplings

		return std::make_tuple(J_Q, J_Q0);
	}


	/**
	 * get the hamiltonian at the given momentum
	 * (CalcMagneticSites() needs to be called once before this function)
	 * @note implements the formalism given by (Toth 2015)
	 * @note a first version for a simplified ferromagnetic dispersion was based on (Heinsdorf 2021)
	 */
	t_mat CalcHamiltonian(const t_vec_real& Qvec) const
	{
		const t_size num_sites = m_sites.size();

		// no sites given
		if(num_sites == 0)
			return t_mat{};

		if(num_sites != m_sites_calc.size())
		{
			std::cerr << "Error: Sites not yet calculated." << std::endl;
			return t_mat{};
		}

		// build the interaction matrices J(Q) and J(-Q) of
		// equations (12) and (14) from (Toth 2015)
		auto [J_Q, J_Q0] = CalcReciprocalJs(Qvec);

		// create the hamiltonian of equation (25) and (26) from (Toth 2015)
		t_mat A = tl2::create<t_mat>(num_sites, num_sites);
		t_mat A_conj_mQ = tl2::create<t_mat>(num_sites, num_sites);
		t_mat B = tl2::create<t_mat>(num_sites, num_sites);
		t_mat C = tl2::zero<t_mat>(num_sites, num_sites);

		bool use_field = !tl2::equals_0<t_real>(m_field.mag, m_eps)
			&& m_field.dir.size() == 3;

		// iterate magnetic sites
		for(t_size i=0; i<num_sites; ++i)
		{
			// get the pre-calculated u and v vectors for the commensurate case
			const t_vec& u_i = m_sites_calc[i].u;
			const t_vec& u_conj_i = m_sites_calc[i].u_conj;
			const t_vec& v_i = m_sites_calc[i].v;
			t_real S_i = m_sites[i].spin_mag;

			for(t_size j=0; j<num_sites; ++j)
			{
				// get the pre-calculated u and v vectors for the commensurate case
				const t_vec& u_j = m_sites_calc[j].u;
				const t_vec& u_conj_j = m_sites_calc[j].u_conj;
				const t_vec& v_j = m_sites_calc[j].v;
				t_real S_j = m_sites[j].spin_mag;

				// get the pre-calculated exchange matrices for the (i, j) coupling
				const t_indices indices_ij = std::make_pair(i, j);
				const t_mat* J_Q33 = nullptr;
				const t_mat* J_Q033 = nullptr;
				if(auto iter = J_Q.find(indices_ij); iter != J_Q.end())
					J_Q33 = &iter->second;
				if(auto iter = J_Q0.find(indices_ij); iter != J_Q0.end())
					J_Q033 = &iter->second;

				if(J_Q33 && J_Q033)
				{
					t_real SiSj = 0.5 * std::sqrt(S_i*S_j);

					// equation (26) from (Toth 2015)
					A(i, j) = SiSj * tl2::inner_noconj<t_vec>(u_i, (*J_Q33) * u_conj_j);
					A_conj_mQ(i, j) = SiSj * tl2::inner_noconj<t_vec>(u_conj_i, (*J_Q33) * u_j);
					B(i, j) = SiSj * tl2::inner_noconj<t_vec>(u_i, (*J_Q33) * u_j);
					C(i, i) += S_j * tl2::inner_noconj<t_vec>(v_i, (*J_Q033) * v_j);
				}
			}  // end of iteration over j sites

			// include external field, equation (28) from (Toth 2015)
			if(use_field)
			{
				t_vec B = tl2::convert<t_vec>(-m_field.dir) * m_field.mag;

				t_vec gv = m_sites[i].g * v_i;
				t_cplx Bgv = tl2::inner_noconj<t_vec>(B, gv);

				// bohr magneton in [meV/T]
				constexpr const t_real muB = tl2::mu_B<t_real>
					/ tl2::meV<t_real> * tl2::tesla<t_real>;

				A(i, i) -= muB * Bgv;
				A_conj_mQ(i, i) -= std::conj(muB * Bgv);
			}
		}  // end of iteration over i sites

		// equation (25) from (Toth 2015)
		t_mat H = tl2::zero<t_mat>(num_sites*2, num_sites*2);
		tl2::set_submat(H, A - C, 0, 0);
		tl2::set_submat(H, B, 0, num_sites);
		tl2::set_submat(H, tl2::herm(B), num_sites, 0);
		tl2::set_submat(H, A_conj_mQ - C, num_sites, num_sites);

		return H;
	}


	/**
	 * get the energies from a hamiltonian
	 * @note implements the formalism given by (Toth 2015)
	 */
	std::vector<EnergyAndWeight> CalcEnergiesFromHamiltonian(
		t_mat _H, const t_vec_real& Qvec,
		bool only_energies = false) const
	{
		const t_size num_sites = m_sites.size();
		if(num_sites == 0 || _H.size1() == 0)
			return {};

		// equation (30) from (Toth 2015)
		t_mat g_sign = tl2::zero<t_mat>(num_sites*2, num_sites*2);
		for(t_size i=0; i<num_sites; ++i)
			g_sign(i, i) = 1.;
		for(t_size i=num_sites; i<2*num_sites; ++i)
			g_sign(i, i) = -1.;

		// equation (31) from (Toth 2015)
		t_mat C_mat;
		t_size chol_try = 0;
		for(; chol_try<m_tries_chol; ++chol_try)
		{
			auto [chol_ok, _C] = tl2_la::chol<t_mat>(_H);

			if(chol_ok)
			{
				C_mat = std::move(_C);
				break;
			}
			else
			{
				if(chol_try >= m_tries_chol-1)
				{
					using namespace tl2_ops;
					std::cerr << "Warning: Cholesky decomposition failed for Q = "
						<< Qvec << "." << std::endl;
					C_mat = std::move(_C);
					break;
				}

				// try forcing the hamilton to be positive definite
				for(t_size i=0; i<2*num_sites; ++i)
					_H(i, i) += m_delta_chol;
			}
		}

		if(chol_try > 0)
		{
			using namespace tl2_ops;
			std::cerr << "Warning: Needed " << chol_try
				<< " corrections for cholesky decomposition for Q = "
				<< Qvec << "." << std::endl;
		}

		// see p. 5 in (Toth 2015)
		t_mat H_mat = C_mat * g_sign * tl2::herm<t_mat>(C_mat);

		bool is_herm = tl2::is_symm_or_herm<t_mat, t_real>(H_mat, m_eps);
		if(!is_herm)
		{
			using namespace tl2_ops;
			std::cerr << "Warning: Hamiltonian is not hermitian for Q = "
				<< Qvec << "." << std::endl;
		}

		// eigenvalues of the hamiltonian correspond to the energies
		// eigenvectors correspond to the spectral weights
		auto [evecs_ok, evals, evecs] =
			tl2_la::eigenvec<t_mat, t_vec, t_cplx, t_real>(
				H_mat, only_energies, is_herm, true);
		if(!evecs_ok)
		{
			using namespace tl2_ops;
			std::cerr << "Warning: Eigensystem calculation failed for Q = "
				<< Qvec << "." << std::endl;
		}


		std::vector<EnergyAndWeight> energies_and_correlations{};
		energies_and_correlations.reserve(evals.size());

		// register energies
		for(const auto& eval : evals)
		{
			EnergyAndWeight EandS { .E = eval.real(), };
			energies_and_correlations.emplace_back(std::move(EandS));
		}

		// weight factors
		if(!only_energies)
		{
			CalcCorrelationsFromHamiltonian(energies_and_correlations,
				H_mat, C_mat, g_sign, Qvec, evecs);
		}

		return energies_and_correlations;
	}


	/**
	 * get the dynamical structure factor from a hamiltonian
	 * @note implements the formalism given by (Toth 2015)
	 */
	void CalcCorrelationsFromHamiltonian(std::vector<EnergyAndWeight>& energies_and_correlations,
		const t_mat& H_mat, const t_mat& C_mat, const t_mat& g_sign,
		const t_vec_real& Qvec, const std::vector<t_vec>& evecs) const
	{
		const t_size num_sites = m_sites.size();
		if(num_sites == 0)
			return;

		// get the sorting of the energies
		std::vector<t_size> sorting = tl2::get_perm(
			energies_and_correlations.size(),
			[&energies_and_correlations](t_size idx1, t_size idx2) -> bool
		{
			return energies_and_correlations[idx1].E >=
				energies_and_correlations[idx2].E;
		});

		t_mat evec_mat = tl2::create<t_mat>(tl2::reorder(evecs, sorting));
		t_mat evec_mat_herm = tl2::herm(evec_mat);

		// equation (32) from (Toth 2015)
		t_mat L_mat = evec_mat_herm * H_mat * evec_mat; // energies
		t_mat E_sqrt = g_sign * L_mat;                  // abs. energies
		for(t_size i=0; i<E_sqrt.size1(); ++i)
			E_sqrt(i, i) = std::sqrt(E_sqrt/*L_mat*/(i, i)); // sqrt. of abs. energies

		// re-create energies, to be consistent with the weights
		energies_and_correlations.clear();
		for(t_size i=0; i<L_mat.size1(); ++i)
		{
			EnergyAndWeight EandS
			{
				.E = L_mat(i, i).real(),
				.S = tl2::zero<t_mat>(3, 3),
				.S_perp = tl2::zero<t_mat>(3, 3),
			};

			energies_and_correlations.emplace_back(std::move(EandS));
		}

		auto [C_inv, inv_ok] = tl2::inv(C_mat);
		if(!inv_ok)
		{
			using namespace tl2_ops;
			std::cerr << "Warning: Inversion failed for Q = "
				<< Qvec << "." << std::endl;
		}

		// equation (34) from (Toth 2015)
		t_mat trafo = C_inv * evec_mat * E_sqrt;
		t_mat trafo_herm = tl2::herm(trafo);

#ifdef __TLIBS2_MAGDYN_DEBUG_OUTPUT__
		t_mat D_mat = trafo_herm * H_mat * trafo;
		std::cout << "D = \n";
		tl2::niceprint(std::cout, D_mat, 1e-4, 4);
		std::cout << "\nE = \n";
		tl2::niceprint(std::cout, E_sqrt, 1e-4, 4);
		std::cout << "\nL = \n";
		tl2::niceprint(std::cout, L_mat, 1e-4, 4);
		std::cout << std::endl;
#endif

#ifdef __TLIBS2_MAGDYN_DEBUG_PY_OUTPUT__
		std::cout
			<< "# --------------------------------------------------------------------------------\n";
		std::cout << "Y = np.zeros(3*3*4*4, dtype=complex).reshape((4,4,3,3))" << std::endl;
		std::cout << "V = np.zeros(3*3*4*4, dtype=complex).reshape((4,4,3,3))" << std::endl;
		std::cout << "Z = np.zeros(3*3*4*4, dtype=complex).reshape((4,4,3,3))" << std::endl;
		std::cout << "W = np.zeros(3*3*4*4, dtype=complex).reshape((4,4,3,3))" << std::endl;
#endif

		// building the spin correlation functions of equation (47) from (Toth 2015)
		for(int x_idx=0; x_idx<3; ++x_idx)
		for(int y_idx=0; y_idx<3; ++y_idx)
		{
			// equations (44) from (Toth 2015)
			auto create_matrices = [num_sites](
				t_mat& V, t_mat& W, t_mat& Y, t_mat& Z)
			{
				V = tl2::create<t_mat>(num_sites, num_sites);
				W = tl2::create<t_mat>(num_sites, num_sites);
				Y = tl2::create<t_mat>(num_sites, num_sites);
				Z = tl2::create<t_mat>(num_sites, num_sites);
			};

			t_mat V, W, Y, Z;
			create_matrices(V, W, Y, Z);

			for(t_size i=0; i<num_sites; ++i)
			for(t_size j=0; j<num_sites; ++j)
			{
				auto calc_mat_elems = [this, i, j, x_idx, y_idx](
					const t_vec_real& Qvec,
					t_mat& Y, t_mat& V, t_mat& Z, t_mat& W)
				{
					// get the sites and spins
					const t_vec_real& pos_i = m_sites[i].pos;
					const t_vec_real& pos_j = m_sites[j].pos;
					t_real S_i = m_sites[i].spin_mag;
					t_real S_j = m_sites[j].spin_mag;

					// get the pre-calculated u vectors
					const t_vec& u_i = m_sites_calc[i].u;
					const t_vec& u_j = m_sites_calc[j].u;
					const t_vec& u_conj_i = m_sites_calc[i].u_conj;
					const t_vec& u_conj_j = m_sites_calc[j].u_conj;

					// pre-factors of equation (44) from (Toth 2015)
					t_real SiSj = 4. * std::sqrt(S_i*S_j);
					t_cplx phase = std::exp(-m_phase_sign * s_imag * s_twopi *
						tl2::inner<t_vec_real>(pos_j - pos_i, Qvec));

					// matrix elements of equation (44) from (Toth 2015)
					Y(i, j) = phase * SiSj * u_i[x_idx] * u_conj_j[y_idx];
					V(i, j) = phase * SiSj * u_conj_i[x_idx] * u_conj_j[y_idx];
					Z(i, j) = phase * SiSj * u_i[x_idx] * u_j[y_idx];
					W(i, j) = phase * SiSj * u_conj_i[x_idx] * u_j[y_idx];

#ifdef __TLIBS2_MAGDYN_DEBUG_PY_OUTPUT__
					std::cout << "Y[" << i << ", " << j << ", "
						<< x_idx << ", " << y_idx << "] = "
						<< Y(i, j).real() << " + " << Y(i, j).imag() << "j\n"
						<< "V[" << i << ", " << j << ", "
						<< x_idx << ", " << y_idx << "] = "
						<< V(i, j).real() << " + " << V(i, j).imag() << "j\n"
						<< "Z[" << i << ", " << j << ", "
						<< x_idx << ", " << y_idx << "] = "
						<< Z(i, j).real() << " + " << Z(i, j).imag() << "j\n"
						<< "W[" << i << ", " << j << ", "
						<< x_idx << ", " << y_idx << "] = "
						<< W(i, j).real() << " + " << W(i, j).imag() << "j"
						<< std::endl;
#endif
				};

				calc_mat_elems(Qvec, Y, V, Z, W);
			} // end of iteration over sites

			auto calc_S = [num_sites, x_idx, y_idx, &trafo, &trafo_herm, &energies_and_correlations]
				(t_mat EnergyAndWeight::*S, const t_mat& Y, const t_mat& V, const t_mat& Z, const t_mat& W)
			{
				// equation (47) from (Toth 2015)
				t_mat M = tl2::create<t_mat>(num_sites*2, num_sites*2);
				tl2::set_submat(M, Y, 0, 0);
				tl2::set_submat(M, V, num_sites, 0);
				tl2::set_submat(M, Z, 0, num_sites);
				tl2::set_submat(M, W, num_sites, num_sites);

				t_mat M_trafo = trafo_herm * M * trafo;

#ifdef __TLIBS2_MAGDYN_DEBUG_OUTPUT__
				std::cout << "M_trafo for x=" << x_idx << ", y=" << y_idx << ":\n";
				tl2::niceprint(std::cout, M_trafo, 1e-4, 4);
				std::cout << std::endl;
#endif

				for(t_size i=0; i<energies_and_correlations.size(); ++i)
					(energies_and_correlations[i].*S)(x_idx, y_idx) += M_trafo(i, i) / t_real(2*num_sites);
			};

			calc_S(&EnergyAndWeight::S, Y, V, Z, W);
		} // end of coordinate iteration

#ifdef __TLIBS2_MAGDYN_DEBUG_PY_OUTPUT__
			std::cout
				<< "# --------------------------------------------------------------------------------\n"
				<< std::endl;
#endif
	}


	/**
	 * applies projectors and weight factors to get neutron intensities
	 * @note implements the formalism given by (Toth 2015)
	 */
	void CalcIntensities(const t_vec_real& Qvec, std::vector<EnergyAndWeight>&
		energies_and_correlations) const
	{
		for(EnergyAndWeight& E_and_S : energies_and_correlations)
		{
			// apply bose factor
			if(m_temperature >= 0.)
				E_and_S.S *= tl2::bose_cutoff(E_and_S.E, m_temperature, m_bose_cutoff);

			// apply orthogonal projector for magnetic neutron scattering,
			// see (Shirane 2002), p. 37, equation (2.64)
			//t_vec bragg_rot = use_field ? m_rot_field * m_bragg : m_bragg;
			//proj_neutron = tl2::ortho_projector<t_mat, t_vec>(bragg_rot, false);
			t_mat proj_neutron = tl2::ortho_projector<t_mat, t_vec>(Qvec, false);
			E_and_S.S_perp = proj_neutron * E_and_S.S * proj_neutron;

			// weights
			E_and_S.weight = std::abs(tl2::trace<t_mat>(E_and_S.S_perp).real());
			E_and_S.weight_full = std::abs(tl2::trace<t_mat>(E_and_S.S).real());

			// polarisation channels
			for(int i=0; i<3; ++i)
			{
				t_mat pol = get_polarisation<t_mat>(i, false);
				t_mat Sperp = pol * E_and_S.S_perp;
				t_mat S = pol * E_and_S.S;

				E_and_S.weight_channel[i] = std::abs(tl2::trace<t_mat>(Sperp).real());
				E_and_S.weight_channel_full[i] = std::abs(tl2::trace<t_mat>(S).real());
			}
		}
	}


	/**
	 * unite degenerate energies and their corresponding eigenstates
	 */
	std::vector<EnergyAndWeight> UniteEnergies(const std::vector<EnergyAndWeight>&
		energies_and_correlations) const
	{
		std::vector<EnergyAndWeight> new_energies_and_correlations{};
		new_energies_and_correlations.reserve(energies_and_correlations.size());

		for(const auto& curState : energies_and_correlations)
		{
			t_real curE = curState.E;

			auto iter = std::find_if(
				new_energies_and_correlations.begin(),
				new_energies_and_correlations.end(),
				[curE, this](const auto& E_and_S) -> bool
			{
				t_real E = E_and_S.E;
				return tl2::equals<t_real>(E, curE, m_eps);
			});

			if(iter == new_energies_and_correlations.end())
			{
				// energy not yet seen
				new_energies_and_correlations.push_back(curState);
			}
			else
			{
				// energy already seen: add correlation matrices and weights
				iter->S += curState.S;
				iter->S_perp += curState.S_perp;

				iter->weight += curState.weight;
				iter->weight_full += curState.weight_full;

				for(int i=0; i<3; ++i)
				{
					iter->weight_channel[i] += curState.weight_channel[i];
					iter->weight_channel_full[i] += curState.weight_channel_full[i];
				}
			}
		}

		return new_energies_and_correlations;
	}


	/**
	 * get the energies and the spin-correlation at the given momentum
	 * (also calculates incommensurate contributions and applies weight factors)
	 * @note implements the formalism given by (Toth 2015)
	 */
	std::vector<EnergyAndWeight> CalcEnergies(const t_vec_real& Qvec,
		bool only_energies = false) const
	{
		std::vector<EnergyAndWeight> EandWs;
		if(m_calc_H)
		{
			t_mat H = CalcHamiltonian(Qvec);
			EandWs = CalcEnergiesFromHamiltonian(H, Qvec, only_energies);
		}

		if(IsIncommensurate())
		{
			// equations (39) and (40) from (Toth 2015)
			t_mat proj_norm = tl2::convert<t_mat>(
				tl2::projector<t_mat_real, t_vec_real>(
					m_rotaxis, true));

			t_mat rot_incomm = tl2::unit<t_mat>(3);
			rot_incomm -= s_imag * m_phase_sign * tl2::skewsymmetric<t_mat, t_vec>(m_rotaxis);
			rot_incomm -= proj_norm;
			rot_incomm *= 0.5;

			t_mat rot_incomm_conj = tl2::conj(rot_incomm);

			std::vector<EnergyAndWeight> EandWs_p, EandWs_m;

			if(m_calc_Hp)
			{
				t_mat H_p = CalcHamiltonian(Qvec + m_ordering);
				EandWs_p = CalcEnergiesFromHamiltonian(
					H_p, Qvec + m_ordering, only_energies);
			}

			if(m_calc_Hm)
			{
				t_mat H_m = CalcHamiltonian(Qvec - m_ordering);
				EandWs_m = CalcEnergiesFromHamiltonian(
					H_m, Qvec - m_ordering, only_energies);
			}

			if(!only_energies)
			{
				// formula 40 from (Toth 2015)
				for(EnergyAndWeight& EandW : EandWs)
					EandW.S = EandW.S * proj_norm;
				for(EnergyAndWeight& EandW : EandWs_p)
					EandW.S = EandW.S * rot_incomm;
				for(EnergyAndWeight& EandW : EandWs_m)
					EandW.S = EandW.S * rot_incomm_conj;
			}

			// unite energies and weights
			for(EnergyAndWeight& EandW : EandWs_p)
				EandWs.emplace_back(std::move(EandW));
			for(EnergyAndWeight& EandW : EandWs_m)
				EandWs.emplace_back(std::move(EandW));
		}

		if(!only_energies)
			CalcIntensities(Qvec, EandWs);

		if(m_unite_degenerate_energies)
			EandWs = UniteEnergies(EandWs);

		return EandWs;
	}


	std::vector<EnergyAndWeight> CalcEnergies(t_real h, t_real k, t_real l,
		bool only_energies = false) const
	{
		// momentum transfer
		const t_vec_real Qvec = tl2::create<t_vec_real>({ h, k, l });
		return CalcEnergies(Qvec, only_energies);
	}


	/**
	 * get the energy minimum
	 * @note a first version for a simplified ferromagnetic dispersion was based on (Heinsdorf 2021)
	 */
	t_real CalcMinimumEnergy() const
	{
		auto energies_and_correlations = CalcEnergies(0., 0., 0., true);
		auto min_iter = std::min_element(
			energies_and_correlations.begin(), energies_and_correlations.end(),
			[](const auto& E_and_S_1, const auto& E_and_S_2) -> bool
		{
			return std::abs(E_and_S_1.E) < std::abs(E_and_S_2.E);
		});

		if(min_iter == energies_and_correlations.end())
			return 0.;
		return min_iter->E;
	}


	/**
	 * get the ground-state energy
	 * @note zero-operator term in expansion of equation (20) in (Toth 2015)
	 */
	t_real CalcGroundStateEnergy() const
	{
		t_real E = 0.;

		for(const ExchangeTerm& term : m_exchange_terms)
		{
			t_mat J = CalcRealJ(term);  // Q=0 -> no rotation needed

			t_vec Si = m_sites[term.site1].spin_mag * m_sites_calc[term.site1].v;
			t_vec Sj = m_sites[term.site2].spin_mag * m_sites_calc[term.site2].v;

			E += tl2::inner_noconj<t_vec>(Si, J * Sj).real();
		}

		return E;
	}
	// --------------------------------------------------------------------


	/**
	 * generates the dispersion plot along the given q path
	 */
	void SaveDispersion(const std::string& filename,
		t_real h_start, t_real k_start, t_real l_start,
		t_real h_end, t_real k_end, t_real l_end,
		t_size num_qs = 128) const
	{
		std::ofstream ofstr{filename};
		SaveDispersion(ofstr, h_start, k_start, l_start, h_end, k_end, l_end, num_qs);
	}


	/**
	 * generates the dispersion plot along the given q path
	 */
	void SaveDispersion(std::ostream& ostr,
		t_real h_start, t_real k_start, t_real l_start,
		t_real h_end, t_real k_end, t_real l_end,
		t_size num_qs = 128) const
	{
		ostr.precision(m_prec);

		ostr
			<< std::setw(m_prec*2) << std::left << "# h"
			<< std::setw(m_prec*2) << std::left << "k"
			<< std::setw(m_prec*2) << std::left << "l"
			<< std::setw(m_prec*2) << std::left << "E"
			<< std::setw(m_prec*2) << std::left << "w"
			<< std::setw(m_prec*2) << std::left << "w_sf1"
			<< std::setw(m_prec*2) << std::left << "w_sf2"
			<< std::setw(m_prec*2) << std::left << "w_nsf"
			<< std::endl;

		for(t_size i=0; i<num_qs; ++i)
		{
			t_real h = std::lerp(h_start, h_end, t_real(i)/t_real(num_qs-1));
			t_real k = std::lerp(k_start, k_end, t_real(i)/t_real(num_qs-1));
			t_real l = std::lerp(l_start, l_end, t_real(i)/t_real(num_qs-1));

			auto energies_and_correlations = CalcEnergies(h, k, l, false);
			for(const auto& E_and_S : energies_and_correlations)
			{
				ostr
					<< std::setw(m_prec*2) << std::left << h
					<< std::setw(m_prec*2) << std::left << k
					<< std::setw(m_prec*2) << std::left << l
					<< std::setw(m_prec*2) << E_and_S.E
					<< std::setw(m_prec*2) << E_and_S.weight
					<< std::setw(m_prec*2) << E_and_S.weight_channel[0]
					<< std::setw(m_prec*2) << E_and_S.weight_channel[1]
					<< std::setw(m_prec*2) << E_and_S.weight_channel[2]
					<< std::endl;
			}
		}
	}


	/**
	 * load a configuration from a file
	 */
	bool Load(const std::string& filename)
	{
		// properties tree
		boost::property_tree::ptree node;

		// read xml file
		std::ifstream ifstr{filename};
		boost::property_tree::read_xml(ifstr, node);

		// check signature
		if(auto optInfo = node.get_optional<std::string>("magdyn.meta.info");
			!optInfo || !(*optInfo==std::string{"magdyn_tool"}))
		{
			return false;
		}

		const auto &magdyn = node.get_child("magdyn");
		return Load(magdyn);
	}


	/**
	 * save a configuration to a file
	 */
	bool Save(const std::string& filename) const
	{
		// properties tree
		boost::property_tree::ptree node;

		// write signature
		node.put<std::string>("meta.info", "magdyn_tool");
		node.put<std::string>("meta.date",
			tl2::epoch_to_str<t_real>(tl2::epoch<t_real>()));

		boost::property_tree::ptree root_node;
		root_node.put_child("magdyn", node);

		// write xml file
		std::ofstream ofstr{filename};
		if(!ofstr)
			return false;

		ofstr.precision(m_prec);
		boost::property_tree::write_xml(ofstr, node,
			boost::property_tree::xml_writer_make_settings(
				'\t', 1, std::string{"utf-8"}));
		return true;
	}


	/**
	 * load a configuration from a property tree
	 */
	bool Load(const boost::property_tree::ptree& node)
	{
		Clear();

		// variables
		if(auto vars = node.get_child_optional("variables"); vars)
		{
			for(const auto &var : *vars)
			{
				auto name = var.second.get_optional<std::string>("name");
				if(!name)
					continue;

				Variable variable;
				variable.name = *name;
				variable.value = var.second.get<t_cplx>("value", 0.);

				AddVariable(std::move(variable));
			}
		}

		// magnetic sites
		if(auto sites = node.get_child_optional("atom_sites"); sites)
		{
			for(const auto &site : *sites)
			{
				MagneticSite magnetic_site;

				magnetic_site.name = site.second.get<std::string>("name", "n/a");
				magnetic_site.index = m_sites.size();

				magnetic_site.pos = tl2::create<t_vec_real>(
				{
					site.second.get<t_real>("position_x", 0.),
					site.second.get<t_real>("position_y", 0.),
					site.second.get<t_real>("position_z", 0.),
				});

				magnetic_site.spin_dir[0] = site.second.get<std::string>("spin_x", "0");
				magnetic_site.spin_dir[1] = site.second.get<std::string>("spin_y", "0");
				magnetic_site.spin_dir[2] = site.second.get<std::string>("spin_z", "1");

				magnetic_site.spin_ortho[0] = site.second.get<std::string>("spin_ortho_x", "");
				magnetic_site.spin_ortho[1] = site.second.get<std::string>("spin_ortho_y", "");
				magnetic_site.spin_ortho[2] = site.second.get<std::string>("spin_ortho_z", "");

				magnetic_site.spin_mag = site.second.get<t_real>("spin_magnitude", 1.);
				magnetic_site.g = -2. * tl2::unit<t_mat>(3);

				AddMagneticSite(std::move(magnetic_site), false);
			}
		}

		// exchange terms / couplings
		if(auto terms = node.get_child_optional("exchange_terms"); terms)
		{
			for(const auto &term : *terms)
			{
				ExchangeTerm exchange_term;

				exchange_term.name = term.second.get<std::string>("name", "n/a");
				exchange_term.index = m_exchange_terms.size();
				exchange_term.site1 = term.second.get<t_size>("atom_1_index", 0);
				exchange_term.site2 = term.second.get<t_size>("atom_2_index", 0);

				// alternatively get the magnetic site indices via the names
				if(auto name1 = term.second.get_optional<std::string>("atom_1_name"); name1)
				{
					t_size site1_old = exchange_term.site1;
					if(auto sites1 = FindMagneticSites(*name1); sites1.size()==1)
						exchange_term.site1 = sites1[0]->index;
					if(exchange_term.site1 != site1_old)
					{
						std::cerr
							<< "Error in coupling " << exchange_term.index
							<< ": Index of site 1 (" << site1_old
							<< ") does not correspond to the selected name (" << *name1
							<< ")." << std::endl;
						exchange_term.site1 = site1_old;
					}
				}
				if(auto name2 = term.second.get_optional<std::string>("atom_2_name"); name2)
				{
					t_size site2_old = exchange_term.site2;
					if(auto sites2 = FindMagneticSites(*name2); sites2.size()==1)
						exchange_term.site2 = sites2[0]->index;
					if(exchange_term.site2 != site2_old)
					{
						std::cerr
							<< "Error in coupling " << exchange_term.index
							<< ": Index of site 2 (" << site2_old
							<< ") does not correspond to the selected name (" << *name2
							<< ")." << std::endl;
						exchange_term.site2 = site2_old;
					}
				}

				exchange_term.dist = tl2::create<t_vec_real>(
				{
					term.second.get<t_real>("distance_x", 0.),
					term.second.get<t_real>("distance_y", 0.),
					term.second.get<t_real>("distance_z", 0.),
				});

				exchange_term.J = term.second.get<std::string>("interaction", "0");

				static const std::array<std::string, 3> comps{{"x", "y", "z"}};
				for(t_size i=0; i<3; ++i)
				{
					exchange_term.dmi[i] = term.second.get<std::string>(
						std::string("dmi_") + comps[i], "0");

					for(t_size j=0; j<3; ++j)
					{
						exchange_term.Jgen[i][j] = term.second.get<std::string>(
							std::string("gen_") + comps[i] + comps[j], "0");
					}
				}

				AddExchangeTerm(std::move(exchange_term), false);
			}
		}

		// external field
		if(auto field = node.get_child_optional("field"); field)
		{
			ExternalField thefield;

			thefield.mag = 0.;
			thefield.align_spins = false;

			thefield.dir = tl2::create<t_vec_real>(
			{
				field->get<t_real>("direction_h", 0.),
				field->get<t_real>("direction_k", 0.),
				field->get<t_real>("direction_l", 1.),
			});

			if(auto optVal = field->get_optional<t_real>("magnitude"))
				thefield.mag = *optVal;

			if(auto optVal = field->get_optional<bool>("align_spins"))
				thefield.align_spins = *optVal;

			SetExternalField(thefield);
		}

		// temperature
		m_temperature = node.get<t_real>("temperature", -1.);

		// ordering vector
		if(auto ordering = node.get_child_optional("ordering"); ordering)
		{
			t_vec_real ordering_vec = tl2::create<t_vec_real>(
			{
				ordering->get<t_real>("h", 0.),
				ordering->get<t_real>("k", 0.),
				ordering->get<t_real>("l", 0.),
			});

			SetOrderingWavevector(ordering_vec);
		}

		// rotation axis
		if(auto axis = node.get_child_optional("rotation_axis"); axis)
		{
			t_vec_real rotaxis = tl2::create<t_vec_real>(
			{
				axis->get<t_real>("h", 1.),
				axis->get<t_real>("k", 0.),
				axis->get<t_real>("l", 0.),
			});

			SetRotationAxis(rotaxis);
		}

		CalcMagneticSites();
		CalcExchangeTerms();
		return true;
	}


	/**
	 * save a configuration to a property tree
	 */
	bool Save(boost::property_tree::ptree& node) const
	{
		// external field
		node.put<t_real>("field.direction_h", m_field.dir[0]);
		node.put<t_real>("field.direction_k", m_field.dir[1]);
		node.put<t_real>("field.direction_l", m_field.dir[2]);
		node.put<t_real>("field.magnitude", m_field.mag);
		node.put<bool>("field.align_spins", m_field.align_spins);

		// ordering vector
		if(m_ordering.size() == 3)
		{
			node.put<t_real>("ordering.h", m_ordering[0]);
			node.put<t_real>("ordering.k", m_ordering[1]);
			node.put<t_real>("ordering.l", m_ordering[2]);
		}

		// rotation axis
		if(m_rotaxis.size() == 3)
		{
			node.put<t_real>("rotation_axis.h", m_rotaxis[0]);
			node.put<t_real>("rotation_axis.k", m_rotaxis[1]);
			node.put<t_real>("rotation_axis.l", m_rotaxis[2]);
		}

		// temperature
		node.put<t_real>("temperature", m_temperature);

		// variables
		for(const auto& var : GetVariables())
		{
			boost::property_tree::ptree itemNode;
			itemNode.put<std::string>("name", var.name);
			itemNode.put<t_cplx>("value", var.value);

			node.add_child("variables.variable", itemNode);
		}

		// magnetic sites
		for(const auto& site : GetMagneticSites())
		{
			boost::property_tree::ptree itemNode;
			itemNode.put<std::string>("name", site.name);

			itemNode.put<t_real>("position_x", site.pos[0]);
			itemNode.put<t_real>("position_y", site.pos[1]);
			itemNode.put<t_real>("position_z", site.pos[2]);

			itemNode.put<std::string>("spin_x", site.spin_dir[0]);
			itemNode.put<std::string>("spin_y", site.spin_dir[1]);
			itemNode.put<std::string>("spin_z", site.spin_dir[2]);

			itemNode.put<std::string>("spin_ortho_x", site.spin_ortho[0]);
			itemNode.put<std::string>("spin_ortho_y", site.spin_ortho[1]);
			itemNode.put<std::string>("spin_ortho_z", site.spin_ortho[2]);

			itemNode.put<t_real>("spin_magnitude", site.spin_mag);

			node.add_child("atom_sites.site", itemNode);
		}

		// exchange terms
		for(const auto& term : GetExchangeTerms())
		{
			boost::property_tree::ptree itemNode;
			itemNode.put<std::string>("name", term.name);

			itemNode.put<t_size>("atom_1_index", term.site1);
			itemNode.put<t_size>("atom_2_index", term.site2);

			itemNode.put<t_real>("distance_x", term.dist[0]);
			itemNode.put<t_real>("distance_y", term.dist[1]);
			itemNode.put<t_real>("distance_z", term.dist[2]);

			itemNode.put<std::string>("interaction", term.J);

			static const std::array<std::string, 3> comps{{"x", "y", "z"}};
			for(t_size i=0; i<3; ++i)
			{
				itemNode.put<std::string>(std::string("dmi_") +
					comps[i], term.dmi[i]);

				for(t_size j=0; j<3; ++j)
				{
					itemNode.put<std::string>(std::string("gen_") +
						comps[i] + comps[j], term.Jgen[i][j]);
				}
			}

			// also save the magnetic site names
			const auto& sites = GetMagneticSites();
			if(term.site1 < sites.size())
				itemNode.put<std::string>("atom_1_name", sites[term.site1].name);
			if(term.site2 < sites.size())
				itemNode.put<std::string>("atom_2_name", sites[term.site2].name);

			node.add_child("exchange_terms.term", itemNode);
		}

		return true;
	}


protected:
	/**
	 * converts the rotation matrix rotating the local spins to ferromagnetic
	 * [001] directions into the vectors comprised of the matrix columns
	 * @see equation (9) and (51) from (Toth 2015)
	 */
	std::tuple<t_vec, t_vec> R_to_uv(const t_mat& R)
	{
		t_vec u = tl2::col<t_mat, t_vec>(R, 0)
			+ s_imag * tl2::col<t_mat, t_vec>(R, 1);
		t_vec v = tl2::col<t_mat, t_vec>(R, 2);

		return std::make_tuple(u, v);
	}


	/**
	 * rotate local spin to ferromagnetic [001] direction
	 * @see equations (7) and (9) from (Toth 2015)
	 */
	std::tuple<t_vec, t_vec> spin_to_uv(const t_vec& spin_dir)
	{
		auto [spin_re, spin_im] = tl2::split_cplx<t_vec, t_vec_real>(spin_dir);

		if(!tl2::equals_0<t_vec_real>(spin_im, m_eps))
			std::cerr << "Warning: Spin vector should be purely real." << std::endl;

		// only use real part, imaginary part should be zero
		//spin_re /= tl2::norm<t_vec_real>(spin_re);
		t_mat_real _rot = tl2::rotation<t_mat_real, t_vec_real>(
			spin_re, m_zdir, &m_rotaxis, m_eps);

		t_mat rot = tl2::convert<t_mat, t_mat_real>(_rot);
		return R_to_uv(rot);
	}


private:
	// magnetic sites
	std::vector<MagneticSite> m_sites{};
	std::vector<MagneticSiteCalc> m_sites_calc{};

	// magnetic couplings
	std::vector<ExchangeTerm> m_exchange_terms{};
	std::vector<ExchangeTermCalc> m_exchange_terms_calc{};

	// open variables
	std::vector<Variable> m_variables{};

	// external field
	ExternalField m_field{};
	// matrix to rotate field into the [001] direction
	t_mat m_rot_field = tl2::unit<t_mat>(3);

	// ordering wave vector for incommensurate structures
	t_vec_real m_ordering = tl2::zero<t_vec_real>(3);
	t_vec_real m_rotaxis = tl2::create<t_vec_real>({1., 0., 0.});

	// calculate the hamiltonian for Q, Q+ordering, and Q-ordering
	bool m_calc_H{true};
	bool m_calc_Hp{true};
	bool m_calc_Hm{true};

	// direction to rotation spins into, usually [001]
	t_vec_real m_zdir = tl2::create<t_vec_real>({0., 0., 1.});

	// temperature (-1: disable bose factor)
	t_real m_temperature{-1};

	// bose cutoff energy to avoid infinities
	t_real m_bose_cutoff{0.025};

	// settings
	bool m_is_incommensurate{false};
	bool m_force_incommensurate{false};
	bool m_unite_degenerate_energies{true};

	// settings for cholesky decomposition
	t_size m_tries_chol{50};
	t_real m_delta_chol{0.0025};

	// precisions
	t_real m_eps{1e-6};
	int m_prec{6};

	// conventions
	t_real m_phase_sign{-1.};

	// constants
	static constexpr const t_cplx s_imag {t_real(0), t_real(1)};
	static constexpr const t_real s_twopi {t_real(2)*tl2::pi<t_real>};
};

}
#endif
