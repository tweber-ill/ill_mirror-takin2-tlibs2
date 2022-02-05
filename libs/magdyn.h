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
 *   - (Heinsdorf 2021) N. Heinsdorf, manual example calculation for a simple
 *                      ferromagnetic case, personal communications, 2021/2022.
 *
 * This file implements the formalism given by (Toth 2015).
 */

#ifndef __TLIBS2_MAGDYN_H__
#define __TLIBS2_MAGDYN_H__

#include <vector>
#include <tuple>
#include <string>

#include <algorithm>
#include <numeric>

#include <iostream>
#include <fstream>
#include <iomanip>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#define USE_LAPACK 1
#include "tlibs2/libs/maths.h"
#include "units.h"
#include "phys.h"
#include "helper.h"


namespace tl2_mag {

// ----------------------------------------------------------------------------
// data types
// ----------------------------------------------------------------------------
using t_size = std::size_t;
using t_real = double;
using t_cplx = std::complex<t_real>;

using t_vec_real = tl2::vec<t_real, std::vector>;
using t_mat_real = tl2::mat<t_real, std::vector>;

using t_vec = tl2::vec<t_cplx, std::vector>;
using t_mat = tl2::mat<t_cplx, std::vector>;
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// input- and output structs
// ----------------------------------------------------------------------------
/**
 * magnetic atom sites
 */
struct AtomSite
{
	std::string name{}; // identifier
	t_size index{};     // index

	t_vec_real pos{};   // atom position

	t_vec spin_dir{};   // spin direction
	t_real spin_mag{};  // spin magnitude
	t_mat g{};          // g factor
};


/**
 * storage for temporary per-site results
 */
struct AtomSiteCalc
{
	t_vec u{}, u_conj{};
	t_vec v{};
};


/**
 * couplings between magnetic atoms
 */
struct ExchangeTerm
{
	std::string name{}; // identifier
	t_size index{};     // index

	t_size atom1{};     // atom 1 index
	t_size atom2{};     // atom 2 index
	t_vec_real dist{};  // distance between unit cells

	t_cplx J{};         // Heisenberg interaction
	t_vec dmi{};        // Dzyaloshinskij-Moriya interaction
};


/**
 * terms related to an external magnetic field
 */
struct ExternalField
{
	t_vec_real dir{};   // field direction
	t_real mag{};       // field magnitude
	bool align_spins{}; // align spins along external field
};


/**
 * eigenenergies and spin-spin correlation matrix
 */
struct EnergyAndWeight
{
	t_real E{};
	t_mat S{};
	t_mat S_p{}, S_m;
	t_mat S_perp{};
	t_real weight{};
	t_real weight_spinflip[2] = {0., 0.};
	t_real weight_nonspinflip{};
};
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// helper functions
// ----------------------------------------------------------------------------
/**
 * converts the rotation matrix rotating the local spins to ferromagnetic
 * [001] directions into the vectors comprised of the matrix columns
 * @see equation (9) from (Toth 2015).
 */
template<class t_mat, class t_vec, class t_cplx>
requires tl2::is_mat<t_mat> && tl2::is_vec<t_vec>
std::tuple<t_vec, t_vec> R_to_uv(const t_mat& R)
{
	// imaginary unit
	const t_cplx imag{0., 1.};

	t_vec u = tl2::col<t_mat, t_vec>(R, 0)
		+ imag*tl2::col<t_mat, t_vec>(R, 1);
	t_vec v = tl2::col<t_mat, t_vec>(R, 2);

	return std::make_tuple(u, v);
}


/**
 * print a matrix
 */
template<class t_mat> requires tl2::is_mat<t_mat>
void dbg_print(const t_mat& mat)
{
	using t_elem = typename t_mat::value_type;

	t_real eps = 1e-6;
	int prec = 3;
	std::cout.precision(prec);

	for(t_size row=0; row<mat.size1(); ++row)
	{
		for(t_size col=0; col<mat.size2(); ++col)
		{
			t_elem elem = mat(row, col);
			tl2::set_eps_0(elem, eps);
			std::cout << std::setw(prec*3) << elem;
		}

		std::cout << std::endl;
	}
}


/**
 * print a vector
 */
template<class t_vec> requires tl2::is_vec<t_vec>
void dbg_print(const t_vec& vec)
{
	using t_elem = typename t_vec::value_type;

	t_real eps = 1e-6;
	int prec = 3;
	std::cout.precision(prec);

	for(t_size row=0; row<vec.size(); ++row)
	{
		t_elem elem = vec[row];
		tl2::set_eps_0(elem, eps);
		std::cout << elem << std::endl;
	}
}
// ----------------------------------------------------------------------------



/**
 * calculates magnon dynamics,
 * implementing the formalism given in (Toth 2015)
 */
class MagDyn
{
public:
	MagDyn() = default;
	~MagDyn() = default;


	// --------------------------------------------------------------------
	// cleanup functions
	// --------------------------------------------------------------------
	/**
	 * clear all
	 */
	void Clear()
	{
		ClearAtomSites();
		ClearExchangeTerms();
		ClearExternalField();

		// clear temperature, -1: don't use
		m_temperature = -1.;

		// clear ordering wave vector
		m_ordering = tl2::zero<t_vec_real>(3);
	}


	/**
	 * clear all atom sites
	 */
	void ClearAtomSites()
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
	const std::vector<AtomSite>& GetAtomSites() const
	{
		return m_sites;
	}


	const std::vector<ExchangeTerm>& GetExchangeTerms() const
	{
		return m_exchange_terms;
	}


	const ExternalField& GetExternalField() const
	{
		return m_field;
	}


	const t_vec& GetBraggPeak() const
	{
		return m_bragg;
	}


	t_real GetTemperature() const
	{
		return m_temperature;
	}


	t_real GetBoseCutoffEnergy() const
	{
		return m_bose_cutoff;
	}
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	// setter
	// --------------------------------------------------------------------
	void SetEpsilon(t_real eps)
	{
		m_eps = eps;
	}


	void SetPrecision(int prec)
	{
		m_prec = prec;
	}


	void SetExternalField(const ExternalField& field)
	{
		m_field = field;
	}


	void SetOrderingWavevector(const t_vec_real& ordering)
	{
		m_ordering = ordering;
		m_is_incommensurate = !tl2::equals_0<t_vec_real>(m_ordering, m_eps);
	}


	const t_vec_real& GetOrderingWavevector() const
	{
		return m_ordering;
	}


	void SetRotationAxis(const t_vec_real& axis)
	{
		m_rotaxis = axis;
	}


	void AddAtomSite(AtomSite&& site)
	{
		site.index = GetAtomSites().size();
		m_sites.emplace_back(
			std::forward<AtomSite&&>(site));
	}


	void AddExchangeTerm(ExchangeTerm&& term)
	{
		term.index = GetExchangeTerms().size();
		m_exchange_terms.emplace_back(
			std::forward<ExchangeTerm&&>(term));
	}


	void AddExchangeTerm(t_size atom1, t_size atom2,
		const t_vec_real& dist, const t_cplx& J)
	{
		ExchangeTerm term
		{
			.atom1 = atom1, // index of first atom
			.atom2 = atom2, // index of second atom
			.dist = dist,   // distance between the atom's unit cells (not the atoms)
			.J = J,         // interaction
		};

		AddExchangeTerm(std::move(term));
	}


	void SetBraggPeak(t_real h, t_real k, t_real l)
	{
		m_bragg = tl2::create<t_vec>({h, k, l});

		// call CalcSpinRotation() afterwards to calculate projector
	}


	void SetTemperature(t_real T)
	{
		m_temperature = T;
	}


	void SetBoseCutoffEnergy(t_real E)
	{
		m_bose_cutoff = E;
	}


	void SetUniteDegenerateEnergies(bool b)
	{
		m_unite_degenerate_energies = b;
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


	/**
	 * calculate the spin rotation trafo
	 */
	void CalcSpinRotation()
	{
		const t_size num_sites = m_sites.size();
		if(num_sites == 0)
			return;

		const t_vec_real zdir = tl2::create<t_vec_real>({0., 0., 1.});

		m_sites_calc.clear();
		m_sites_calc.reserve(num_sites);

		bool use_field = !tl2::equals_0<t_real>(m_field.mag, m_eps)
			&& m_field.dir.size() == 3;

		if(use_field)
		{
			// rotate field to [001] direction
			m_rot_field = tl2::convert<t_mat>(
				tl2::rotation<t_mat_real, t_vec_real>(
					-m_field.dir, zdir));
		}

		if(m_bragg.size() == 3)
		{
			// calculate orthogonal projector for magnetic neutron scattering
			// see (Shirane 2002), p. 37, eq. (2.64)
			t_vec bragg_rot = use_field ? m_rot_field * m_bragg : m_bragg;

			m_proj_neutron = tl2::ortho_projector<t_mat, t_vec>(
				bragg_rot, false);
		}
		else
		{
			// no bragg peak given -> don't project
			m_proj_neutron = tl2::unit<t_mat>(3);
		}


		for(t_size site_idx=0; site_idx<m_sites.size(); ++site_idx)
		{
			const AtomSite& site = m_sites[site_idx];

			// rotate local spin to ferromagnetic [001] direction
			auto [spin_re, spin_im] =
				tl2::split_cplx<t_vec, t_vec_real>(site.spin_dir);
			t_mat rot = tl2::convert<t_mat>(
				tl2::rotation<t_mat_real, t_vec_real>(spin_re, zdir));

			// spin rotation of equation (9) from (Toth 2015)
			t_vec u, v;

			if(m_field.align_spins)
				std::tie(u, v) = R_to_uv<t_mat, t_vec, t_cplx>(m_rot_field);
			else
				std::tie(u, v) = R_to_uv<t_mat, t_vec, t_cplx>(rot);

			AtomSiteCalc site_calc{};
			site_calc.u_conj = std::move(tl2::conj(u));
			site_calc.u = std::move(u);
			site_calc.v = std::move(v);
			m_sites_calc.emplace_back(std::move(site_calc));
		}
	}


	/**
	 * update the site and term indices
	 */
	void CalcIndices()
	{
		for(t_size site_idx=0; site_idx<m_sites.size(); ++site_idx)
		{
			AtomSite& site = m_sites[site_idx];
			site.index = site_idx;
		}

		for(t_size term_idx=0; term_idx<m_exchange_terms.size(); ++term_idx)
		{
			ExchangeTerm& term = m_exchange_terms[term_idx];
			term.index = term_idx;
		}
	}


	/**
	 * get the hamiltonian at the given momentum
	 * (CalcSpinRotation() needs to be called once before this function)
	 * @note implements the formalism given by (Toth 2015)
	 * @note a first version for a simplified ferromagnetic dispersion was based on (Heinsdorf 2021)
	 */
	t_mat GetHamiltonian(t_real h, t_real k, t_real l) const
	{
		const t_size num_sites = m_sites.size();
		if(num_sites == 0)
			return {};

		// momentum
		const t_vec_real Q = tl2::create<t_vec_real>({h, k, l});

		// constants: imaginary unit and 2pi
		constexpr const t_cplx imag{0., 1.};
		constexpr const t_real twopi = t_real(2)*tl2::pi<t_real>;
		// bohr magneton in [meV/T]
		constexpr const t_real muB = tl2::muB<t_real>
			/ tl2::meV<t_real> * tl2::tesla<t_real>;

		// build the interaction matrices J(Q) and J(-Q) of
		// equations (12) and (14) from (Toth 2015)
		t_mat J_Q = tl2::zero<t_mat>(num_sites*3, num_sites*3);
		t_mat J_mQ = tl2::zero<t_mat>(num_sites*3, num_sites*3);
		t_mat J_Q0 = tl2::zero<t_mat>(num_sites*3, num_sites*3);
		for(const ExchangeTerm& term : m_exchange_terms)
		{
			if(term.atom1 >= num_sites || term.atom2 >= num_sites)
				continue;

			// exchange interaction matrix with dmi as anti-symmetric part,
			// see (Toth 2015) p. 2
			t_mat J = tl2::diag<t_mat>(
				tl2::create<t_vec>({term.J, term.J, term.J}));

			if(term.dmi.size() == 3)
			{
				// cross product matrix
				J += tl2::skewsymmetric<t_mat, t_vec>(-term.dmi);
			}

			// incommensurate case: rotation wrt magnetic unit cell
			// equations (21), (6) and (2) from (Toth 2015)
			if(m_is_incommensurate)
			{
				t_real rot_UC_angle = tl2::inner<t_vec_real>(m_ordering, term.dist);
				if(!tl2::equals_0<t_real>(rot_UC_angle, m_eps))
				{
					t_mat rot_UC = tl2::convert<t_mat>(
						 tl2::rotation<t_mat_real, t_vec_real>(
							m_rotaxis, rot_UC_angle));

					J = J * rot_UC;
				}
			}

			// equation (14) from (Toth 2015)
			t_cplx phase_Q = std::exp(-imag * twopi *
				tl2::inner<t_vec_real>(term.dist, Q));
			t_cplx phase_mQ = std::exp(-imag * twopi *
				tl2::inner<t_vec_real>(-term.dist, Q));

			t_mat J_T = tl2::trans(J);
			t_real factor = 1.; //0.5;

			// include these two terms to fulfill
			// equation (11) from (Toth 2015)
			tl2::add_submat<t_mat>(J_Q, factor * J * phase_Q,
				term.atom1*3, term.atom2*3);
			tl2::add_submat<t_mat>(J_Q, factor * J_T * phase_mQ,
				term.atom2*3, term.atom1*3);

			tl2::add_submat<t_mat>(J_mQ, factor * J * phase_mQ,
				term.atom1*3, term.atom2*3);
			tl2::add_submat<t_mat>(J_mQ, factor * J_T * phase_Q,
				term.atom2*3, term.atom1*3);

			tl2::add_submat<t_mat>(J_Q0, factor*J,
				term.atom1*3, term.atom2*3);
			tl2::add_submat<t_mat>(J_Q0, factor*J_T,
				term.atom2*3, term.atom1*3);
		}


		// create the hamiltonian of equation (25) and (26) from (Toth 2015)
		t_mat A = tl2::create<t_mat>(num_sites, num_sites);
		t_mat A_mQ = tl2::create<t_mat>(num_sites, num_sites);
		t_mat B = tl2::create<t_mat>(num_sites, num_sites);
		t_mat C = tl2::zero<t_mat>(num_sites, num_sites);

		bool use_field = !tl2::equals_0<t_real>(m_field.mag, m_eps)
			&& m_field.dir.size() == 3;

		for(t_size i=0; i<num_sites; ++i)
		for(t_size j=0; j<num_sites; ++j)
		{
			t_mat J_sub_Q = tl2::submat<t_mat>(J_Q, i*3, j*3, 3, 3);
			t_mat J_sub_mQ = tl2::submat<t_mat>(J_mQ, i*3, j*3, 3, 3);

			// TODO: check units of S_i and S_j
			t_real S_i = m_sites[i].spin_mag;
			t_real S_j = m_sites[j].spin_mag;
			const t_vec& u_i = m_sites_calc[i].u;
			const t_vec& u_j = m_sites_calc[j].u;
			const t_vec& u_conj_j = m_sites_calc[j].u_conj;
			const t_vec& v_i = m_sites_calc[i].v;

			t_real factor = 0.5 * std::sqrt(S_i*S_j);
			A(i, j) = factor *
				tl2::inner_noconj<t_vec>(u_i, J_sub_Q * u_conj_j);
			A_mQ(i, j) = factor *
				tl2::inner_noconj<t_vec>(u_i, J_sub_mQ * u_conj_j);
			B(i, j) = factor *
				tl2::inner_noconj<t_vec>(u_i, J_sub_Q * u_j);

			if(i == j)
			{
				for(t_size k=0; k<num_sites; ++k)
				{
					// TODO: check unit of S_k
					t_real S_k = m_sites[k].spin_mag;
					const t_vec& v_k = m_sites_calc[k].v;

					t_mat J_sub_Q0 = tl2::submat<t_mat>(J_Q0, i*3, k*3, 3, 3);
					C(i, j) += S_k * tl2::inner_noconj<t_vec>(v_i, J_sub_Q0 * v_k);
				}
			}

			// include external field, equation (28) from (Toth 2015)
			if(use_field && i == j)
			{
				t_vec B = tl2::convert<t_vec>(
					m_field.dir / tl2::norm<t_vec_real>(m_field.dir));
				B = B * m_field.mag;

				t_vec gv = m_sites[i].g * v_i;
				t_cplx Bgv = tl2::inner_noconj<t_vec>(B, gv);

				A(i, j) -= 0.5*muB * Bgv;
				A_mQ(i, j) -= 0.5*muB * Bgv;
			}
		}


		// test matrix block
		//return A_conj - C;

		t_mat H = tl2::zero<t_mat>(num_sites*2, num_sites*2);
		tl2::set_submat(H, A - C, 0, 0);
		tl2::set_submat(H, B, 0, num_sites);
		tl2::set_submat(H, tl2::herm(B), num_sites, 0);
		tl2::set_submat(H, tl2::conj(A_mQ) - C, num_sites, num_sites);

		return H;
	}


	/**
	 * get the energies and the spin-correlation at the given momentum
	 * @note implements the formalism given by (Toth 2015)
	 */
	std::vector<EnergyAndWeight> GetEnergies(
		t_mat _H, t_real h, t_real k, t_real l,
		bool only_energies = false) const
	{
		const t_size num_sites = m_sites.size();
		if(num_sites == 0 || _H.size1() == 0)
			return {};

		// constants: imaginary unit and 2pi
		constexpr const t_cplx imag{0., 1.};
		constexpr const t_real twopi = t_real(2)*tl2::pi<t_real>;

		// equation (30) from (Toth 2015)
		t_mat g = tl2::zero<t_mat>(num_sites*2, num_sites*2);
		for(t_size i=0; i<num_sites; ++i)
			g(i, i) = 1.;
		for(t_size i=num_sites; i<2*num_sites; ++i)
			g(i, i) = -1.;

		// equation (31) from (Toth 2015)
		t_mat C;
		for(t_size retry=0; retry<m_retries_chol; ++retry)
		{
			auto [chol_ok, _C] = tl2_la::chol<t_mat>(_H);

			if(chol_ok)
			{
				C = _C;
				break;
			}
			else
			{
				// try forcing the hamilton to be positive definite
				for(t_size i=0; i<2*num_sites; ++i)
					_H(i, i) += m_eps_chol;
			}

			if(!chol_ok && retry == m_retries_chol-1)
			{
				std::cerr << "Warning: Cholesky decomposition failed"
					<< " for Q = (" << h << ", " << k << ", " << l << ")"
					<< "." << std::endl;
				C = _C;
			}
		}

		t_mat C_herm = tl2::herm<t_mat>(C);

		// see p. 5 in (Toth 2015)
		t_mat H = C * g * C_herm;
		//dbg_print(_H);
		//dbg_print(H);

		bool is_herm = tl2::is_symm_or_herm<t_mat, t_real>(H, m_eps);
		if(!is_herm)
		{
			std::cerr << "Warning: Hamiltonian is not hermitian"
				<< " for Q = (" << h << ", " << k << ", " << l << ")"
				<< "." << std::endl;
		}

		// eigenvalues of the hamiltonian correspond to the energies
		// eigenvectors correspond to the spectral weights
		auto [evecs_ok, evals, evecs] =
			tl2_la::eigenvec<t_mat, t_vec, t_cplx, t_real>(
				H, only_energies, is_herm, true);
		if(!evecs_ok)
		{
			std::cerr << "Warning: Eigensystem calculation failed"
				<< " for Q = (" << h << ", " << k << ", " << l << ")"
				<< "." << std::endl;
		}


		std::vector<EnergyAndWeight> energies_and_correlations{};
		energies_and_correlations.reserve(evals.size());

		// register energies
		for(const auto& eval : evals)
		{
			EnergyAndWeight EandS
			{
				.E = eval.real(),
			};

			energies_and_correlations.emplace_back(std::move(EandS));
		}


		// weight factors
		if(!only_energies)
		{
			// get the sorting of the energies
			std::vector<t_size> sorting(energies_and_correlations.size());
			std::iota(sorting.begin(), sorting.end(), 0);

			std::stable_sort(sorting.begin(), sorting.end(),
				[&energies_and_correlations](t_size idx1, t_size idx2) -> bool
				{
					return energies_and_correlations[idx1].E >=
						energies_and_correlations[idx2].E;
				});

			//for(t_size idx=0; idx<sorting.size(); ++idx)
			//	std::cout << idx << " -> " << sorting[idx] << std::endl;

			//energies_and_correlations = tl2::reorder(energies_and_correlations, sorting);
			evecs = tl2::reorder(evecs, sorting);
			evals = tl2::reorder(evals, sorting);

			/*for(t_vec& evec : evecs)
			{
				dbg_print<t_vec>(evec);
				std::cout << std::endl;
			}*/

			t_mat evec_mat = tl2::create<t_mat>(evecs);
			t_mat evec_mat_herm = tl2::herm(evec_mat);
			//dbg_print(evec_mat);
			//std::cout << std::endl;

			// equation (32) from (Toth 2015)
			t_mat L = evec_mat_herm * H * evec_mat;
			t_mat E = g*L;

			// re-create energies, to be consistent with the weights
			energies_and_correlations.clear();
			for(t_size i=0; i<L.size1(); ++i)
			{
				EnergyAndWeight EandS
				{
					.E = L(i,i).real(),
					.S = tl2::zero<t_mat>(3, 3),
					.S_p = tl2::zero<t_mat>(3, 3),
					.S_m = tl2::zero<t_mat>(3, 3),
					.S_perp = tl2::zero<t_mat>(3, 3),
				};

				energies_and_correlations.emplace_back(
					std::move(EandS));
			}

			t_mat E_sqrt = E;
			for(t_size i=0; i<E.size1(); ++i)
			{
				t_real energy = E_sqrt(i, i).real();
				E_sqrt(i, i) = std::sqrt(std::abs(energy));
			}

			auto [C_inv, inv_ok] = tl2::inv(C);
			if(!inv_ok)
			{
				std::cerr << "Warning: Inversion failed"
					<< " for Q = (" << h << ", " << k << ", " << l << ")"
					<< "." << std::endl;
			}

			// equation (34) from (Toth 2015)
			t_mat trafo = C_inv * evec_mat * E_sqrt;
			t_mat trafo_herm = tl2::herm(trafo);

			//t_mat D = trafo_herm * _H * trafo;
			//dbg_print(D);
			//dbg_print(E);
			//dbg_print(L);
			//std::cout << std::endl;

			// momentum
			const t_vec_real Q = tl2::create<t_vec_real>({h, k, l});

			/*std::cout << "Y = np.zeros(3*3*4*4, dtype=complex).reshape((4,4,3,3))" << std::endl;
			std::cout << "V = np.zeros(3*3*4*4, dtype=complex).reshape((4,4,3,3))" << std::endl;
			std::cout << "Z = np.zeros(3*3*4*4, dtype=complex).reshape((4,4,3,3))" << std::endl;
			std::cout << "W = np.zeros(3*3*4*4, dtype=complex).reshape((4,4,3,3))" << std::endl;*/

			// building the spin correlation functions of equation (47) from (Toth 2015)
			for(int x_idx=0; x_idx<3; ++x_idx)
			for(int y_idx=0; y_idx<3; ++y_idx)
			{
				// equations (44) from (Toth 2015)
				t_mat V = tl2::create<t_mat>(num_sites, num_sites);
				t_mat W = tl2::create<t_mat>(num_sites, num_sites);
				t_mat Y = tl2::create<t_mat>(num_sites, num_sites);
				t_mat Z = tl2::create<t_mat>(num_sites, num_sites);

				// incommensurate case
				t_mat V_p, W_p, Y_p, Z_p;
				t_mat V_m, W_m, Y_m, Z_m;
				if(m_is_incommensurate)
				{
					V_p = tl2::create<t_mat>(num_sites, num_sites);
					W_p = tl2::create<t_mat>(num_sites, num_sites);
					Y_p = tl2::create<t_mat>(num_sites, num_sites);
					Z_p = tl2::create<t_mat>(num_sites, num_sites);

					V_m = tl2::create<t_mat>(num_sites, num_sites);
					W_m = tl2::create<t_mat>(num_sites, num_sites);
					Y_m = tl2::create<t_mat>(num_sites, num_sites);
					Z_m = tl2::create<t_mat>(num_sites, num_sites);
				}

				for(t_size i=0; i<num_sites; ++i)
				for(t_size j=0; j<num_sites; ++j)
				{
					const t_vec_real& pos_i = m_sites[i].pos;
					const t_vec_real& pos_j = m_sites[j].pos;

					t_real S_i = m_sites[i].spin_mag;
					t_real S_j = m_sites[j].spin_mag;

					const t_vec& u_i = m_sites_calc[i].u;
					const t_vec& u_j = m_sites_calc[j].u;
					const t_vec& u_conj_i = m_sites_calc[i].u_conj;
					const t_vec& u_conj_j = m_sites_calc[j].u_conj;

					// TODO: check units of S
					t_real factor = 4. * std::sqrt(S_i*S_j);

					t_cplx phase = std::exp(+imag * twopi *
						tl2::inner<t_vec_real>(pos_j - pos_i, Q));

					// TODO: check phase factors
					Y(i, j) = factor * phase * u_i[x_idx] * u_conj_j[y_idx];
					V(i, j) = factor * phase * u_conj_i[x_idx] * u_conj_j[y_idx];
					Z(i, j) = factor * phase * u_i[x_idx] * u_j[y_idx];
					W(i, j) = factor * phase * u_conj_i[x_idx] * u_j[y_idx];

					// incommensurate case
					if(m_is_incommensurate)
					{
						t_cplx phase_p = std::exp(+imag * twopi *
							tl2::inner<t_vec_real>(pos_j - pos_i,
								Q + m_ordering));

						t_cplx phase_m = std::exp(+imag * twopi *
							tl2::inner<t_vec_real>(pos_j - pos_i,
								Q - m_ordering));

						Y_p(i, j) = factor * phase_p * u_i[x_idx] * u_conj_j[y_idx];
						V_p(i, j) = factor * phase_p * u_conj_i[x_idx] * u_conj_j[y_idx];
						Z_p(i, j) = factor * phase_p * u_i[x_idx] * u_j[y_idx];
						W_p(i, j) = factor * phase_p * u_conj_i[x_idx] * u_j[y_idx];

						Y_m(i, j) = factor * phase_m * u_i[x_idx] * u_conj_j[y_idx];
						V_m(i, j) = factor * phase_m * u_conj_i[x_idx] * u_conj_j[y_idx];
						Z_m(i, j) = factor * phase_m * u_i[x_idx] * u_j[y_idx];
						W_m(i, j) = factor * phase_m * u_conj_i[x_idx] * u_j[y_idx];
					}

					/*std::cout
						<< "Y[" << i << ", " << j << ", "
						<< x_idx << ", " << y_idx << "] = "
						<< Y(i, j).real() << " + " << Y(i, j).imag() << "j"
						<< std::endl;
					std::cout
						<< "V[" << i << ", " << j << ", "
						<< x_idx << ", " << y_idx << "] = "
						<< V(i, j).real() << " + " << V(i, j).imag() << "j"
						<< std::endl;
					std::cout
						<< "Z[" << i << ", " << j << ", "
						<< x_idx << ", " << y_idx << "] = "
						<< Z(i, j).real() << " + " << Z(i, j).imag() << "j"
						<< std::endl;
					std::cout
						<< "W[" << i << ", " << j << ", "
						<< x_idx << ", " << y_idx << "] = "
						<< W(i, j).real() << " + " << W(i, j).imag() << "j"
						<< std::endl;*/
				} // end of iteration over sites

				// equation (47) from (Toth 2015)
				t_mat M = tl2::create<t_mat>(num_sites*2, num_sites*2);
				tl2::set_submat(M, Y, 0, 0);
				tl2::set_submat(M, V, num_sites, 0);
				tl2::set_submat(M, Z, 0, num_sites);
				tl2::set_submat(M, W, num_sites, num_sites);

				t_mat M_trafo = trafo_herm * M * trafo;

				// incommensurate case
				t_mat M_p, M_m;
				t_mat M_trafo_p, M_trafo_m;
				if(m_is_incommensurate)
				{
					M_p = tl2::create<t_mat>(num_sites*2, num_sites*2);
					tl2::set_submat(M_p, Y_p, 0, 0);
					tl2::set_submat(M_p, V_p, num_sites, 0);
					tl2::set_submat(M_p, Z_p, 0, num_sites);
					tl2::set_submat(M_p, W_p, num_sites, num_sites);

					M_m = tl2::create<t_mat>(num_sites*2, num_sites*2);
					tl2::set_submat(M_m, Y_m, 0, 0);
					tl2::set_submat(M_m, V_m, num_sites, 0);
					tl2::set_submat(M_m, Z_m, 0, num_sites);
					tl2::set_submat(M_m, W_m, num_sites, num_sites);

					M_trafo_p = trafo_herm * M_p * trafo;
					M_trafo_m = trafo_herm * M_m * trafo;
				}

				for(t_size i=0; i<num_sites*2; ++i)
				{
					t_mat& S = energies_and_correlations[i].S;
					S(x_idx, y_idx) += M_trafo(i, i) / t_real(2*num_sites);

					if(m_is_incommensurate)
					{
						t_mat& S_p = energies_and_correlations[i].S_m;
						t_mat& S_m = energies_and_correlations[i].S_p;

						S_p(x_idx, y_idx) += M_trafo_p(i, i) / t_real(2*num_sites);
						S_m(x_idx, y_idx) += M_trafo_m(i, i) / t_real(2*num_sites);
					}
				}

				/*using namespace tl2_ops;
				tl2::set_eps_0<t_mat, t_real>(V, m_eps);
				tl2::set_eps_0<t_mat, t_real>(W, m_eps);
				tl2::set_eps_0<t_mat, t_real>(Y, m_eps);
				tl2::set_eps_0<t_mat, t_real>(Z, m_eps);
				std::cout << "x_idx=" << x_idx << ", y_idx=" << y_idx;
				std::cout << ", Q = (" << h << ", " << k << ", " << l << ")." << std::endl;
				std::cout << "V=" << V << std::endl;
				std::cout << "W=" << W << std::endl;
				std::cout << "Y=" << Y << std::endl;
				std::cout << "Z=" << Z << std::endl;*/
			} // end of coordinate iteration


			t_mat proj_norm, rot_incomm, rot_incomm_conj;
			if(m_is_incommensurate)
			{
				// equations (39) and (40) from (Toth 2015)

				proj_norm = tl2::convert<t_mat>(
					tl2::projector<t_mat_real, t_vec_real>(
						m_rotaxis, true));

				rot_incomm = tl2::unit<t_mat>(3);
				rot_incomm -= imag * tl2::skewsymmetric<t_mat, t_vec>(m_rotaxis);
				rot_incomm -= proj_norm;
				rot_incomm *= 0.5;

				rot_incomm_conj = tl2::conj(rot_incomm);
			}


			for(t_size i=0; i<num_sites*2; ++i)
			{
				//t_size site_idx = i % num_sites;
				auto& E_and_S = energies_and_correlations[i];
				const t_real& E = E_and_S.E;
				t_mat& S = E_and_S.S;
				t_mat& S_p = E_and_S.S_p;
				t_mat& S_m = E_and_S.S_m;
				t_mat& S_perp = E_and_S.S_perp;
				t_real& w = E_and_S.weight;
				t_real& w_SF1 = E_and_S.weight_spinflip[0];
				t_real& w_SF2 = E_and_S.weight_spinflip[1];
				t_real& w_NSF = E_and_S.weight_nonspinflip;

				if(m_is_incommensurate)
				{
					// TODO: formula 40 from (Toth 2015)
					S = S*proj_norm +
						S_p*rot_incomm +
						S_m*rot_incomm_conj;
				}

				if(m_temperature >= 0.)
				{
					// apply bose factor
					S *= tl2::bose_cutoff(E, m_temperature,
						m_bose_cutoff);
				}

				// apply the orthogonal projector for magnetic neutron scattering
				S_perp = (m_proj_neutron * tl2::herm(S)) * (S * m_proj_neutron);
				//S_perp = m_proj_neutron * S;

				// weights
				w_SF1 = std::abs(S_perp(0, 0).real());
				w_SF2 = std::abs(S_perp(1, 1).real());
				w_NSF = std::abs(S_perp(2, 2).real());
				//w = tl2::trace<t_mat>(S_perp).real();
				w = w_SF1 + w_SF2 + w_NSF;
			}
		}


		// unite degenerate energies and their corresponding eigenstates
		if(m_unite_degenerate_energies)
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

				// energy not yet seen
				if(iter == new_energies_and_correlations.end())
				{
					new_energies_and_correlations.push_back(curState);
				}

				// energy already seen
				else
				{
					// add correlation matrices and weights
					iter->S += curState.S;
					iter->S_p += curState.S_p;
					iter->S_m += curState.S_m;;
					iter->S_perp += curState.S_perp;
					iter->weight += curState.weight;
					iter->weight_spinflip[0] += curState.weight_spinflip[0];
					iter->weight_spinflip[1] += curState.weight_spinflip[1];
					iter->weight_nonspinflip += curState.weight_nonspinflip;
				}
			}

			energies_and_correlations = std::move(new_energies_and_correlations);
		}

		return energies_and_correlations;
	}


	/**
	 * get the energies and the spin-correlation at the given momentum
	 * @note implements the formalism given by (Toth 2015)
	 */
	std::vector<EnergyAndWeight> GetEnergies(
		t_real h, t_real k, t_real l,
		bool only_energies = false) const
	{
		t_mat H = GetHamiltonian(h, k, l);
		return GetEnergies(H, h, k, l, only_energies);
	}


	/**
	 * get the energy of the goldstone mode
	 * @note a first version for a simplified ferromagnetic dispersion was based on (Heinsdorf 2021)
	 */
	t_real GetGoldstoneEnergy() const
	{
		auto energies_and_correlations = GetEnergies(0., 0., 0., true);
		auto min_iter = std::min_element(
			energies_and_correlations.begin(), energies_and_correlations.end(),
			[](const auto& E_and_S_1, const auto& E_and_S_2) -> bool
			{
				return E_and_S_1.E < E_and_S_2.E;
			});

		if(min_iter != energies_and_correlations.end())
			return min_iter->E;

		return 0.;
	}


	/**
	 * generates the dispersion plot along the given q path
	 */
	void SaveDispersion(const std::string& filename,
		t_real h_start, t_real k_start, t_real l_start,
		t_real h_end, t_real k_end, t_real l_end,
		t_size num_qs = 128) const
	{
		std::ofstream ofstr{filename};
		ofstr.precision(m_prec);

		ofstr << std::setw(m_prec*2) << std::left << "# h"
			<< std::setw(m_prec*2) << std::left << "k"
			<< std::setw(m_prec*2) << std::left << "l"
			<< std::setw(m_prec*2) << std::left << "E"
			<< std::setw(m_prec*2) << std::left << "w"
			<< std::setw(m_prec*2) << std::left << "w_sf1"
			<< std::setw(m_prec*2) << std::left << "w_sf2"
			<< std::setw(m_prec*2) << std::left << "w_nsf"
			<< "\n";

		for(t_size i=0; i<num_qs; ++i)
		{
			t_real h = std::lerp(h_start, h_end, t_real(i)/t_real(num_qs-1));
			t_real k = std::lerp(k_start, k_end, t_real(i)/t_real(num_qs-1));
			t_real l = std::lerp(l_start, l_end, t_real(i)/t_real(num_qs-1));


			auto energies_and_correlations = GetEnergies(h, k, l, false);
			for(const auto& E_and_S : energies_and_correlations)
			{
				t_real E = E_and_S.E;
				t_real w = E_and_S.weight;
				t_real w_sf1 = E_and_S.weight_spinflip[0];
				t_real w_sf2 = E_and_S.weight_spinflip[1];
				t_real w_nsf = E_and_S.weight_nonspinflip;

				ofstr
					<< std::setw(m_prec*2) << std::left << h
					<< std::setw(m_prec*2) << std::left << k
					<< std::setw(m_prec*2) << std::left << l
					<< std::setw(m_prec*2) << E
					<< std::setw(m_prec*2) << w
					<< std::setw(m_prec*2) << w_sf1
					<< std::setw(m_prec*2) << w_sf2
					<< std::setw(m_prec*2) << w_nsf
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
		{
			return false;
		}

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

		// atom sites
		if(auto sites = node.get_child_optional("atom_sites"); sites)
		{
			for(const auto &site : *sites)
			{
				AtomSite atom_site;

				atom_site.name = site.second.get<std::string>("name", "n/a");
				atom_site.index = m_sites.size();

				atom_site.pos = tl2::create<t_vec_real>(
				{
					site.second.get<t_real>("position_x", 0.),
					site.second.get<t_real>("position_y", 0.),
					site.second.get<t_real>("position_z", 0.),
				});

				atom_site.spin_dir = tl2::create<t_vec>(
				{
					site.second.get<t_real>("spin_x", 0.),
					site.second.get<t_real>("spin_y", 0.),
					site.second.get<t_real>("spin_z", 1.),
				});

				atom_site.spin_mag = site.second.get<t_real>("spin_magnitude", 1.);
				atom_site.g = -2. * tl2::unit<t_mat>(3);

				m_sites.emplace_back(std::move(atom_site));
			}
		}

		// exchange terms
		if(auto terms = node.get_child_optional("exchange_terms"); terms)
		{
			for(const auto &term : *terms)
			{
				ExchangeTerm exchange_term;

				exchange_term.name = term.second.get<std::string>("name", "n/a");
				exchange_term.index = m_exchange_terms.size();
				exchange_term.atom1 = term.second.get<t_size>("atom_1_index", 0);
				exchange_term.atom2 = term.second.get<t_size>("atom_2_index", 0);

				exchange_term.dist = tl2::create<t_vec_real>(
				{
					term.second.get<t_real>("distance_x", 0.),
					term.second.get<t_real>("distance_y", 0.),
					term.second.get<t_real>("distance_z", 0.),
				});

				exchange_term.J = term.second.get<t_real>("interaction", 0.);

				exchange_term.dmi = tl2::create<t_vec>(
				{
					term.second.get<t_real>("dmi_x", 0.),
					term.second.get<t_real>("dmi_y", 0.),
					term.second.get<t_real>("dmi_z", 0.),
				});

				m_exchange_terms.emplace_back(std::move(exchange_term));
			}
		}

		// external field
		if(auto field = node.get_child_optional("field"); field)
		{
			m_field.mag = 0.;
			m_field.align_spins = false;

			m_field.dir = tl2::create<t_vec_real>(
			{
				field->get<t_real>("direction_h", 0.),
				field->get<t_real>("direction_k", 0.),
				field->get<t_real>("direction_l", 1.),
			});

			if(auto optVal = field->get_optional<t_real>("magnitude"))
				m_field.mag = *optVal;

			if(auto optVal = field->get_optional<bool>("align_spins"))
				m_field.align_spins = *optVal;
			}

		// temperature
		m_temperature = node.get<t_real>("temperature", -1.);

		// bragg peak
		if(auto bragg = node.get_child_optional("bragg"); bragg)
		{
			m_bragg = tl2::create<t_vec>(
			{
				bragg->get<t_real>("h", 1.),
				bragg->get<t_real>("k", 0.),
				bragg->get<t_real>("l", 0.),
			});
		}

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

		CalcSpinRotation();
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

		// bragg peak
		if(m_bragg.size() == 3)
		{
			node.put<t_real>("bragg.h", m_bragg[0].real());
			node.put<t_real>("bragg.k", m_bragg[1].real());
			node.put<t_real>("bragg.l", m_bragg[2].real());
		}

		// ordering vector
		if(m_ordering.size() == 3)
		{
			node.put<t_real>("ordering.h", m_ordering[0]);
			node.put<t_real>("ordering.k", m_ordering[1]);
			node.put<t_real>("ordering.l", m_ordering[2]);
		}

		// temperature
		node.put<t_real>("temperature", m_temperature);

		// atom sites
		for(const auto& site : GetAtomSites())
		{
			boost::property_tree::ptree itemNode;
			itemNode.put<std::string>("name", site.name);
			itemNode.put<t_real>("position_x", site.pos[0]);
			itemNode.put<t_real>("position_y", site.pos[1]);
			itemNode.put<t_real>("position_z", site.pos[2]);
			itemNode.put<t_real>("spin_x", site.spin_dir[0].real());
			itemNode.put<t_real>("spin_y", site.spin_dir[1].real());
			itemNode.put<t_real>("spin_z", site.spin_dir[2].real());
			itemNode.put<t_real>("spin_magnitude", site.spin_mag);

			node.add_child("atom_sites.site", itemNode);
		}

		// exchange terms
		for(const auto& term : GetExchangeTerms())
		{
			boost::property_tree::ptree itemNode;
			itemNode.put<std::string>("name", term.name);
			itemNode.put<t_size>("atom_1_index", term.atom1);
			itemNode.put<t_size>("atom_2_index", term.atom2);
			itemNode.put<t_real>("distance_x", term.dist[0]);
			itemNode.put<t_real>("distance_y", term.dist[1]);
			itemNode.put<t_real>("distance_z", term.dist[2]);
			itemNode.put<t_real>("interaction", term.J.real());
			itemNode.put<t_real>("dmi_x", term.dmi[0].real());
			itemNode.put<t_real>("dmi_y", term.dmi[1].real());
			itemNode.put<t_real>("dmi_z", term.dmi[2].real());

			node.add_child("exchange_terms.term", itemNode);
		}

		return true;
	}


private:
	std::vector<AtomSite> m_sites{};
	std::vector<ExchangeTerm> m_exchange_terms{};
	std::vector<AtomSiteCalc> m_sites_calc{};

	// external field
	ExternalField m_field{};
	// matrix to rotate field into the [001] direction
	t_mat m_rot_field = tl2::unit<t_mat>(3);

	// ordering wave vector for incommensurate structures
	t_vec_real m_ordering = tl2::zero<t_vec_real>(3);
	t_vec_real m_rotaxis = tl2::create<t_vec_real>({1., 0., 0.});
	bool m_is_incommensurate{false};

	bool m_unite_degenerate_energies{true};

	// bragg peak needed for calculating projector
	t_vec m_bragg{};

	// temperature (-1: disable bose factor)
	t_real m_temperature{-1};

	// bose cutoff energy to avoid infinities
	t_real m_bose_cutoff = 0.025;

	// orthogonal projector for magnetic neutron scattering
	// see (Shirane 2002), p. 37, eq. (2.64)
	t_mat m_proj_neutron = tl2::unit<t_mat>(3);

	t_size m_retries_chol = 10;
	t_real m_eps_chol = 0.05;

	t_real m_eps = 1e-6;
	int m_prec = 6;
};

}
#endif
