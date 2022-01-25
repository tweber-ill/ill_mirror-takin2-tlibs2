/**
 * tlibs2
 * magnon dynamics
 * @author Tobias Weber <tweber@ill.fr>
 * @date january-2022
 * @license GPLv3, see 'LICENSE' file
 *
 * References:
 *   - (Toth 2015) S. Toth and B. Lake, J. Phys.: Condens. Matter 27 166002 (2015):
 *     https://doi.org/10.1088/0953-8984/27/16/166002
 *   - (Heinsdorf 2021) N. Heinsdorf, example ferromagnetic calculation, personal communications, 2021/2022.
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
	t_vec pos{};        // atom position
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
	t_size atom1{};     // atom 1 index
	t_size atom2{};     // atom 2 index
	t_vec dist{};       // distance between unit cells
	t_cplx J{};         // Heisenberg interaction
	t_vec dmi{};        // Dzyaloshinskij-Moriya interaction
};


/**
 * terms related to an external magnetic field
 */
struct ExternalField
{
	t_vec dir{};        // field direction
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


template<class t_mat> requires tl2::is_mat<t_mat>
void dbg_print(const t_mat& mat)
{
	using t_elem = typename t_mat::value_type;

	t_real eps = 1e-6;
	int prec = 3;
	std::cout.precision(prec);

	for(std::size_t row=0; row<mat.size1(); ++row)
	{
		for(std::size_t col=0; col<mat.size2(); ++col)
		{
			t_elem elem = mat(row, col);
			tl2::set_eps_0(elem, eps);
			std::cout << std::setw(prec*3) << elem;
		}

		std::cout << std::endl;
	}
}


template<class t_vec> requires tl2::is_vec<t_vec>
void dbg_print(const t_vec& vec)
{
	using t_elem = typename t_vec::value_type;

	t_real eps = 1e-6;
	int prec = 3;
	std::cout.precision(prec);

	for(std::size_t row=0; row<vec.size(); ++row)
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

	/**
	 * clear all
	 */
	void Clear()
	{
		ClearAtomSites();
		ClearExchangeTerms();
		ClearExternalField();
		ClearTemperature();
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


	void ClearTemperature()
	{
		// -1: don't use
		m_temperature = -1.;
	}


	void SetEpsilon(t_real eps)
	{
		m_eps = eps;
	}


	void SetPrecision(int prec)
	{
		m_prec = prec;
	}


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

	void SetExternalField(const ExternalField& field)
	{
		m_field = field;
	}


	void AddAtomSite(AtomSite&& site)
	{
		m_sites.emplace_back(std::forward<AtomSite&&>(site));
	}


	void AddExchangeTerm(ExchangeTerm&& term)
	{
		m_exchange_terms.emplace_back(std::forward<ExchangeTerm&&>(term));
	}


	void AddExchangeTerm(t_size atom1, t_size atom2, const t_vec& dist, const t_cplx& J)
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
			auto [field_re, field_im] =
				tl2::split_cplx<t_vec, t_vec_real>(m_field.dir);
			m_rot_field = tl2::convert<t_mat>(
				tl2::rotation<t_mat_real, t_vec_real>(
					field_re, zdir, zdir));
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


		for(const AtomSite& site : m_sites)
		{
			// rotate local spin to ferromagnetic [001] direction
			auto [spin_re, spin_im] =
				tl2::split_cplx<t_vec, t_vec_real>(site.spin_dir);
			t_mat rot = tl2::convert<t_mat>(
				tl2::rotation<t_mat_real, t_vec_real>(
					spin_re, zdir, zdir));

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
		const t_vec Q = tl2::create<t_vec>({h, k, l});

		// constants: imaginary unit and 2pi
		constexpr const t_cplx imag{0., 1.};
		constexpr const t_real twopi = t_real(2)*tl2::pi<t_real>;
		// bohr magneton in [meV/T]
		constexpr const t_real muB = tl2::muB<t_real>
			/ tl2::meV<t_real> * tl2::tesla<t_real>;

		// build the interaction matrices J(Q) and J(-Q) of
		// equations (12) and (14) from (Toth 2015)
		t_mat J_Q = tl2::zero<t_mat>(num_sites*3, num_sites*3);
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

			t_mat J_T = tl2::trans(J);

			t_cplx phase_Q = std::exp(-imag * twopi *
				tl2::inner<t_vec>(term.dist, Q));
			t_cplx phase_mQ = std::exp(-imag * twopi *
				tl2::inner<t_vec>(term.dist, -Q));

			t_real factor = 1.; //0.5;
			add_submat<t_mat>(J_Q, factor*J*phase_Q,
				term.atom1*3, term.atom2*3);
			add_submat<t_mat>(J_Q, factor*J_T*phase_mQ,
				term.atom2*3, term.atom1*3);

			add_submat<t_mat>(J_Q0, factor*J,
				term.atom1*3, term.atom2*3);
			add_submat<t_mat>(J_Q0, factor*J_T,
				term.atom2*3, term.atom1*3);
		}


		// create the hamiltonian of equation (25) and (26) from (Toth 2015)
		t_mat A = tl2::create<t_mat>(num_sites, num_sites);
		t_mat B = tl2::create<t_mat>(num_sites, num_sites);
		t_mat C = tl2::zero<t_mat>(num_sites, num_sites);

		bool use_field = !tl2::equals_0<t_real>(m_field.mag, m_eps)
			&& m_field.dir.size() == 3;

		for(t_size i=0; i<num_sites; ++i)
		{
			for(t_size j=0; j<num_sites; ++j)
			{
				t_mat J_sub_Q = submat<t_mat>(J_Q, i*3, j*3, 3, 3);

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
				B(i, j) = factor *
					tl2::inner_noconj<t_vec>(u_i, J_sub_Q * u_j);

				if(i == j)
				{
					for(t_size k=0; k<num_sites; ++k)
					{
						// TODO: check unit of S_k
						t_real S_k = m_sites[k].spin_mag;
						const t_vec& v_k = m_sites_calc[k].v;

						t_mat J_sub_Q0 = submat<t_mat>(J_Q0, i*3, k*3, 3, 3);
						C(i, j) += S_k * tl2::inner_noconj<t_vec>(v_i, J_sub_Q0 * v_k);
					}
				}

				// include external field, equation (28) from (Toth 2015)
				if(use_field && i == j)
				{
					t_vec B = m_field.dir / tl2::norm<t_vec>(m_field.dir);
					B = B * m_field.mag;

					t_vec gv = m_sites[i].g * v_i;
					t_cplx Bgv = tl2::inner_noconj<t_vec>(B, gv);

					A(i, j) -= 0.5*muB * Bgv;
				}
			}
		}


		// test matrix block
		//return A_conj - C;

		t_mat H = tl2::zero<t_mat>(num_sites*2, num_sites*2);
		set_submat(H, A - C, 0, 0);
		set_submat(H, B, 0, num_sites);
		set_submat(H, tl2::herm(B), num_sites, 0);
		set_submat(H, tl2::conj(A) - C, num_sites, num_sites);

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

		// constants: imaginary unit and 2pi
		constexpr const t_cplx imag{0., 1.};
		constexpr const t_real twopi = t_real(2)*tl2::pi<t_real>;

		// formula 40 from (Toth 2015)
		t_mat proj_norm = tl2::projector<t_mat, t_vec>(
			tl2::create<t_vec>({1., 0., 0.}), true);

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

		// if we're not interested in the spectral weights, we can ignore duplicates
		bool remove_duplicates = only_energies;
		for(const auto& eval : evals)
		{
			if(remove_duplicates &&
				std::find_if(energies_and_correlations.begin(), energies_and_correlations.end(),
					[&eval, this](const auto& E_and_S) -> bool
				{
					t_real E = E_and_S.E;
					return tl2::equals<t_real>(E, eval.real(), m_eps);
				}) != energies_and_correlations.end())
			{
				continue;
			}

			EnergyAndWeight EandS
			{
				.E = eval.real(),
			};

			energies_and_correlations.emplace_back(std::move(EandS));
		}


		// weight factors
		if(!only_energies)
		{
			// momentum
			const t_vec Q = tl2::create<t_vec>({h, k, l});

			// get the sorting of the energies
			std::vector<t_size> sorting(energies_and_correlations.size());
			std::iota(sorting.begin(), sorting.end(), 0);

			std::stable_sort(sorting.begin(), sorting.end(),
				[&energies_and_correlations](t_size idx1, t_size idx2) -> bool
				{
					return energies_and_correlations[idx1].E >=
						energies_and_correlations[idx2].E;
				});

			//for(std::size_t idx=0; idx<sorting.size(); ++idx)
			//	std::cout << idx << " -> " << sorting[idx] << std::endl;

			//energies_and_correlations = tl2::reorder(energies_and_correlations, sorting);
			evecs = tl2::reorder(evecs, sorting);
			evals = tl2::reorder(evals, sorting);

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


			// building the spin correlation functions of equation (47) from (Toth 2015)
			for(int x_idx=0; x_idx<3; ++x_idx)
			{
				for(int y_idx=0; y_idx<3; ++y_idx)
				{
					// equations (44) from (Toth 2015)
					t_mat V = tl2::create<t_mat>(num_sites, num_sites);
					t_mat W = tl2::create<t_mat>(num_sites, num_sites);
					t_mat Y = tl2::create<t_mat>(num_sites, num_sites);
					t_mat Z = tl2::create<t_mat>(num_sites, num_sites);

					for(t_size i=0; i<num_sites; ++i)
					{
						for(t_size j=0; j<num_sites; ++j)
						{
							const t_vec& pos_i = m_sites[i].pos;
							const t_vec& pos_j = m_sites[j].pos;
							t_real S_i = m_sites[i].spin_mag;
							t_real S_j = m_sites[j].spin_mag;

							const t_vec& u_i = m_sites_calc[i].u;
							const t_vec& u_j = m_sites_calc[j].u;
							const t_vec& u_conj_i = m_sites_calc[i].u_conj;
							const t_vec& u_conj_j = m_sites_calc[j].u_conj;

							t_cplx phase_pos = std::sqrt(S_i*S_j) * std::exp(
								+imag * twopi * tl2::inner<t_vec>(pos_j - pos_i, Q));
							t_cplx phase_neg = std::sqrt(S_i*S_j) * std::exp(
								-imag * twopi * tl2::inner<t_vec>(pos_j - pos_i, Q));

							// TODO: check phase factors
							Y(i, j) = phase_pos * u_i[x_idx] * u_conj_j[y_idx];
							V(i, j) = phase_pos * u_conj_i[x_idx] * u_conj_j[y_idx];
							Z(i, j) = phase_neg * u_i[x_idx] * u_j[y_idx];
							W(i, j) = phase_neg * u_conj_i[x_idx] * u_j[y_idx];
						}
					}

					// equation (47) from (Toth 2015)
					t_mat M = tl2::create<t_mat>(num_sites*2, num_sites*2);
					set_submat(M, Y, 0, 0);
					set_submat(M, V, num_sites, 0);
					set_submat(M, Z, 0, num_sites);
					set_submat(M, W, num_sites, num_sites);

					t_mat M_trafo = trafo_herm * M * trafo;

					for(t_size i=0; i<num_sites*2; ++i)
					{
						t_mat& S = energies_and_correlations[i].S;
						S(x_idx, y_idx) += M_trafo(i, i) / t_real(2*num_sites);
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
				}
			}


			for(t_size i=0; i<num_sites*2; ++i)
			{
				//t_size site_idx = i % num_sites;
				auto& E_and_S = energies_and_correlations[i];
				const t_real& E = E_and_S.E;
				t_mat& S = E_and_S.S;
				t_mat& S_perp = E_and_S.S_perp;
				t_real& w = E_and_S.weight;
				t_real& w_SF1 = E_and_S.weight_spinflip[0];
				t_real& w_SF2 = E_and_S.weight_spinflip[1];
				t_real& w_NSF = E_and_S.weight_nonspinflip;

				// formula 40 from (Toth 2015)
				// TODO: add incommensurate cases
				//S = S * proj_norm;

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
				w = tl2::trace<t_mat>(S_perp).real();
				w_SF1 = S_perp(0, 0).real();
				w_SF2 = S_perp(1, 1).real();
				w_NSF = S_perp(2, 2).real();
			}
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
			<< std::setw(m_prec*2) << std::left << "energies"
			<< "\n";

		for(t_size i=0; i<num_qs; ++i)
		{
			t_real h = std::lerp(h_start, h_end, t_real(i)/t_real(num_qs-1));
			t_real k = std::lerp(k_start, k_end, t_real(i)/t_real(num_qs-1));
			t_real l = std::lerp(l_start, l_end, t_real(i)/t_real(num_qs-1));


			auto energies_and_correlations = GetEnergies(h, k, l, true);
			for(const auto& E_and_S : energies_and_correlations)
			{
				t_real E = E_and_S.E;

				ofstr
					<< std::setw(m_prec*2) << std::left << h
					<< std::setw(m_prec*2) << std::left << k
					<< std::setw(m_prec*2) << std::left << l
					<< std::setw(m_prec*2) << E
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

				atom_site.pos = tl2::create<t_vec>(
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
				exchange_term.atom1 = term.second.get<t_size>("atom_1_index", 0);
				exchange_term.atom2 = term.second.get<t_size>("atom_2_index", 0);

				exchange_term.dist = tl2::create<t_vec>(
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
			m_field.dir = tl2::zero<t_vec>(3);
			m_field.mag = 0.;
			m_field.align_spins = false;

			m_field.dir = tl2::create<t_vec>(
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

		CalcSpinRotation();
		return true;
	}


	/**
	 * save a configuration to a property tree
	 */
	bool Save(boost::property_tree::ptree& node) const
	{
		// external field
		node.put<t_real>("field.direction_h", m_field.dir[0].real());
		node.put<t_real>("field.direction_k", m_field.dir[1].real());
		node.put<t_real>("field.direction_l", m_field.dir[2].real());
		node.put<t_real>("field.magnitude", m_field.mag);
		node.put<bool>("field.align_spins", m_field.align_spins);

		// bragg peak
		if(m_bragg.size() == 3)
		{
			node.put<t_real>("bragg.h", m_bragg[0].real());
			node.put<t_real>("bragg.k", m_bragg[1].real());
			node.put<t_real>("bragg.l", m_bragg[2].real());
		}

		// temperature
		node.put<t_real>("temperature", m_temperature);

		// atom sites
		for(const auto& site : GetAtomSites())
		{
			boost::property_tree::ptree itemNode;
			itemNode.put<std::string>("name", site.name);
			itemNode.put<t_real>("position_x", site.pos[0].real());
			itemNode.put<t_real>("position_y", site.pos[1].real());
			itemNode.put<t_real>("position_z", site.pos[2].real());
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
			itemNode.put<t_real>("distance_x", term.dist[0].real());
			itemNode.put<t_real>("distance_y", term.dist[1].real());
			itemNode.put<t_real>("distance_z", term.dist[2].real());
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
