/**
 * tlibs2 -- swig interface for magdyn
 * @author Tobias Weber <tweber@ill.fr>
 * @date 12-oct-2023
 * @license see 'LICENSE' file
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

%module tl2_magdyn
%{
	#include "magdyn.h"
%}


%include "std_vector.i"
%include "std_array.i"
%include "std_string.i"
%include "std_shared_ptr.i"
%include "std_complex.i"


%template(CplxD) std::complex<double>;

%template(VecD) std::vector<double>;
%template(VecCplx) std::vector<std::complex<double>>;

%template(ArrD3) std::array<double, 3>;
%template(ArrD4) std::array<double, 4>;
%template(ArrCplx3) std::array<std::complex<double>, 3>;
%template(ArrCplx4) std::array<std::complex<double>, 4>;


// ----------------------------------------------------------------------------
// matrix and vector containers
// ----------------------------------------------------------------------------
//%include "maths.h"

//%template(VectorD) tl2::vec<double>;
//%template(VectorCplx) tl2::vec<std::complex<double>>;
//%template(MatrixD) tl2::mat<double>;
//%template(MatrixCplx) tl2::mat<std::complex<double>>;
// ----------------------------------------------------------------------------


%include "magdyn.h"

// ----------------------------------------------------------------------------
// input- and output structs (and vectors of them)
// ----------------------------------------------------------------------------
%template(MagneticSiteD) tl2_mag::t_MagneticSite<
	tl2::vec<double>,
	tl2::mat<std::complex<double>>,
	double,
	std::size_t>;
%template(VecMagneticSite) std::vector<
	tl2_mag::t_MagneticSite<
		tl2::vec<double>,
		tl2::mat<std::complex<double>>,
		double,
		std::size_t>
	>;
%template(VecMagneticSitePtr) std::vector<
	const tl2_mag::t_MagneticSite<
		tl2::vec<double>,
		tl2::mat<std::complex<double>>,
		double,
		std::size_t>*
	>;

%template(MagneticSiteCalcD) tl2_mag::t_MagneticSiteCalc<
	tl2::vec<std::complex<double>>>;
%template(VecMagneticSiteCalc) std::vector<
	tl2_mag::t_MagneticSiteCalc<
		tl2::vec<std::complex<double>>>
	>;

%template(ExchangeTermD) tl2_mag::t_ExchangeTerm<
	tl2::vec<double>,
	std::size_t>;
%template(VecExchangeTerm) std::vector<
	tl2_mag::t_ExchangeTerm<
		tl2::vec<double>,
		std::size_t>
	>;

%template(ExchangeTermCalcD) tl2_mag::t_ExchangeTermCalc<
	tl2::mat<std::complex<double>>,
	tl2::vec<std::complex<double>>,
	std::complex<double>>;
%template(VecExchangeTermCalc) std::vector<
	tl2_mag::t_ExchangeTermCalc<
		tl2::mat<std::complex<double>>,
		tl2::vec<std::complex<double>>,
		std::complex<double>>
	>;

%template(ExternalFieldD) tl2_mag::t_ExternalField<
	tl2::vec<double>,
	double>;

%template(EnergyAndWeightD) tl2_mag::t_EnergyAndWeight<
	tl2::mat<std::complex<double>>,
	double>;
%template(VecEnergyAndWeight) std::vector<
	tl2_mag::t_EnergyAndWeight<
		tl2::mat<std::complex<double>>,
		double>
	>;

%template(VariableD) tl2_mag::t_Variable<
	std::complex<double>>;
%template(VecVariable) std::vector<
	tl2_mag::t_Variable<
		std::complex<double>>
	>;
// ----------------------------------------------------------------------------


/**
 * main magdyn class
 */
%template(MagDynD) tl2_mag::MagDyn<
	tl2::mat<std::complex<double>>,
	tl2::vec<std::complex<double>>,
	tl2::mat<double>,
	tl2::vec<double>,
	std::complex<double>,
	double,
	std::size_t>;
