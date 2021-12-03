/**
 * peak finder test
 * @author Tobias Weber <tweber@ill.fr>
 * @date dec-2021
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#include "libs/maths.h"
#include "libs/fit.h"


using t_real = double;


int main()
{
	const int spline_order = 3;
	const int spline_pts = 512;

	std::vector<t_real> x = tl2::linspace<t_real>(0, 360., 256);

	std::vector<t_real> y, y_neg;
	y.reserve(x.size());
	y_neg.reserve(x.size());

	for(t_real val : x)
	{
		t_real val_fkt = std::sin(val / 180. * tl2::pi<t_real>);

		y.push_back(val_fkt);
		y_neg.push_back(-val_fkt);
	}


	std::vector<t_real> maxima_x, maxima_sizes, maxima_widths;
	tl2::find_peaks(x.size(), x.data(), y.data(), spline_order,
		maxima_x, maxima_sizes, maxima_widths,
		spline_pts, 1e-8);

	std::vector<t_real> minima_x, minima_sizes, minima_widths;
	tl2::find_peaks(x.size(), x.data(), y_neg.data(), spline_order,
		minima_x, minima_sizes, minima_widths,
		spline_pts, 1e-8);

	std::cout << "local maxima:" << std::endl;
	for(std::size_t i=0; i<maxima_x.size(); ++i)
		std::cout << "\tx = " << maxima_x[i] << std::endl;

	std::cout << "local minima:" << std::endl;
	for(std::size_t i=0; i<minima_x.size(); ++i)
		std::cout << "\tx = " << minima_x[i] << std::endl;

	return 0;
}
