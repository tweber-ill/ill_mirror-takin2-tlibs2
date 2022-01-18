/**
 * tlibs2 -- helpers
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2019-2021
 * @note Forked on 2018 from my privately developed "misc" project (https://github.com/t-weber/misc).
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 * "misc" project
 * Copyright (C) 2017-2021  Tobias WEBER (privately developed).
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

#ifndef __TL2_HELPERS_H__
#define __TL2_HELPERS_H__

#include <QtCore/QLocale>
#include <locale>


namespace tl2 {


/**
 * set the "C" locale
 */
static inline void set_locales()
{
	std::ios_base::sync_with_stdio(false);

	::setlocale(LC_ALL, "C");
	std::locale::global(std::locale("C"));
	QLocale::setDefault(QLocale::C);
}


/**
 * reorder a vector according to a permutation
 */
template<class t_vec, class t_perm = std::vector<std::size_t>>
t_vec reorder(const t_vec& vec, const t_perm& perm)
{
	t_vec vec_new;
	vec_new.reserve(vec.size());

	for(std::size_t i=0; i<vec.size(); ++i)
		vec_new.push_back(vec[perm[i]]);

	return vec_new;
}


}
#endif
