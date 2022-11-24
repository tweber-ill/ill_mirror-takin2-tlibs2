/**
 * tlibs test file
 * @author Tobias Weber <tweber@ill.fr>
 * @date 18/oct/2022
 * @license GPLv3
 *
 * ----------------------------------------------------------------------------
 * tlibs -- a physical-mathematical C++ template library
 * Copyright (C) 2017-2022  Tobias WEBER (Institut Laue-Langevin (ILL),
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

// g++ -std=c++20 -I.. -I /usr/include/hdf5/serial -o h5file h5file.cpp -L /usr/lib/x86_64-linux-gnu/hdf5/serial/ -lhdf5_cpp -lhdf5

#include <iostream>
#include <string>
#include "libs/h5file.h"


void read_tests()
{
	H5::H5File h5file("/users/tw/Downloads/mail_tmp/065006.nxs", H5F_ACC_RDONLY);

	std::vector<std::string> entries;
	tl2::get_h5_entries(h5file, "/", entries);
	for(const std::string& entry : entries)
		std::cout << entry << " ";
	std::cout << std::endl;

	std::cout << std::endl;

	int steps = 0;
	bool ok = tl2::get_h5_scalar(h5file, "entry0/data_scan/total_steps", steps);
	std::cout << std::boolalpha << ok << std::endl;
	std::cout << "steps: " << steps << std::endl;

	std::cout << std::endl;

	std::string ident;
	ok = tl2::get_h5_string(h5file, "entry0/experiment_identifier", ident);
	std::cout << std::boolalpha << ok << std::endl;
	std::cout << "ident: " << ident << std::endl;

	std::cout << std::endl;

	std::vector<int> scanned;
	ok = tl2::get_h5_vector(h5file, "entry0/data_scan/scanned_variables/variables_names/scanned", scanned);
	std::cout << std::boolalpha << ok << std::endl;
	for(int i : scanned)
	std::cout << i << " ";
	std::cout << std::endl;

	std::cout << std::endl;

	std::vector<std::string> property;
	ok = tl2::get_h5_string_vector(h5file, "entry0/data_scan/scanned_variables/variables_names/property", property);
	std::cout << std::boolalpha << ok << std::endl;
	for(const std::string& str : property)
	std::cout << str << " ";
	std::cout << std::endl;

	std::cout << std::endl;

	std::vector<std::vector<double>> data;
	ok = tl2::get_h5_matrix(h5file, "entry0/data_scan/scanned_variables/data", data);
	std::cout << std::boolalpha << ok << std::endl;
	for(const auto& row : data)
	{
		for(double val : row)
			std::cout << val << " ";
		std::cout << std::endl;
	}

	h5file.close();
}


void write_tests()
{
	H5::H5File h5file("test.hdf", H5F_ACC_TRUNC);
	h5file.createGroup("test_group");

	tl2::set_h5_scalar(h5file, "test_group/test_scalar", 1234.567);

	std::vector<double> vec {{ 123., 456., 789.  }};
	tl2::set_h5_vector(h5file, "test_group/test_vector", vec);

	tl2::set_h5_string<std::string>(h5file, "test_group/test_string", "abc123");

	std::vector<std::string> strvec {{ "123", "abc", "xyz"  }};
	tl2::set_h5_string_vector(h5file, "test_group/test_string_vector", strvec);

	h5file.close();
}


int main()
{
	read_tests();
	std::cout << "\n--------------------------------------------------------------------------------\n" << std::endl;
	write_tests();

	return 0;
}
