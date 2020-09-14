/**
 * minimal minimisation/fitting interface to minuit
 * @author Tobias Weber <tweber@ill.fr>
 * @date 14-sep-20
 * @license GPLv3, see 'LICENSE' file
 *
 * g++-10 -std=c++20 -I.. -o min min1.cpp -lMinuit2 -lgomp
 */

#include <vector>
#include <string>
#include <mutex>

#include <Minuit2/FCNBase.h>
#include <Minuit2/MnFcn.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnPrint.h>


using t_real_min = double;


class MinuitFunc : public ROOT::Minuit2::FCNBase
{
public:
	virtual t_real_min operator()(const std::vector<t_real_min>& params) const override
	{
		// in case the function cannot be called in multiple threads
		// remove the mutex otherwise
		std::lock_guard<std::mutex> _lock{m_mtx};

		std::cout << "called func with params: ";
		for(t_real_min param : params)
			std::cout << param << ", ";
		std::cout << std::endl;

		//return chi^2 here
		return 0.;
	}


	virtual t_real_min Up() const override
	{
		// sigma^2
		return 1.;
	}


private:
	mutable std::mutex m_mtx{};
};



int main()
{
	// add parameters
	std::vector<std::string> names{{std::string{"param1"}, std::string{"param2"}}};
	std::vector<t_real_min> vals{{0., 0.}};
	std::vector<t_real_min> errs{{1., 1.}};

	ROOT::Minuit2::MnUserParameters params;
	for(std::size_t param=0; param<names.size(); ++param)
	{
		params.Add(names[param], 
			static_cast<t_real_min>(vals[param]), 
			static_cast<t_real_min>(errs[param]));
	}


	// minimise
	MinuitFunc fkt;
	ROOT::Minuit2::MnMigrad migrad{fkt, params, 2};
	ROOT::Minuit2::FunctionMinimum mini = migrad();
	bool mini_valid = mini.IsValid() && mini.HasValidParameters() && mini.UserState().IsValid();


	// get back minimised parameters
	for(std::size_t param=0; param<names.size(); ++param)
	{
		vals[param] = static_cast<t_real_min>(mini.UserState().Value(names[param]));
		errs[param] = static_cast<t_real_min>(std::abs(mini.UserState().Error(names[param])));
	}

	std::cout << mini << std::endl;


	return 0;
}
