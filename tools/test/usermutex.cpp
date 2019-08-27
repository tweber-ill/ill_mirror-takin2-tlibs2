/**
 * user (not kernel) based mutex test
 * @author Tobias Weber <tweber@ill.fr>
 * @date aug-19
 * @license GPLv3, see 'LICENSE' file
 *
 * @see https://en.wikipedia.org/wiki/Peterson%27s_algorithm
 */

#include <iostream>
#include <thread>
#include <chrono>


bool want_to_enter[] = {0, 0};
int active_thread = 0;


void func(int threadno)
{
	const int other_thread = (threadno==0 ? 1 : 0);

	for(int i=0; i<1000; ++i)
	{
		want_to_enter[threadno] = 1;
		active_thread = other_thread;

		// busy wait
		while(want_to_enter[other_thread] && active_thread==other_thread)
			std::this_thread::sleep_for(std::chrono::microseconds{1000});

		// critical section
		{
			std::cout << "In thread " << threadno << ", iteration: " << i << "." << std::endl;
		}

		want_to_enter[threadno] = 0;		
	}
}


int main()
{
	std::thread th0{func, 0};
	std::thread th1{func, 1};

	th0.join();
	th1.join();

	return 0;
}
