#ifndef __CINT__

#include "SaschaPhysics.h"
#include <time.h>

//extern template
//struct enum_traits<particle_id>;
extern constexpr particle_id enum_traits<particle_id>::enumerators[];
extern constexpr particle_id enum_traits<particle_id>::charged[];

using namespace std;

/**
 * @brief the main routine
 * @param argc number of parameters
 * @param argv the parameters as strings
 * @return exit code
 */
int main(int argc, char **argv)
{
	clock_t start, end;
	start = clock();

	// Create instance of analysis class
	SaschaPhysics* analysis = new SaschaPhysics;

	// Perform basic configuration
	if (!analysis->BaseConfig(argc, argv, "GoAT", "Physics")) {
		system("man ./documents/goat.man");
		return 1;
	}

	// Perform full initialisation
	if (!analysis->GConfigFile::Init()) {
		cout << "ERROR: Init failed!" << endl;
		return 1;
	}

	// Run over files
	analysis->TraverseFiles();

	end = clock();
	cout << endl
	<< "Time required for execution: "
	<< (double)(end-start)/CLOCKS_PER_SEC
	<< " seconds." << endl << endl;

	if (analysis)
		delete analysis;

	return 0;
}

#endif
