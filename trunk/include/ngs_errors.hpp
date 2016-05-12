#include <iostream>
#include <string>
#include <stdlib.h>

#include "./ngs_global.hpp"

#ifndef NGS_ERRORS
#define NGS_ERRORS

namespace ngs {
    //! generic error function -> will issue a retrieved message and abort the program
    inline void error(const std::string &msg) {
	    std::cerr << std::endl << "ERROR: " << msg << std::endl;
	    std::exit(EXIT_FAILURE);
    }

    //! generic warning function -> will write retrieved warning message to std::cerr
    inline void warning(const std::string &msg) {
	    std::cerr << std::endl << "WARNING: " << msg << std::endl;
    }

    //! Unknown file / path name -> will issue a message and abort the program
    inline void ioerror(const std::string &fname) {
	    std::cerr << std::endl << "ERROR: Unable to open '" << fname << "'." << std::endl;
	    std::cerr << "Please check this file's name and path."  << std::endl;
	    std::exit(EXIT_FAILURE);
    }

    //! upon format sniffer errors
    inline void formating_error(const std::string &fname,
		                const std::string &reason) {
            std::cerr << std::endl << "ERROR: wrong format in '" << fname << "'." << std::endl;
	    std::cerr << "Message was: '" << reason << "'." << std::endl;
	    std::exit(EXIT_FAILURE);
    }
    
    //! asking for a feature which is not yet implemented
    inline void not_implemented_error(const std::string &feature) {
        std::cerr << std::endl << "ERROR: The feature '" << feature << "' is not yet implemented." << std::endl;
	    std::cerr << "Please inquire the runtime at the developer." << std::endl;
	    std::exit(EXIT_FAILURE);
    }
 
	inline void commandline_error(const std::string &msg) {
		std::cerr << std::endl << "ERROR parsing the command line: " << msg << std::endl;
		std::exit(EXIT_FAILURE);
	}
} // end namespace ngs

#endif // NGS_ERRORS
