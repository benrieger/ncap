#include <vector>
#include <string>
#include <iostream>

#include <boost/program_options.hpp>
using namespace boost;
namespace po = program_options;

#include "../include/ngs_confread.hpp"

namespace ngs {
    namespace conffile {
        void conf_extractor(boost::program_options::variables_map &vm,
                        std::vector<std::string> &fadaptors,
                        std::vector<std::string> &radaptors,
                        std::vector<std::string> &uadaptors) {

        ngs::conffile::convert(vm["forward"].as<std::string>(),   fadaptors);
        ngs::conffile::convert(vm["mate"].as<std::string>(),      radaptors);
        ngs::conffile::convert(vm["universal"].as<std::string>(), uadaptors);


        }
    }
    
} //~ngs
