#include <string>
#include <vector>
#include <limits>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "ngs_errors.hpp"

#ifndef NGS_CONF_READ
#define NGS_CONF_READ

namespace ngs {
    namespace conffile {
    //! custom NaN-type for catching 'wrong' numbers
    const unsigned short int WRONG = std::numeric_limits<unsigned short int >::max();

    /*! will split the given string at ',', ';' and ' ' and push the
        obtained fields into a vector of type T
    */
    template<typename T>
    inline void convert(const std::string &s, std::vector<T> &v) {
        std::list<std::string> tmp;
        boost::algorithm::split(tmp, s, boost::is_any_of(",; "),
                                boost::token_compress_on);
        for (auto&& p: tmp) {
            try {
                v.push_back(boost::lexical_cast<T>(p));
            }
            catch (boost::bad_lexical_cast &msg) {
                v.push_back(boost::lexical_cast<T>(WRONG));
            }
        }
    }

    /*! will extract information from the configuration file
    */
    void conf_extractor(boost::program_options::variables_map &vm,
                        std::vector<std::string> &fadaptors,
                        std::vector<std::string> &radaptors); 
    } //~ conffile
        
} //~ ngs

#endif 
