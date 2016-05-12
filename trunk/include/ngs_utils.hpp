// $Rev: 868 $
// $LastChangedDate: 2012-01-24 16:17:38 +0100 (Di, 24 Jan 2012) $
// $Author: meesters $

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <limits>
#include <map>
#include <iostream>
#include <stdio.h>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/filesystem.hpp>

#ifndef NGS_UTILS
#define NGS_UTILS

namespace ngs {
    //! custom NaN-type for catching 'wrong' numbers
    const unsigned short int WRONG = std::numeric_limits<unsigned short int>::max();

    //! writes a vector of arbitrary type to the selected stream
    template<typename T>
    inline void write_vector(const std::vector<T> &vec, std::ostream &out, const std::string separator = " ") {
        for (size_t i = 0; i < vec.size() - 1; i++) {
            out << vec[i] << separator;
        }
        // output last element w/o separator
        out << vec[vec.size() - 1];
    }

    //! tests for infinity
    template<typename T>
    inline bool isinf(const T value) {
        return std::numeric_limits<T>::has_infinity &&
               value == std::numeric_limits<T>::infinity();
    }

    /*! will split the given input string at ',', ';', or ' '
        and push the obtained fields into a vector of type T
    */
    template<typename T>
    inline void convert(const std::string &s, std::vector<T> &v) {
        std::list<std::string> tmp;
        //bool flag = false;
        boost::algorithm::split(tmp, s, boost::is_any_of(",; "),
                                boost::token_compress_on);
        for (std::list<std::string>::const_iterator p = tmp.begin();
            p != tmp.end(); ++p) {
                //flag = false;
                try {
                    v.push_back(boost::lexical_cast<T>(*p));
                }
                catch (boost::bad_lexical_cast &msg) {
                    //flag = true;
                    v.push_back(boost::lexical_cast<T>(WRONG));
                }
                /*if (flag) {
                    v.clear();
                    for (size_t i = 0; i < tmp.size(); ++i) {
                        v.push_back(boost::lexical_cast<T>(WRONG));
                    }
                }*/
           }
    }

    /*! will remove empty items from a string-vector
    */
    inline std::vector<std::string> compress(const std::vector<std::string> &v) {
       std::vector<std::string> tmp;
       for (std::vector<std::string>::const_iterator it = v.begin();
            it != v.end(); it++) {
           if (!it->empty())
               tmp.push_back(*it);
       }
       return tmp;
    }

    /*! check whether or not a strings ends on another string
    */
    inline bool endswith(const std::string &in, const std::string &ending) {
        unsigned int lastmatchpos = in.rfind(ending); // Find the last occurrence of ending
        bool isending = lastmatchpos != std::string::npos;
        // If the string was found, make sure that any characters that follow 
        // it are the ones we're trying to ignore
        for(size_t i = lastmatchpos + ending.length(); (i < in.length()) && isending; i++)
            if( (in[i] != '\n') && (in[i] != '\r') ) {
               isending = false;
               break;
            }
        return isending;
    }

    /*! will reset all elements of a vector to zero
    */
    template<typename T>
    inline void reset(std::vector<T> &conv, const unsigned int &size) {
        conv.clear();
        for (size_t i = 0; i < size; i++) {
            conv.push_back((T)0);
        }
    }

    /*! converts the recieved item to a string
    */
    template<typename T>
    inline std::string to_string(const T &t) {
        std::stringstream ss;
        ss << t;
        return ss.str();
    }

    //! returns the complement allele as a character
    inline char complement(const char &in) {
       char rvalue = 'X';
       if (in == 'A') rvalue = 'T';
       else if (in == 'T') rvalue = 'A';
       else if (in == 'C') rvalue = 'G';
       else if (in == 'G') rvalue = 'C';
       return rvalue;
    }
    //! returns a default value, if a map-key is not present
    template<typename K, typename V>
    inline V getwithdefault(const std::map<K, V> &m, const K &key, const V &defaultvalue) {
       typename std::map<K,V>::const_iterator it = m.find(key);
       if (it == m.end()) {
           return defaultvalue;
       }
       else {
           return it->second;
       }
    }

    //! returns a filename prefix and a suffix
    inline boost::tuple<std::string, std::string> split_filename(const std::string& filename) {
       boost::tuple<std::string, std::string> rvalue("", "");
       size_t lastdot = filename.find_last_of(".");
       // not a vali1d filename? put everything in the prefix!
       if (lastdot == std::string::npos) rvalue.get<0>() = filename;
       else {
           rvalue.get<0>() = filename.substr(0, lastdot);
           rvalue.get<1>() = filename.substr(lastdot + 1, std::string::npos);
       }
       return rvalue; 
    }

	//! will carry out a command and return stdout
	inline std::string exec(const char* cmd) {
        FILE* pipe = popen(cmd, "r");
        if (!pipe) return "ERROR";
        char buffer[128];
        std::string result = "";
        while (!feof(pipe)) {
            if (fgets(buffer, 128, pipe) != NULL)
                result += buffer;
       }
       pclose(pipe);
       return result;
    }
} // end namespace ngs
    
#endif // NGS_UTILS
