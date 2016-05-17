#include <vector>
#include <string>

#include "ngs_read.hpp"

#ifndef NGS_ALIGNER
#define NGS_ALINGER

namespace ngs {
    // forward declaration:
    class Read;

    //! a vector of strings to hold adaptors
    typedef std::vector<std::string> v_adaptors;
    //! a map of strings, key is the starting position
    typedef std::map<size_t, std::string> m_adaptors;

    //! a null function to point to, if a certain step shall be ommitted
    inline void null_func(Read &r) {}
    //! a null function to point to, if a certain step shall be ommitted - accepts 2nd arg
    template<typename T>
    inline void null_func_t(Read &r, const T &t) {}

    //not in use: void trim_adaptor(Read &r, const std::string &adaptor);
    void trim_adaptors(Read &r, const v_adaptors &adaptors);
    void trim_primer(Read &r);
    void trim_primers(Read &r);
    void quality_filter(Read &r);

    void NWScore(char* X, char* Y) {
        size_t sY = sizeof(Y);
        size_t sX = sizeof(X);
        // int **mat = (int **)malloc(rows * sizeof(int*));
        // for(int i = 0; i < rows; i++) mat[i] = (int *)malloc(cols * sizeof(int));
        int **score = (int **)malloc(sX * sizeof(int*));
        for (size_t i = 0; i < sizeof(Y); i++) score[i] = (int *)malloc(sY * sizeof(int));

        score[0][0] = 0;
        for (size_t j = 1; j < sY; ++j) score[0][j] = score[0][j-1] - 2;



        // free the memory again
        for(size_t i = 0; i < sizeof(Y); i++) free(score[i]);
        free(score);
    }
} //~ngs

#endif
