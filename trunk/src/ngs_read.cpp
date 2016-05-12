#include <memory>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <stdlib.h>
#include <assert.h>
#include <algorithm>

#include <omp.h>

// include own software parts
#include "../include/ngs_global.hpp"
#include "../include/ngs_read.hpp"

ngs::Read::Read() {
    std::string tmp("");
    this->id      = std::make_shared<std::string>(std::move(tmp));
    this->seq     = std::make_shared<std::string>(std::move(tmp));
    this->comment = std::make_shared<std::string>(std::move(tmp));
    this->quality = std::make_shared<std::string>(std::move(tmp));
    valid   = true; 
    l       = 0;
    r       = 0;
}

// copy-constructor
ngs::Read::Read(const Read &r) :
                valid(r.valid),
                id(r.id),
                seq(r.seq),
                comment(r.comment),
                quality(r.quality),
                l(r.l),
                r(r.r)
                {
                }
                
// explicit constructor
ngs::Read::Read(std::string const id_,
                std::string const seq_,
                std::string const comment_,
                std::string const quality_
               ) {
    id      = std::make_shared<std::string>(id_);
    seq     = std::make_shared<std::string>(seq_);
    comment = std::make_shared<std::string>(comment_);
    quality = std::make_shared<std::string>(quality_);
    valid   = true;
    l       = 0;
    r       = seq_.size() - 1;
}

ngs::Read::~Read() {
   id.reset();
   seq.reset();
   comment.reset();
   quality.reset();
}

void ngs::Read::reverse_complement() {
    std::reverse(this->seq->begin(), this->seq->end());
    std::reverse(this->quality->begin(), this->quality->end());
    std::replace(this->seq->begin(), this->seq->end(), 'A', 'T');
    std::replace(this->seq->begin(), this->seq->end(), 'T', 'A');
    std::replace(this->seq->begin(), this->seq->end(), 'G', 'C');
    std::replace(this->seq->begin(), this->seq->end(), 'C', 'G');
}

void ngs::Read::complement() {
    std::replace(this->seq->begin(), this->seq->end(), 'A', 'T');
    std::replace(this->seq->begin(), this->seq->end(), 'T', 'A');
    std::replace(this->seq->begin(), this->seq->end(), 'G', 'C');
    std::replace(this->seq->begin(), this->seq->end(), 'C', 'G');
}

namespace ngs {
    std::ostream& operator<< (std::ostream& out, const Read& r) {
        out << *r.id      << '\n' //std::endl
            << *r.seq     << '\n' //std::endl 
            << *r.comment << '\n' //std::endl
            << *r.quality << '\n'; //std::endl;
        return out;
    }
 
}



