#include <iostream>
#include <memory>
#include <string.h>
#include <stdio.h>

#include <boost/tuple/tuple.hpp>

#include "ngs_global.hpp"
#include "ngs_interface.hpp"

#ifndef NGS_READ
#define NGS_READ

namespace ngs {
    /* \brief Class for holding NGS short reads in fastq format
       
       \class Reads
       Will hold exactely 1 read for instance
    */
    class Read {
        friend std::ostream& operator<< (std::ostream& out, const Read& read);
            // id as number for internal use
            static unsigned int next_id;
        public:
            Read();
            Read(const Read& r);
            Read(std::string const id_,       //= ""
                 std::string const seq_,      //= ""
                 std::string const comment_,  //= ""
                 std::string const quality_); //= ""
            Read(std::shared_ptr<std::string> id_,
                 std::shared_ptr<std::string> seq_,
                 std::shared_ptr<std::string> comment_,
                 std::shared_ptr<std::string> quality_) {
                 this->id      = id_;
                 this->seq     = seq_;
                 this->comment = comment_;
                 this->quality = quality_;
                 this->l       = 0;
                 this->r       = seq_->size();
                 this->valid   = true;
            }
            virtual inline const Read& operator=(const Read &r) {
                id      = r.id;
                seq     = r.seq;
                comment = r.comment;
                quality = r.quality;
                valid   = r.valid;
                return *this;
            };
            virtual ~Read();

            // defining getters
            //! retrieve id of the read
            inline std::string get_id()      const {return *id;};
            //! retrieve the sequence
            inline std::string get_seq()     const {return *seq;};
            //! retrieve length of the sequence
            inline size_t get_length()       const {return seq->size();};
            //! retrieve the comment
            inline std::string get_comment() const {return *comment;};
            //! retrieve the quality string
            inline std::string get_quality() const {return *quality;};
            //! retrieve PHRED score
            unsigned int get_phred_score();
	    //! get left limit
	    inline size_t get_l() const {return l;};
            //! get right limit
	    inline size_t get_r() const {return r;};
            //! retrieve validity information
            inline bool is_valid() const {return valid;};
            //! retrieve the current size of the object
            inline size_t size() const {return this->r - this->l;};

            //! set sequence
            inline void set_seq(std::string _seq) {this->seq = std::make_shared<std::string>(std::move(_seq));};
            //! set the quality string
            inline void set_quality(std::string _qual) {this->quality = std::make_shared<std::string>(std::move(_qual));};
            //! set validity
            inline void set_validity(bool const &val) {this->valid = val;};
            //! set validity with comment
            
            inline void set_validity(bool const &val,
                                     std::string const c) {
                this->valid = val;
                this->comment = std::make_shared<std::string>(c);
            }
            
            //! set left limit
            inline void set_left(size_t const _l) {this->l = _l;};
            //! set right limit
            inline void set_right(size_t const _r) {this->r = _r;};

            //! reverse complement
            void reverse_complement();
            //! complement
            void complement();
            //! for the reverse complement we need a nullfunc as well:
            inline void nullfunc() {};

            inline const char& operator[] (const size_t i) {return seq->at(i);};

            // public members to access the sequence string faster
            inline std::string::const_iterator s_begin() {return seq->begin();};
            inline std::string::const_iterator s_end()   {return seq->end();};
            inline std::string::const_iterator q_begin() {return quality->begin();};
            inline std::string::const_iterator q_end()   {return quality->end();};
            inline void adapt_limits() {
                if ( l > 0 || r < this->seq->size() ) {
                    if (__builtin_expect((this->l > 0),0)) {
                        this->seq->erase(seq->begin() + this->l, seq->end());
                        this->quality->erase(quality->begin() + this->l, quality->end());
                    }
                    this->seq->erase(this->seq->begin() + this->r, this->seq->end());
                    this->quality->erase(this->quality->begin() + this->r, this->quality->end());
                }
            };

	friend inline bool operator==(const Read& lhs, const Read& rhs) {
            return strncmp(lhs.seq->c_str(), rhs.seq->c_str(), lhs.seq->size()) == 0 &&
		   lhs.seq->size() == rhs.seq->size();
	};

        friend inline bool operator>=(const Read& lhs, const Read& rhs) {
            unsigned int lhss = 0;
	    unsigned int rhss = 0;
	    for (size_t i = 0; i < lhs.seq->size(); i++) {
		 lhss += lhs.seq->at(i);
		 rhss += rhs.seq->at(i);
	    }
            return lhss >= rhss;
        };
			
        private:
            //! validity indicator
            bool valid;

        protected:
            std::shared_ptr<std::string> id;
            std::shared_ptr<std::string> seq;
            std::shared_ptr<std::string> comment;
            std::shared_ptr<std::string> quality;

        public:
            // we need two limits, to be set upon each change
            size_t l, r;
    };
} //~ngs namespace

#endif //NGS_READ
