#ifndef NGS_STATS
#define NGS_STATS

#include <map>

namespace ngs {
    //! forward declaration
    class Read;

    //! left-over length per read
    typedef std::map<size_t, size_t> ReadLength;
    //! stats per read
    typedef std::map<size_t, char> ReadStat;
    //! summary stats for data
    typedef std::map<size_t, std::map<char, unsigned int > > dstats; // map cannot be beaten for its performance??
    //! summary stats for quality
    typedef std::map<unsigned int, unsigned int> qstats;

    /** counter class for thread-local quality stats
      *
      */
    class Counter {
        public:
            Counter();
            //Counter(size_t readlength);
            static Counter* getCounter() { return new Counter();};
            
            // half-private getter functions for GlobalStats instances
            inline const dstats _get_full()     const {return full_;};
            inline const dstats _get_trimmed()  const {return trimmed_;};

            //! a quality filter and gather function
            void qual_fil(Read &r,
                          std::vector<Read> *o_r,
                          size_t position);
            //! a quality filter and gather function - paired end scenario
            void qual_fil_paired_end(Read &r,
                                     std::vector<Read> *o_r,
                                     size_t position);
            // getter for size
            inline size_t get_size() const { return size;};

        private:
            // the size for the current read type
            size_t size;
            // summary for full data
            dstats full_;
            // summary for trimmed data
            dstats trimmed_;
            // read lengths
            ReadLength readlength_;

            //! this should never be used or linked
            Counter& operator=(const Counter &c);
    };

    /** a statistics overview class
     *
     */
    class GlobalStats final : public Counter {
        public:
            /** retrieve a pointer to a Statistics instance
               invoke like: Statistics* pointer = Statistics::get_Statistics();
            */
            static GlobalStats* get_Statistics();

            GlobalStats& operator+=(const Counter c);

        private:
            //! the singleton instance
            static GlobalStats* inst_; 
            GlobalStats();
            GlobalStats(const GlobalStats&);
            GlobalStats& operator=(const GlobalStats&);
            // summary for full data
            dstats full_;
            // summary for trimmed data
            dstats trimmed_;
    };

} //~ngs namespace

#endif //NGS_STATS
