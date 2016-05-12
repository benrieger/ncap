#include <ctime>
#include <string>
#include <iostream>
using namespace std;

#ifndef TIMERHPP
#define TIMERHPP

class Timer {
    public:
        Timer(const string &name) : iname(name), 
                                start(clock()) {}
        ~Timer();
    private:
        string iname;
        clock_t start;
};
 
inline Timer::~Timer() {
    cerr << "timer " << iname << " took "
         << double(clock() - start) / CLOCKS_PER_SEC
         << " seconds" << endl;
}
#endif //TIMER

/* The following part shows a brief, but complete usage
   example:

int main() {
   vector<int> vi;
   {
       Timer t0("some name");
       for (int i = 0; i < 50000; ++i) {
           vi.insert(vi.begin(), rand());
       }
   } // implicit call of ~Timer
   {
       Timer t1("other name");
       vector<int>::const_iterator imax = max_element(vi.begin(), vi.end());
       cout << "The largest element is " << *imax << endl;
       vector<int>::const_iterator imin = min_element(vi.begin(), vi.end());
       cout << "The smallest element is " << *imin << endl;
   } // implicit call of ~Timer
}
*/
