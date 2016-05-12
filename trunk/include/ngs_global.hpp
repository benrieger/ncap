#include "ngs_interface.hpp"

#ifndef NGS_GLOBALS
#define NGS_GLOBALS

#ifdef MAIN
    #define EXTERN
#else
    #define EXTERN extern
#endif

// ngs user data as kept in the UsageData struct
EXTERN ngs::UsageData ud;

#endif // NGS_GLOBALS
