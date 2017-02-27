#include "config.h"
#include <MTToolBox/version.h>

extern "C" {
    const char * get_mttoolbox_version()
    {
        return VERSION;
    }
}
