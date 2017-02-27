#include <MTToolBox/version.h>
#include <iostream>
#include <cstring>

int main()
{
    using namespace std;
    const char * v = get_mttoolbox_version();
    cout << v << endl;
    if (strlen(v) > 0) {
        return 0;
    } else {
        return 1;
    }
}
