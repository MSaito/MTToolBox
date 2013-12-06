#include <sstream>
#include <NTL/GF2XFactoring.h>
#include <NTL/ZZ.h>
#include <NTL/vector.h>
#include <MTToolBox/AlgorithmPrimitivity.hpp>
#include <MTToolBox/period.hpp>


namespace MTToolBox {
    using namespace std;
    using namespace NTL;

    AlgorithmPrimitivity::AlgorithmPrimitivity(const char * prime_factors[])
    {
        mersenne = false;
        stringstream ss;
        long length = 0;
        for (int i = 0; prime_factors[i]; i++) {
            length = i + 1;
        }
        this->primes = new Vec<ZZ>;
        this->primes->SetLength(length);
        for (int i = 0; prime_factors[i]; i++) {
            ss << prime_factors[i];
            ss << " ";
            ZZ w;
            ss >> w;
            (*primes)[i] = w;
        }
    }

    bool AlgorithmPrimitivity::operator()(int max_degree,
                                          const NTL::GF2X& poly) const
    {
        if (max_degree != deg(poly)) {
            return false;
        }
        if (!IterIrredTest(poly)) {
            return false;
        }
        if (mersenne) {
            return true;
        }
        if (isPrime(poly, max_degree, *primes)) {
            return true;
        } else {
            return false;
        }
    }

    const AlgorithmPrimitivity MersennePrimitivity;

    const char * prime_factors2_128_1[] = {
        "3", "5", "17", "257", "641", "65537", "274177", "6700417",
        "67280421310721", NULL};

    const char * prime_factors2_160_1[] = {
        "3", "5", "11", "17", "31", "41", "257", "61681", "65537", "414721",
        "4278255361", "44479210368001", NULL};

    const char * prime_factors2_192_1[] = {
        "3", "5", "7", "13", "17", "97", "193", "241", "257", "641", "673",
        "65537", "6700417", "22253377", "18446744069414584321", NULL};

    const char * prime_factors2_224_1[] = {
        "3", "5", "17", "29", "43", "113", "127", "257", "449", "2689", "5153",
        "65537", "15790321", "183076097", "54410972897", "358429848460993",
        NULL};

    const char * prime_factors2_256_1[] = {
        "3", "5", "17", "257", "641", "65537", "274177", "6700417",
        "67280421310721", "59649589127497217", "5704689200685129054721",
        NULL};

    const char * prime_factors2_288_1[] = {
        "3", "5", "7", "13", "17", "19", "37", "73", "97", "109",
        "193", "241", "257", "433", "577", "673", "1153", "6337",
        "38737", "65537", "22253377", "38941695937", "278452876033",
        "487824887233", NULL};

    const char * prime_factors2_320_1[] = {
        "3", "5", "11", "17", "31", "41", "257",
        "641", "61681", "65537", "414721", "3602561",
        "6700417", "4278255361", "44479210368001",
        "94455684953484563055991838558081", NULL};

    const char * prime_factors2_352_1[] = {
        "3", "5", "17", "23", "89", "257", "353",
        "397", "683", "2113", "65537", "229153",
        "5304641", "119782433", "2931542417", "43872038849",
        "275509565477848842604777623828011666349761", NULL};

    const char * prime_factors2_384_1[] = {
        "3", "5", "7", "13", "17", "97", "193",
        "241", "257", "641", "673", "769", "65537",
        "274177", "6700417", "22253377", "67280421310721",
        "18446744069414584321",
        "442499826945303593556473164314770689", NULL};

    const char * prime_factors2_416_1[] = {
        "3", "5", "17", "53", "157", "257", "1613",
        "2731", "8191", "65537", "858001", "928513",
        "308761441", "18558466369", "23877647873",
        "21316654212673", "715668470267111297",
        "78919881726271091143763623681", NULL};

    const char * prime_factors2_448_1[] = {
        "3", "5", "17", "29", "43", "113", "127",
        "257", "449", "641", "2689", "5153", "65537",
        "6700417", "15790321", "183076097", "54410972897",
        "358429848460993", "167773885276849215533569",
        "37414057161322375957408148834323969", NULL};

    const char * prime_factors2_480_1[] = {
        "3", "5", "7", "11", "13", "17", "31",
        "41", "61", "97", "151", "193", "241", "257",
        "331", "673", "1321", "23041", "61681", "65537",
        "414721", "22253377", "394783681", "4278255361",
        "4562284561", "46908728641", "44479210368001",
        "14768784307009061644318236958041601", NULL};

    const char * prime_factors2_512_1[] = {
        "3", "5", "17", "257", "641", "65537",
        "274177", "6700417", "67280421310721",
        "1238926361552897", "59649589127497217",
        "5704689200685129054721",
        "93461639715357977769163558199606896584051237541638188580280321",
        NULL};

    const char * prime_factors2_544_1[] = {
        "3", "5", "17", "137", "257", "953",
        "5441", "26317", "43691", "65537", "131071",
        "354689", "383521", "2368179743873", "2879347902817",
        "373200722470799764577",
        "335631827046798245410603730138717057",
        "63406006407727721042109834220642811713",
        NULL};
}
