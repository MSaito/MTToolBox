int main() {
    XorShift xs(1);
    GF2X poly;
    minpoly<uint32_t>(poly, xs);
    cout << "degree = " << deg(poly) << endl;
    // 2^128 -1 の素因子
    const char * factors128_1[] = {
        "3", "5", "17", "257", "641", "65537", "274177", "6700417",
        "67280421310721", NULL};
    if (isPrime(poly, 128, factors128_1)) {
        cout << "period is 2^128 -1." << endl;
    } else {
        cout << "period is unknown." << endl;
    }
    return 0;
}
