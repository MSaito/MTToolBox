int main() {
    XorShift xs(1);
    AlgorithmEquidsitribution<uint32_t> eq(xs, 32);
    int veq[32];
    int delta = eq.get_all_equidist(veq);
    int bitSize = xs.bitSize();
    for (int i = 0; i < 32; i++) {
        cout << "k("<< dec << setw(2) << (i + 1) << "):"
             << setw(3) << veq[i] << "  d(" << setw(2) << (i + 1) << "):"
             << setw(3) << (bitSize / (i + 1)) - veq[i] << endl;
    }
    cout << "delta:" << delta << endl;
    return 0;
}
