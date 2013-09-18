class XorShift : public AbstractGenerator<uint32_t> {
public:
    XorShift(uint32_t v) { seed(v); }
    uint32_t generate() {
        uint32_t t = (x ^ (x << a));
        x = y; y = z; z = w;
        return w = (w ^ (w >> c)) ^ (t ^ (t >> b));
    }
    void seed(uint32_t v) {
        w = ~v;
        x = y = z = v;
    }
    int bitSize() const { return 128; }
private:
    enum {a = 5, b = 14, c = 1};
    uint32_t x, y, z, w;
};
