class XorShift : public  EquidistributionCalculatable<uint32_t> {
public:
    XorShift(uint32_t v) { seed(v); }
    XorShift(const XorShift& src) :
        EquidistributionCalculatable<uint32_t>() {
        x = src.x; y = src.y; z = src.z; w = src.w;
    }
    XorShift * clone() const { return new XorShift(*this);}
    uint32_t generate();// 省略
    uint32_t generate(int outBitLen) {
        uint32_t mask = 0;
        mask = (~mask) << (32 - outBitLen);
        return generate() & mask;
    }
    bool isZero() const { return x == 0 && y == 0 &&
            z == 0 && w == 0; }
    void setZero() { x = y = z = w = 0; }
