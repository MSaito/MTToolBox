    void add(EquidistributionCalculatable<uint32_t>& other) {
        XorShift* that = dynamic_cast<XorShift *>(&other);
        if (that == 0) {
            throw std::invalid_argument(
                "the adder should have the same type as the addee.");
        }
        x ^= that->x; y ^= that->y; z ^= that->z; w ^= that->w;
    }
    void seed(uint32_t v);
    int bitSize() const { return 128; }
    void setUpParam(AbstractGenerator<uint32_t>& generator){}
    const std::string getHeaderString(){return ""; }
    const std::string getParamString(){return ""; }
private:
    enum {a = 5, b = 14, c = 1};
    uint32_t x, y, z, w;
};
