#ifndef HELPER_CLASSES
#define HELPER_CLASSES

class NumBits {
    public:
        uint8_t higher;
        uint8_t lower;
        NumBits() {}
        NumBits(uint8_t _higher, uint8_t _lower) : higher(_higher), lower(_lower) {}
        NumBits(const NumBits & nbits) : higher(nbits.higher), lower(nbits.lower) {}
        void print() {
            std::cout << std::setw(STR_LEN) << "Number of bits in Bh" << " : " << std::setw(NUM_LEN) << (int)higher << std::endl;
            std::cout << std::setw(STR_LEN) << "Number of bits in Bl" << " : " << std::setw(NUM_LEN) << (int)lower << std::endl;
        }
};

template<typename CSROrdinal, typename Value>
class UpperBounds {
    public:
        CSROrdinal max_widthComp;
        CSROrdinal max_widthH;
        CSROrdinal max_width;
        CSROrdinal max_collisionsH;
        CSROrdinal max_degA;
        long flops;

        UpperBounds (CSROrdinal _max_widthComp, CSROrdinal _max_widthH, CSROrdinal _max_width, CSROrdinal _max_collisionsH, CSROrdinal _max_degA, long _flops) :
            max_widthComp(_max_widthComp), max_widthH(_max_widthH), max_width(_max_width), max_collisionsH(_max_collisionsH), max_degA(_max_degA), flops(_flops) {}
        UpperBounds (const UpperBounds & that) :
            max_widthComp(that.max_widthComp), max_widthH(that.max_widthH), max_width(that.max_width), max_collisionsH(that.max_collisionsH), max_degA(that.max_degA), flops(that.flops) {}

        UpperBounds & operator=(const UpperBounds & that) {
            max_widthComp = that.max_widthComp;
            max_widthH = that.max_widthH;
            max_width = that.max_width;
            max_collisionsH = that.max_collisionsH;
            max_degA = that.max_degA;
            flops = that.flops;
            return *this;
        }

        void print() {
            std::cout << std::setw(STR_LEN) << "max_widthComp" << " : " << std::setw(NUM_LEN) << max_widthComp << std::endl;
            std::cout << std::setw(STR_LEN) << "max_widthH" << " : " << std::setw(NUM_LEN) << max_widthH << std::endl;
            std::cout << std::setw(STR_LEN) << "max_width" << " : " << std::setw(NUM_LEN) << max_width << std::endl;
            std::cout << std::setw(STR_LEN) << "max_collisionsH" << " : " << std::setw(NUM_LEN) << max_collisionsH << std::endl;
            std::cout << std::setw(STR_LEN) << "max_degA" << " : " << std::setw(NUM_LEN) << max_degA << std::endl;
            std::cout << std::setw(STR_LEN) << "flops" << " : " << std::setw(NUM_LEN) << flops << std::endl;
        }
};

template<typename CSROrdinal, typename Value, typename HOrdType, typename LOrdType>
class TwoLevelMatrix {
    public:
        Matrix<HOrdType, CSROrdinal, CSROrdinal> H;
        LOrdType * L;
        Value * values;
        NumBits num_bits;

        // For Compression usage in symbolic phase
        CSROrdinal * row_ptr_compress;
        CSROrdinal compressed_edge_num;
        CSROrdinal * seg_ptr_compress;
        uint16_t * LA_compress;
        LOrdType * L_compress;
        uint32_t * values_compress;

        void setNULL() {
            L=NULL;
            values=NULL;
            row_ptr_compress=NULL;
            seg_ptr_compress=NULL;
            LA_compress=NULL;
            L_compress=NULL;
            values_compress=NULL;
        }

        TwoLevelMatrix(NumBits _num_bits) : num_bits(_num_bits), compressed_edge_num(0) { setNULL(); }

        ~TwoLevelMatrix() {
            if (L)
                free(L);
            if (row_ptr_compress)
                free(row_ptr_compress);
            if (seg_ptr_compress)
                free(seg_ptr_compress);
            if (LA_compress)
                free(LA_compress);
            if (L_compress)
                free(L_compress);
            if (values_compress)
                free(values_compress);
            setNULL();
        }

        TwoLevelMatrix(const TwoLevelMatrix &) = delete;
        TwoLevelMatrix(TwoLevelMatrix &&) = delete;
        TwoLevelMatrix& operator=(const TwoLevelMatrix &) = delete;
        TwoLevelMatrix& operator=(TwoLevelMatrix &&) = delete;
};

#endif // HELPER_CLASSES
