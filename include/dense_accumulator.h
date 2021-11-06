#ifndef DENSE_ACCUMULATOR_H

template<typename CSROrdinal, typename Value, typename DHMValue, typename LOrdType>
class DenseHashMap {
    public:
        DHMValue * densearray;
        LOrdType * keys;

        CSROrdinal densearray_size;
        CSROrdinal keys_size;

        CSROrdinal used_densearray_size;

        CSROrdinal count;
        DHMValue initial_value;

        DenseHashMap(CSROrdinal width, CSROrdinal max_insertion, DHMValue _initial_value=0)
            : densearray(NULL), keys(NULL), densearray_size(width), initial_value(_initial_value) {
                keys_size = std::min(width, max_insertion);

                densearray = new DHMValue[densearray_size];
                keys = new LOrdType[keys_size];

                for ( CSROrdinal i=0; i<densearray_size; ++i)
                    densearray[i] = initial_value;

                used_densearray_size = 0;
                count = 0;
            }

        ~DenseHashMap() {
            delete [] densearray;
            delete [] keys;
        }

        DenseHashMap(const DenseHashMap&) = delete;
        DenseHashMap(DenseHashMap&&) = delete;
        DenseHashMap& operator=(const DenseHashMap&) = delete;
        DenseHashMap& operator=(DenseHashMap&&) = delete;

        void insert(LOrdType key) {
            if (densearray[key] == initial_value) {
                keys[used_densearray_size++] = key;
                densearray[key] = 1;
            }
        }

        void insertOr(LOrdType key, DHMValue val) {
            assert(val != initial_value);
            assert(key < densearray_size);
            if (densearray[key] == initial_value) {
                assert(used_densearray_size < keys_size);
                keys[used_densearray_size++] = key;
                densearray[key] = val;
            } else
                densearray[key] |= val;
        }

        void insertInc(LOrdType key, Value val) {
            assert(val != initial_value);
            if (densearray[key] == initial_value) {
                keys[used_densearray_size++] = key;
                densearray[key] = val;
            } else
                densearray[key] += val;
        }

        CSROrdinal getUsedSize() {
            return used_densearray_size;
        }

        template<typename DataType>
            uint8_t bitcounts(DataType x) {
                uint8_t b = sizeof(DataType)*8;
                uint8_t c = 0;
                for(uint8_t i=0; i<b && x!=0; i++) {
                    c += (x & 1);
                    x = x>>1;
                }
                return c;
            }

        CSROrdinal getUsedSizeOr(bool reset=true) {
            CSROrdinal size = 0;
            for ( CSROrdinal i=0; i<used_densearray_size; ++i ) {
                size += bitcounts(densearray[ keys[i] ]);
                if (reset)
                    densearray[ keys[i] ] = initial_value;
            }
            return size;
        }

        bool getKeyValue(LOrdType & key, Value & val) {
            if (count >= used_densearray_size) return false;
            key = keys[count];
            val = densearray[key];
            densearray[key] = initial_value;
            count++;
            return true;
        }

        void resetSize() {
            used_densearray_size = 0;
            count = 0;
        }

        void reset() {
            for ( CSROrdinal i=0; i<used_densearray_size; ++i)
                densearray[ keys[i] ] = initial_value;
            resetSize();
        }
};

#endif
