#ifndef SPARSE_HASHMAP_H

#include <iostream>
#include "asserter.h"
#include "spgemm_defs.h"

template<typename CSROrdinal, typename Value, bool NUMERIC, typename HOrdType>
class SparseHashMap {
    private:
        CSROrdinal * densearray;
        HOrdType * keys;
        CSROrdinal * collisions;

        Value * values;
        Value * values_collisions;

        const CSROrdinal initial_value;

        CSROrdinal densearray_size;
        CSROrdinal keys_size;
        CSROrdinal collisions_size;

        CSROrdinal used_densearray_size;
        CSROrdinal used_collisions_size;

    public:
        SparseHashMap (CSROrdinal width, UpperBounds<CSROrdinal, Value> upb) : initial_value(inf<CSROrdinal>()) {

            densearray_size = width;
            keys_size = std::min(upb.max_widthH, width);
            collisions_size = upb.max_collisionsH;

            // each value has two elements included
            densearray = new CSROrdinal[densearray_size<<1];
            for ( CSROrdinal i=0; i<(densearray_size<<1); ++i)
                densearray[i] = initial_value;

            keys = new HOrdType[keys_size];
            // each value has two elements included
            collisions = new CSROrdinal[collisions_size << 1];

            used_densearray_size = 0;
            used_collisions_size = 0;

            if (NUMERIC) {
                values = new Value[keys_size];
                values_collisions = new Value[collisions_size];
            } else {
                values = NULL;
                values_collisions = NULL;
            }
        }

        ~SparseHashMap() {
            delete [] densearray;
            delete [] keys;
            delete [] collisions;

            if (NUMERIC) {
                delete [] values;
                delete [] values_collisions;
            }
        }

        SparseHashMap(const SparseHashMap&) = delete;
        SparseHashMap(SparseHashMap&&) = delete;
        SparseHashMap& operator=(const SparseHashMap&) = delete;
        SparseHashMap& operator=(SparseHashMap&&) = delete;

        void insert(HOrdType key, CSROrdinal pos, Value v) {
            if (densearray[(CSROrdinal)key<<1] == initial_value) {
                densearray[(CSROrdinal)key<<1] = pos;
                densearray[((CSROrdinal)key<<1) + 1] = initial_value;
                assert(used_densearray_size < keys_size);
                values[used_densearray_size] = v;
                keys[used_densearray_size++] = key;
            } else {
                assert(used_collisions_size < collisions_size);
                collisions[used_collisions_size<<1] = pos;
                collisions[(used_collisions_size<<1) + 1] = densearray[((CSROrdinal)key<<1) + 1];
                values_collisions[used_collisions_size] = v;
                densearray[((CSROrdinal)key<<1) + 1] = used_collisions_size;
                used_collisions_size++;
            }
        }

        void insert(HOrdType key, CSROrdinal value) {
            assert(key < densearray_size);
            assert(value != initial_value);

            if (densearray[(CSROrdinal)key<<1] == initial_value) {
                densearray[(CSROrdinal)key<<1] = value;
                densearray[((CSROrdinal)key<<1) + 1] = initial_value;
                assert(used_densearray_size < keys_size);
                keys[used_densearray_size++] = key;
            } else {
                assert(used_collisions_size < collisions_size);
                collisions[used_collisions_size<<1] = value;
                collisions[(used_collisions_size<<1) + 1] = densearray[((CSROrdinal)key<<1) + 1];
                densearray[((CSROrdinal)key<<1) + 1] = used_collisions_size;
                used_collisions_size++;
            }
        }

        void getKey(CSROrdinal id, HOrdType & key, Value& val) {
            assert(NUMERIC);
            assert(id < used_densearray_size);
            val = values[id];
            key = keys[id];
        }

        void getKey(CSROrdinal id, HOrdType & key) {
            assert(id < used_densearray_size);
            key = keys[id];
        }

        CSROrdinal getKeysNum() {
            return used_densearray_size;
        }

        bool getCollisions(CSROrdinal & pos, CSROrdinal & cur, Value & v) {
            assert(NUMERIC);

            if (cur == initial_value)
                return false;

            assert(cur < used_collisions_size);
            assert(collisions[cur << 1] != initial_value);

            v = values_collisions[cur];
            pos = collisions[cur << 1];
            cur = collisions[(cur << 1) + 1];
            return true;
        }

        bool getFirstValue(HOrdType key, CSROrdinal & value, CSROrdinal & next) {
            assert(key < densearray_size);
            assert(densearray[(CSROrdinal)key << 1] != initial_value);

            value = densearray[(CSROrdinal)key << 1];
            next = densearray[((CSROrdinal)key << 1) + 1];
            densearray[(CSROrdinal)key << 1] = initial_value;
            return next != initial_value;
        }

        bool getCollisions(CSROrdinal & value, CSROrdinal & cur) {
            if (cur == initial_value)
                return false;

            assert(cur < used_collisions_size);
            assert(collisions[cur << 1] != initial_value);

            value = collisions[cur << 1];
            cur = collisions[(cur << 1) + 1];
            return true;
        }

        void resetSize() {
            used_densearray_size = 0;
            used_collisions_size = 0;
        }

        void reset() {
            for ( CSROrdinal i=0; i<used_densearray_size; ++i) {
                densearray[ (CSROrdinal)keys[i] << 1 ] = initial_value;
                densearray[ ((CSROrdinal)keys[i] << 1) + 1 ] = initial_value;
            }
            resetSize();
        }

        void printSize() {
            std::cout << "used densearray size: " << used_densearray_size << std::endl;
            std::cout << "used collisions size: " << used_collisions_size << std::endl;
        }

        void print() {
            std::cout << "used densearray size: " << used_densearray_size << std::endl;
            std::cout << "keys are : " << std::endl;
            for (CSROrdinal i=0; i<used_densearray_size; ++i)
                std::cout << keys[i] << " ";
            std::cout << std::endl;
            std::cout << "In densearray"<< densearray_size <<" " << keys_size << " " << collisions_size << std::endl;
            for (CSROrdinal i=0; i<densearray_size; ++i) {
                if (densearray[i<<1] != initial_value) {
                    std::cout << i << std::endl;
                }
            }
            std::cout << "collisions are:" << std::endl;
            for (CSROrdinal i=0; i<used_collisions_size; ++i)
                std::cout << collisions[i<<1] << " "  << collisions[(i<<1)+1] << std::endl;
        }
};

#endif
