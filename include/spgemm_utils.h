#ifndef UTILS_H_
#define UTILS_H_

#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <iostream>
#include <chrono>
#include <map>
#include <sys/stat.h>
#include <omp.h>

#include "asserter.h"
#include "hooks.h"
#include "csr_matrix.h"

#include "spgemm_defs.h"
#include "helper_classes.h"
#include "sparse_accumulator.h"
#include "dense_accumulator.h"

unsigned char MSB(uint32_t v) {
    // start counting from 1 for MSB, which means, for 4, will return 2
    if (v == 0)
        return 0;

    unsigned char msb = 0, t;
    t = ((v >> 16) == 0) ? 0 : 16;
    v = (v >> msb);
    msb += t;

    t = (((v >> 8) == 0) ? 0 : 8);
    v = (v >> t);
    msb += t;

    t = (((v >> 4) == 0) ? 0 : 4);
    v = (v >> t);
    msb += t;

    t = (((v >> 2) == 0) ? 0 : 2);
    v = (v >> t);
    msb += t;

    t = (((v >> 1) == 0) ? 0 : 1);
    v = (v >> t);
    msb += t;

    return msb+v;
}

template<typename CSROrdinal, typename Value>
NumBits getNumBits(CSROrdinal n) {
    // The inital value for lower and higher assumes :
    //    ((long)1<<16) < B.ncols() <= ((long)1<<32)
    uint8_t num_bits_lower = 16;
    uint8_t num_bits_higher = 16;
    if (n <= ((long)1<<16)) {
        num_bits_lower = 8;
        num_bits_higher = 8;
    }
    if (n > ((long)1<<32)) {
        num_bits_lower = 32;
        num_bits_higher = 32;
    }

    // based on L2 cache size, the maximum width we can have for hashmap
    CSROrdinal width = L2_CACHE_SIZE/(sizeof(Value));
    assert((width & (width - 1)) == 0); // check that width is power of 2

    // IF L2 cache can store less than the given range
    if ( width < ((long)1<<(num_bits_lower)) ) {
        if (n <= (((long)1<<num_bits_higher)*width)) // Use the width if possible
            num_bits_lower = std::min((uint8_t)(MSB(width)-1), (uint8_t)MSB(n));
    }

    CSROrdinal higher_range = (n >> num_bits_lower);
    num_bits_higher = MSB(higher_range) - ((n & (n - 1)) == 0);
    // twisting the bits to use small sized data type
    uint8_t total_bits = num_bits_higher + num_bits_lower;
    uint8_t sizes[] = {16, 16+8};
    for (auto s : sizes) {
        if (total_bits <= s) {
            if (num_bits_higher > 8 || num_bits_lower > total_bits - 8) {
                num_bits_higher = 8;
                num_bits_lower = total_bits - 8;
                break;
            }
        }
    }

    // check the correctness
    assert( (long(1) << (num_bits_higher + num_bits_lower)) >= n );
    assert( (long(1) << (num_bits_higher + num_bits_lower - 1)) < n );

    return NumBits(num_bits_higher, num_bits_lower);
}

template<typename CSROrdinal, typename Value, typename HOrdType, typename LOrdType>
void constructHL (
        Matrix<CSROrdinal, CSROrdinal, Value> & B_org,
        TwoLevelMatrix<CSROrdinal, Value, HOrdType, LOrdType> & B,
        size_t chunk_size, bool statistic_verb = false) {

    CSROrdinal nrows = B_org.nrows();
    CSROrdinal ncols = B_org.ncols();
    CSROrdinal nnz = B_org.nnz();
    CSROrdinal * ptr = B_org.Ptr();
    CSROrdinal * colindices = B_org.ColIndices();
    // In this implementation, we assume B_org is sorted within each row,
    // so values pointers are the same.
    B.values = B_org.Values();

    /* Allocate ptr */
    B.H.AllocatePtr(nrows+1);
    CSROrdinal * ptrH = B.H.Ptr();
    ptrH[0] = 0;

    // Used for checking whether to use compression or not
    B.row_ptr_compress = (CSROrdinal *) malloc (sizeof(CSROrdinal) * (nrows+1));
    B.row_ptr_compress[0] = 0;

    /* Calculate the row sizes of B.H
      and generate matrixL within the same round of visiting matrix */
    CSROrdinal lower_bits = (1 << B.num_bits.lower) - 1;
#pragma omp parallel for schedule(dynamic, chunk_size)
    for (CSROrdinal v=0; v<nrows; ++v) {
        CSROrdinal count = 0;

        // For compression
        CSROrdinal count_compress_edge = 0;

        CSROrdinal st=ptr[v];
        CSROrdinal end=ptr[v+1];

        HOrdType edgeH;
        if (st != end) {
            edgeH = ((colindices[st]) >> B.num_bits.lower);
            count++;

            //For compression
            CSROrdinal edge_comp = colindices[st] >> 5;
            count_compress_edge++;

            // generate row size with count
            for (CSROrdinal u_pos=st+1; u_pos<end; ++u_pos) {
                // neighbor list is supporsed to be sorted and no multi-colindices
                assert(colindices[u_pos] > colindices[u_pos-1]);

                count += (edgeH != (colindices[u_pos] >> B.num_bits.lower));
                edgeH = (colindices[u_pos] >> B.num_bits.lower);

                //For compression
                count_compress_edge += (edge_comp != (colindices[u_pos] >> 5));
                edge_comp = colindices[u_pos] >> 5;
            }
        }
        ptrH[v+1] = count;
        //For compression
        B.row_ptr_compress[v+1] = count_compress_edge;
    }
    std::partial_sum(ptrH, ptrH + nrows + 1, ptrH);
    //For compression
    std::partial_sum(B.row_ptr_compress, B.row_ptr_compress + nrows + 1, B.row_ptr_compress);

    float compression_ratio = (float)(B.row_ptr_compress[nrows])/nnz;
    if (compression_ratio < comp_ratio_maxium)
        USE_COMPRESSION = true;
    else
        USE_COMPRESSION = false;

    bool useHL = !(((USE_COMPRESSION && ((ncols >> 5) < (L2_CACHE_SIZE/sizeof(uint32_t)))) || (ncols < L2_CACHE_SIZE/sizeof(CSROrdinal))) && (ncols < L2_CACHE_SIZE/sizeof(Value)));

    /* Allocate colindices and values */
    HOrdType * colindicesH;
    CSROrdinal * valuesH;
    if (useHL) {
        B.H.AllocateColIndices(ptrH[nrows]);
        B.H.AllocateValues(ptrH[nrows]+1);
        colindicesH = B.H.ColIndices();
        valuesH = B.H.Values();
        B.L = (LOrdType *) malloc (sizeof(LOrdType) * nnz);
    }

    if (USE_COMPRESSION) {
        B.seg_ptr_compress = (CSROrdinal *) malloc (sizeof(CSROrdinal) * (ptrH[nrows]+1));
        B.LA_compress = (uint16_t *) malloc (sizeof(uint16_t) * (B.row_ptr_compress[nrows]));
        B.L_compress = (LOrdType *) malloc (sizeof(LOrdType) * (B.row_ptr_compress[nrows]));
        B.values_compress = (uint32_t *) malloc (sizeof(uint32_t) * (B.row_ptr_compress[nrows]));
    }

    /* get all colindices, get the data of B.H: colindices, and values */
#pragma omp parallel for schedule(dynamic, chunk_size)
    for (CSROrdinal v=0; v<nrows; ++v) {

        CSROrdinal st=ptr[v];
        CSROrdinal end=ptr[v+1];
        CSROrdinal stH=ptrH[v];

        //For compression
        if (useHL || USE_COMPRESSION) {
            CSROrdinal stCompress, edge_comp;
            uint32_t compress;
            if (USE_COMPRESSION) stCompress = B.row_ptr_compress[v];

            if (useHL) {
                // generate matrixL
                for (CSROrdinal u_pos=st; u_pos<end; ++u_pos)
                    B.L[u_pos] = (colindices[u_pos] & lower_bits);
            }

            if (st != end) {

                if(USE_COMPRESSION) {
                    edge_comp = colindices[st] >> 5;
                    B.seg_ptr_compress[stH] = stCompress;
                    compress = (1 << (colindices[st] & 31));
                }

                HOrdType edgeH=((colindices[st]) >> B.num_bits.lower);
                if (useHL) {
                    colindicesH[stH] = edgeH;
                    valuesH[stH] = st;
                    stH++;
                }

                for (CSROrdinal u_pos=st+1; u_pos<end; ++u_pos) {
                    // neighbor list is supporsed to be sorted and no multi-colindices
                    assert(colindices[u_pos] > colindices[u_pos-1]);

                    if (USE_COMPRESSION) {
                        if (edge_comp != (colindices[u_pos] >> 5)) {
                            B.values_compress[stCompress] = compress;
                            B.LA_compress[stCompress] = edge_comp;
                            B.L_compress[stCompress] = edge_comp & (lower_bits >> 5);

                            compress=(1 << (colindices[u_pos] & 31));
                            edge_comp=colindices[u_pos] >> 5;

                            stCompress++;
                        } else {
                            compress |= (1<< (colindices[u_pos] & 31));
                        }
                    }

                    if (edgeH != (colindices[u_pos] >> B.num_bits.lower)) {

                        edgeH = (colindices[u_pos] >> B.num_bits.lower);

                        if (useHL) {
                            colindicesH[stH] = edgeH;
                            valuesH[stH] = u_pos;
                            if(USE_COMPRESSION)
                                B.seg_ptr_compress[stH] = stCompress;
                            stH++;
                        }
                    }
                }

                if(USE_COMPRESSION) {
                    B.values_compress[stCompress] = compress;
                    B.LA_compress[stCompress] = edge_comp;
                    B.L_compress[stCompress] = edge_comp & (lower_bits >> 5);
                    stCompress++;
                }
            }
            if(USE_COMPRESSION)
                if (B.row_ptr_compress[v+1] != stCompress) throw -1;
        }
    }
    if (useHL)
        valuesH[ptrH[nrows]] = nnz;

    if(USE_COMPRESSION) {
        B.compressed_edge_num = B.row_ptr_compress[nrows];
        B.seg_ptr_compress[ptrH[nrows]] = B.row_ptr_compress[nrows];
    }

    if (statistic_verb) {
        std::cout << std::setw(STR_LEN) << "Compression Ratio is" << " : " << std::fixed << std::setprecision(6) << std::setw(NUM_LEN) << compression_ratio << ", and compression is";
        if (!USE_COMPRESSION)
            std::cout << " Not";
        std::cout << " applied." << std::endl;
        std::cout << std::setw(STR_LEN) << "S has non-zeros" << " : " << std::setw(NUM_LEN) << B.H.nnz();
        if (USE_COMPRESSION && ((B_org.ncols()>> 5) < (L2_CACHE_SIZE/sizeof(uint32_t))))
            std::cout << ", direct merging on dense array will be applied" << std::endl;
        else
            std::cout << std::endl;
    }
}

template<typename CSROrdinal, typename Value, typename HOrdType, typename LOrdType>
UpperBounds<CSROrdinal, Value> getUpperBounds (
        Matrix<CSROrdinal, CSROrdinal, Value> & A,
        Matrix<CSROrdinal, CSROrdinal, Value> & B_org,
        TwoLevelMatrix<CSROrdinal, Value, HOrdType, LOrdType> & B,
        size_t chunk_size, bool statistic_verb = false) {

    CSROrdinal max_widthComp=0, max_width=0, max_widthH=0, max_degA=0, max_collisions=0;
    long flops=0;

    CSROrdinal * ptrA = A.Ptr();
    CSROrdinal * ptrBH = B.H.Ptr();
    CSROrdinal * ptrBComp = B.row_ptr_compress;

#pragma omp parallel for schedule(dynamic, chunk_size)  reduction(max:max_widthComp) reduction(max:max_width)  reduction(max:max_widthH) reduction(max:max_degA) reduction(max:max_collisions) reduction(+:flops)
    for (CSROrdinal v=0; v<A.nrows(); ++v) {

        CSROrdinal st=ptrA[v];
        CSROrdinal end=ptrA[v+1];

        max_degA = (max_degA > end-st) ? max_degA : end-st;

        CSROrdinal _widthComp = 0;
        CSROrdinal _widthH = 0;
        CSROrdinal _width = 0;
        CSROrdinal _max_deg = 0;
        for (CSROrdinal u_pos=st; u_pos<end; ++u_pos) {
            CSROrdinal u = A.ColIndices()[u_pos];
            CSROrdinal u_degH = ptrBH[u+1] - ptrBH[u];
            _max_deg = (_max_deg > u_degH) ? _max_deg : u_degH;
            _widthComp += (ptrBComp[u+1] - ptrBComp[u]);
            _widthH += u_degH;
            _width += B_org.Ptr()[u+1] - B_org.Ptr()[u];
        }
        flops += _width;
        CSROrdinal _collisionH = _widthH - _max_deg;
        max_widthComp = (max_widthComp > _widthComp) ? max_widthComp : _widthComp;
        max_widthH = (max_widthH > _widthH) ? max_widthH : _widthH;
        max_width = (max_width > _width) ? max_width : _width;
        max_collisions = (max_collisions > _collisionH) ? max_collisions : _collisionH;
    }

    return UpperBounds<CSROrdinal, Value>(max_widthComp, max_widthH, max_width, max_collisions, max_degA, flops);
}
#endif
