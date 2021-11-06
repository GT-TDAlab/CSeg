#ifndef SPGEMM_HPP
#define SPGEMM_HPP

#include "spgemm_utils.h"

template<typename CSROrdinal, typename Value, typename HOrdType, typename LOrdType>
void twoLSymbolic(Matrix<CSROrdinal, CSROrdinal, Value> & A,
        Matrix<CSROrdinal, CSROrdinal, Value> & B_org,
        TwoLevelMatrix<CSROrdinal, Value, HOrdType, LOrdType> & B,
        Matrix<CSROrdinal, CSROrdinal, Value> & C,
        UpperBounds<CSROrdinal, Value> & upb,
        CSROrdinal chunk_size) {

    CSROrdinal nrowsA = A.nrows();
    CSROrdinal * ptrA = A.Ptr();
    CSROrdinal * colindicesA = A.ColIndices();

    CSROrdinal nrowsB = B_org.nrows();
    CSROrdinal ncolsB = B_org.ncols();
    CSROrdinal * ptrB = B_org.Ptr();
    CSROrdinal * colindicesB = B_org.ColIndices();

    CSROrdinal widthH = ((long)1<<B.num_bits.higher);
    CSROrdinal widthL = ((long)1<<B.num_bits.lower);

    CSROrdinal * ptrBH = B.H.Ptr();
    HOrdType * colindicesBH = B.H.ColIndices();
    CSROrdinal * valuesBH = B.H.Values();
    LOrdType * colindicesBL = B.L;

    CSROrdinal * ptrC = C.Ptr();
    ptrC[0] = 0;
    if (USE_COMPRESSION && ((nrowsB >> 5) < (L2_CACHE_SIZE/sizeof(uint32_t)))) {
#pragma omp parallel
        {
            DenseHashMap<CSROrdinal, Value, uint32_t, uint16_t> hashmapLC(((uint32_t)(ncolsB + 31))>>5, upb.max_widthComp);
#pragma omp for schedule(dynamic, chunk_size)
            for (CSROrdinal v=0; v<nrowsA; ++v) {
                // if degree of v in A is less than one, then we can just add the degree
                if (ptrA[v+1] - ptrA[v] == 1) {
                    CSROrdinal u=colindicesA[ptrA[v]];
                    ptrC[v+1] = ptrB[u+1] - ptrB[u];
                    continue;
                } else if (ptrA[v+1] - ptrA[v] == 0) {
                    ptrC[v+1] = 0;
                    continue;
                }

                for (CSROrdinal u_pos=ptrA[v]; u_pos<ptrA[v+1]; ++u_pos) {
                    CSROrdinal u=colindicesA[u_pos];
                    CSROrdinal start=B.row_ptr_compress[u];
                    CSROrdinal end=B.row_ptr_compress[u+1];
                    for (CSROrdinal k=start; k<end; ++k)
                        hashmapLC.insertOr(B.LA_compress[k], B.values_compress[k]);
                }
                ptrC[v+1] = hashmapLC.getUsedSizeOr();
                hashmapLC.resetSize();
            }
        }
    } else if (ncolsB < L2_CACHE_SIZE/sizeof(CSROrdinal)) {
#pragma omp parallel
        {
            DenseHashMap<CSROrdinal, Value, bool, CSROrdinal> hashmap(ncolsB, upb.max_width, false);
#pragma omp for schedule(dynamic, chunk_size)
            for (CSROrdinal v=0; v<nrowsA; ++v) {
                // if degree of v in A is less than one, then we can just add the degree
                if (ptrA[v+1] - ptrA[v] == 1) {
                    CSROrdinal u=colindicesA[ptrA[v]];
                    ptrC[v+1] = ptrB[u+1] - ptrB[u];
                    continue;
                } else if (ptrA[v+1] - ptrA[v] == 0) {
                    ptrC[v+1] = 0;
                    continue;
                }

                for (CSROrdinal u_pos=ptrA[v]; u_pos<ptrA[v+1]; ++u_pos) {
                    CSROrdinal u=colindicesA[u_pos];
                    CSROrdinal start=ptrB[u];
                    CSROrdinal end=ptrB[u+1];
                    for (CSROrdinal w_pos=start; w_pos<end; ++w_pos) {
                        CSROrdinal w=colindicesB[w_pos];
                        hashmap.insert(w);
                    }
                }
                ptrC[v+1] = hashmap.getUsedSizeOr();
                hashmap.reset();
            }
        }
    } else {
#pragma omp parallel
        {
            SparseHashMap<CSROrdinal, Value, false, HOrdType> hashmapH(widthH, upb);

            //For compression
            assert(widthL>>5);
            DenseHashMap<CSROrdinal, Value, uint32_t, LOrdType> hashmapLC(widthL>>5, widthL>>5);

            DenseHashMap<CSROrdinal, Value, bool, LOrdType> hashmapL(widthL, widthL);

#pragma omp for schedule(dynamic, chunk_size)
            for (CSROrdinal v=0; v<nrowsA; ++v) {

                // if degree of v in A is less than one, then we can just add the degree
                if (ptrA[v+1] - ptrA[v] == 1) {
                    CSROrdinal u=colindicesA[ptrA[v]];
                    ptrC[v+1] = valuesBH[ptrBH[u+1]] - valuesBH[ptrBH[u]];

                    continue;
                } else if (ptrA[v+1] - ptrA[v] == 0) {

                    ptrC[v+1] = 0;
                    continue;
                }

                CSROrdinal v_sizeC = 0;
                // gather all the position info
                for (CSROrdinal u_pos=ptrA[v]; u_pos<ptrA[v+1]; ++u_pos) {
                    CSROrdinal u=colindicesA[u_pos];
                    for (CSROrdinal wH_pos=ptrBH[u]; wH_pos<ptrBH[u+1]; ++wH_pos) {
                        HOrdType wH = colindicesBH[wH_pos];
                        hashmapH.insert(wH, wH_pos);
                    }
                }

                CSROrdinal keynum = hashmapH.getKeysNum();

                for (CSROrdinal k=0; k<keynum; ++k) {
                    HOrdType colindex;
                    hashmapH.getKey(k, colindex);

                    CSROrdinal next;
                    CSROrdinal pos;
                    bool moreValue = hashmapH.getFirstValue(colindex, pos, next);
                    if (!moreValue) {
                        v_sizeC += valuesBH[pos+1] - valuesBH[pos];
                    } else {
                        if(USE_COMPRESSION) {
                            while(true) {
                                CSROrdinal start = B.seg_ptr_compress[pos];
                                CSROrdinal end = B.seg_ptr_compress[pos+1];
                                for (CSROrdinal k=start; k<end; ++k)
                                    hashmapLC.insertOr(B.L_compress[k], B.values_compress[k]);
                                if (!hashmapH.getCollisions(pos, next))
                                    break;
                            }

                            v_sizeC += hashmapLC.getUsedSizeOr();
                            hashmapLC.resetSize();
                        } else {
                            while(true) {
                                CSROrdinal start = valuesBH[pos];
                                CSROrdinal end = valuesBH[pos+1];
                                for (CSROrdinal k=start; k<end; ++k)
                                    hashmapL.insert(colindicesBL[k]);
                                if (!hashmapH.getCollisions(pos, next))
                                    break;
                            }

                            v_sizeC += hashmapL.getUsedSize();
                            hashmapL.reset();
                        }
                    }
                }
                hashmapH.resetSize();
                // store the size of row in C
                ptrC[v+1] = v_sizeC;
            }
        }

    }
    std::partial_sum(ptrC, ptrC + nrowsA + 1, ptrC);
}

template<typename CSROrdinal, typename Value, typename HOrdType, typename LOrdType>
void twoLNumeric(Matrix<CSROrdinal, CSROrdinal, Value> & A,
        Matrix<CSROrdinal, CSROrdinal, Value> & B_org,
        TwoLevelMatrix<CSROrdinal, Value, HOrdType, LOrdType> & B,
        Matrix<CSROrdinal, CSROrdinal, Value> & C,
        UpperBounds<CSROrdinal, Value> & upb,
        CSROrdinal chunk_size) {

    CSROrdinal nrowsA = A.nrows();
    CSROrdinal * ptrA = A.Ptr();
    CSROrdinal * colindicesA = A.ColIndices();
    Value * valuesA = A.Values();

    CSROrdinal nrowsB = B_org.nrows();
    CSROrdinal ncolsB = B_org.ncols();
    CSROrdinal * ptrB = B_org.Ptr();
    CSROrdinal * colindicesB = B_org.ColIndices();
    Value * valuesB = B_org.Values();

    CSROrdinal * ptrC = C.Ptr();
    CSROrdinal * colindicesC = C.ColIndices();
    Value     * valuesC = C.Values();

    CSROrdinal widthH = ((long)1<<B.num_bits.higher);
    CSROrdinal widthL = ((long)1<<B.num_bits.lower);

    if (nrowsA < L2_CACHE_SIZE/sizeof(Value)) {
#pragma omp parallel
        {
            DenseHashMap<CSROrdinal, Value, Value, CSROrdinal> hashmap(ncolsB, upb.max_width, inf<Value>());
#pragma omp for schedule(dynamic, chunk_size)
            for (CSROrdinal v=0; v<nrowsA; ++v) {
                // if degree of v in A is less than one, then we can just add the degree
                if (ptrA[v+1] - ptrA[v] == 1) {
                    CSROrdinal u=colindicesA[ptrA[v]];
                    Value valvu=valuesA[ptrA[v]];
                    CSROrdinal pos_vw=ptrC[v];
                    for (CSROrdinal pos_w=ptrB[u]; pos_w<ptrB[u+1]; pos_w++, pos_vw++) {
                        colindicesC[pos_vw]=colindicesB[pos_w];
                        valuesC[pos_vw]=valuesB[pos_w]*valvu;
                    }
                    continue;
                } else if (ptrA[v+1] - ptrA[v] == 0)
                    continue;

                for (CSROrdinal u_pos=ptrA[v]; u_pos<ptrA[v+1]; ++u_pos) {
                    CSROrdinal u=colindicesA[u_pos];
                    Value valvu=valuesA[u_pos];
                    CSROrdinal start=ptrB[u];
                    CSROrdinal end=ptrB[u+1];
                    for (CSROrdinal w_pos=start; w_pos<end; ++w_pos) {
                        CSROrdinal w=colindicesB[w_pos];
                        Value valuw = valuesB[w_pos];
                        hashmap.insertInc(w, valvu*valuw);
                    }
                }
                CSROrdinal nonzerovw;
                Value valuevw;
                CSROrdinal posvw=ptrC[v];
                while(hashmap.getKeyValue(nonzerovw, valuevw)) {
                    colindicesC[posvw] = nonzerovw;
                    valuesC[posvw++] = valuevw;
                }
                hashmap.resetSize();
            }
        }
    } else {
        CSROrdinal * ptrBH = B.H.Ptr();
        HOrdType * colindicesBH = B.H.ColIndices();
        CSROrdinal * valuesBH = B.H.Values();
        Value   * valuesBL = B.values;
        LOrdType * colindicesBL = B.L;

#pragma omp parallel
        {
            SparseHashMap<CSROrdinal, Value, true, HOrdType> hashmapH(widthH, upb);
            DenseHashMap<CSROrdinal, Value, Value, LOrdType> hashmapL(widthL, widthL, inf<Value>());

#pragma omp for schedule(dynamic, chunk_size)
            for (CSROrdinal v=0; v<nrowsA; ++v) {

                // if degree of v in A is less than one, then we can just add the degree
                if (ptrA[v+1] - ptrA[v] == 1) {
                    CSROrdinal u=colindicesA[ptrA[v]];

                    // for each nonzero u-w in B
                    Value val = valuesA[ptrA[v]];
                    CSROrdinal j=ptrC[v];
                    for (CSROrdinal w_pos=ptrBH[u]; w_pos<ptrBH[u+1]; ++w_pos) {
                        CSROrdinal wH = (CSROrdinal)colindicesBH[w_pos];
                        for (CSROrdinal i=valuesBH[w_pos]; i<valuesBH[w_pos+1]; ++i, ++j) {
                            colindicesC[j] = (wH << B.num_bits.lower) + colindicesBL[i];
                            valuesC[j] = val * valuesBL[i];
                        }
                    }
                    continue;
                } else if (ptrA[v+1] - ptrA[v] == 0) {
                    continue;
                }

                // gather all the position info
                CSROrdinal u_start = ptrA[v];
                CSROrdinal u_end = ptrA[v+1];
                for (CSROrdinal u_pos=u_start; u_pos<u_end; ++u_pos) {
                    CSROrdinal u = colindicesA[u_pos];
                    Value val = valuesA[u_pos];
                    const CSROrdinal wH_start = ptrBH[u];
                    const CSROrdinal wH_end = ptrBH[u+1];
                    for (CSROrdinal wH_pos=wH_start; wH_pos<wH_end; ++wH_pos) {
                        HOrdType wH = colindicesBH[wH_pos];
                        hashmapH.insert(wH, wH_pos, val);
                    }
                }

                CSROrdinal v_sizeC = ptrC[v];
                CSROrdinal keynum = hashmapH.getKeysNum();
                for (CSROrdinal k=0; k<keynum; ++k) {
                    HOrdType nonzero;
                    Value      val;
                    hashmapH.getKey(k, nonzero, val);
                    CSROrdinal next;
                    CSROrdinal pos;
                    bool moreValue = hashmapH.getFirstValue(nonzero, pos, next);
                    CSROrdinal nonzero_wH = nonzero << B.num_bits.lower;
                    if (!moreValue) {
                        CSROrdinal start = valuesBH[pos];
                        CSROrdinal end = valuesBH[pos+1];
                        for (CSROrdinal w_pos=start; w_pos<end; ++w_pos, ++v_sizeC) {
                            colindicesC[v_sizeC] = nonzero_wH + colindicesBL[w_pos];
                            valuesC[v_sizeC] = val * valuesBL[w_pos];
                        }
                    } else {
                        do {
                            CSROrdinal start = valuesBH[pos];
                            CSROrdinal end = valuesBH[pos+1];
                            for (CSROrdinal k=start; k<end; ++k) {
                                hashmapL.insertInc(colindicesBL[k], val*valuesBL[k]);
                            }
                        } while (hashmapH.getCollisions(pos, next, val));

                        // read from hashmapL and write to C
                        LOrdType nonzero_wL;
                        Value valuec;
                        while(hashmapL.getKeyValue(nonzero_wL, valuec)) {
                            colindicesC[v_sizeC] = nonzero_wH + nonzero_wL;
                            valuesC[v_sizeC++] = valuec;
                        }
                        hashmapL.resetSize();
                    }
                }
                hashmapH.resetSize();
                assert(v_sizeC == ptrC[v+1]);
            }
        }
    }
}

template<typename CSROrdinal, typename Value, typename HOrdType, typename LOrdType>
void Multiplication(
        Matrix<CSROrdinal, CSROrdinal, Value> & A,
        Matrix<CSROrdinal, CSROrdinal, Value> & B_org,
        NumBits nbits,
        Matrix<CSROrdinal, CSROrdinal, Value> & C,
        bool verb = false,
        bool statistic_verb = false) {

    // Construct S, and check if use compression and if uses, stores compression
    // results in B as well
    if (verb) HooksRegionBegin("Construct S and get upperbound");
    TwoLevelMatrix<CSROrdinal, Value, HOrdType, LOrdType> B(nbits);
    constructHL<CSROrdinal, Value, HOrdType, LOrdType>(B_org, B, 512, statistic_verb);
    UpperBounds<CSROrdinal, Value> upb = getUpperBounds(A, B_org, B, 512, statistic_verb);
    if (verb) HooksRegionEnd("Construct S and get upperbound", STR_LEN, NUM_LEN);
    if (statistic_verb) upb.print();

    // allocate for ptr of C
    if (verb) HooksRegionBegin("MallocPtr");
    C.AllocatePtr(A.nrows() + 1);
    if (verb) HooksRegionEnd("MallocPtr", STR_LEN, NUM_LEN);

    // symbolic: compute ptr of C
    if (verb) HooksRegionBegin("Symbolic");
    CSROrdinal chunk_size = 0;
#pragma omp parallel
    {
        chunk_size = A.nrows() / omp_get_num_threads() ;
    }
    twoLSymbolic(A, B_org, B, C, upb, chunk_size / TASK_PER_THREAD_SYM > 0 ? chunk_size / TASK_PER_THREAD_SYM : 1);
    if (verb) HooksRegionEnd("Symbolic", STR_LEN, NUM_LEN);

    // allocate for col indices and values of C
    if (verb) HooksRegionBegin("MallocColIndices");
    C.AllocateColIndices(C.Ptr()[C.nrows()]);
    C.AllocateValues(C.Ptr()[C.nrows()]);
    if (verb) HooksRegionEnd("MallocColIndices", STR_LEN, NUM_LEN);

    // numeric: compute for col indices and values in C
    if (verb) HooksRegionBegin("Numeric");
    twoLNumeric(A, B_org, B, C, upb, chunk_size / TASK_PER_THREAD_NUM > 0 ? chunk_size / TASK_PER_THREAD_NUM : 1);
    if (verb) HooksRegionEnd("Numeric", STR_LEN, NUM_LEN);
}

// Multiplication where the maximum degree or A/B is 1.
template<typename CSROrdinal, typename Value>
void Multiplication_trival(
        Matrix<CSROrdinal, CSROrdinal, Value> & A,
        Matrix<CSROrdinal, CSROrdinal, Value> & B,
        Matrix<CSROrdinal, CSROrdinal, Value> & C,
        bool verb = false)
{
    if (verb) HooksRegionBegin("Construct S and get upperbound");
    // for trival case, we don't construct S, but to keep
    if (verb) HooksRegionEnd("Construct S and get upperbound", STR_LEN, NUM_LEN);

    if (verb) HooksRegionBegin("MallocPtr");
    // allocate for ptr
    CSROrdinal * ptrA = A.Ptr();
    CSROrdinal * colindicesA = A.ColIndices();
    Value * valuesA = A.Values();

    CSROrdinal nrowsB = B.nrows();
    CSROrdinal * ptrB = B.Ptr();
    CSROrdinal * colindicesB = B.ColIndices();
    Value * valuesB = B.Values();

    C.AllocatePtr(A.nrows() + 1);
    CSROrdinal * ptrC = C.Ptr();
    if (verb) HooksRegionEnd("MallocPtr", STR_LEN, NUM_LEN);

    if (verb) HooksRegionBegin("Symbolic");
    ptrC[0]=0;
    if (A.nrows() == A.nnz() && B.nrows() == B.nnz()) {
#pragma omp parallel for schedule (static)
        for (CSROrdinal v=0; v<A.nrows(); ++v)
            ptrC[v+1]=v+1;
    } else {
#pragma omp parallel for schedule (static)
        for (CSROrdinal v=0; v<A.nrows(); ++v) {
            if (ptrA[v+1] == ptrA[v]) {
                ptrC[v+1]=0;
                continue;
            }
            CSROrdinal u=colindicesA[ptrA[v]];
            ptrC[v+1]=ptrB[u+1]-ptrB[u];
        }
        std::partial_sum(ptrC, ptrC + C.nrows() + 1, ptrC);
    }
    if (verb) HooksRegionEnd("Symbolic", STR_LEN, NUM_LEN);

    if (verb) HooksRegionBegin("MallocColIndices");
    // allocate for col indices and values
    C.AllocateColIndices(C.Ptr()[C.nrows()]);
    C.AllocateValues(C.Ptr()[C.nrows()]);
    CSROrdinal * colindicesC = C.ColIndices();
    Value * valuesC = C.Values();
    if (verb) HooksRegionEnd("MallocColIndices", STR_LEN, NUM_LEN);

    if (verb) HooksRegionBegin("Numeric");
#pragma omp parallel for schedule (static)
    for (CSROrdinal v=0; v<A.nrows(); ++v) {
        if (ptrC[v+1] == ptrC[v]) continue;
        CSROrdinal u=colindicesA[ptrA[v]];
        colindicesC[ptrC[v]]=colindicesB[ptrB[u]];
        valuesC[ptrC[v]]=valuesA[v]*valuesB[ptrB[u]];
    }
    if (verb) HooksRegionEnd("Numeric", STR_LEN, NUM_LEN);
}

template<typename CSROrdinal, typename Value>
void SpGEMM (
        Matrix<CSROrdinal, CSROrdinal, Value> & A,
        Matrix<CSROrdinal, CSROrdinal, Value> & B,
        Matrix<CSROrdinal, CSROrdinal, Value> & C,
        bool verb = false,
        bool statistic_verb = false) {
    if (A.ncols() != B.nrows()) {
        std::cerr << "Error: the number of columns in A is not equal to the number of rows in B : "<< A.ncols() << ", " << B.nrows() << std::endl;
        std::exit( 1 );
    }

    C.init(A.nrows(), B.ncols());

    if (statistic_verb)
        std::cout << std::setw(STR_LEN) << "Statistics for SpGEMM" << " :" << std::endl;
    if ((A.maxDegree() == 1 && B.maxDegree() == 1) && !statistic_verb) {
        Multiplication_trival(A, B, C, verb);
        return;
    }

    NumBits nbits = getNumBits<CSROrdinal, Value>(B.ncols());
    if (statistic_verb) nbits.print();

    // assign different types for higher bits part and lower bits part
    bool uint16ForH = nbits.higher > 8 && nbits.higher <= 16;
    bool uint8ForH = nbits.higher <= 8;
    bool uint16ForL = nbits.lower > 8 && nbits.lower <= 16;
    bool uint8ForL = nbits.lower <= 8;
    if (uint8ForH && uint8ForL)
        Multiplication<CSROrdinal, Value, uint8_t, uint8_t>(A, B, nbits, C, verb, statistic_verb );
    else if (uint16ForH && uint16ForL)
        Multiplication<CSROrdinal, Value, uint16_t, uint16_t>(A, B, nbits, C, verb, statistic_verb );
    else if (uint16ForH && uint8ForL)
        Multiplication<CSROrdinal, Value, uint16_t, uint8_t>(A, B, nbits, C, verb, statistic_verb );
    else if (uint8ForH && uint16ForL)
        Multiplication<CSROrdinal, Value, uint8_t, uint16_t>(A, B, nbits, C, verb, statistic_verb );
    else if (uint8ForH)
        Multiplication<CSROrdinal, Value, uint8_t, CSROrdinal>(A, B, nbits, C, verb, statistic_verb );
    else if (uint16ForH)
        Multiplication<CSROrdinal, Value, uint16_t, CSROrdinal>(A, B, nbits, C, verb, statistic_verb );
    else if (uint8ForL)
        Multiplication<CSROrdinal, Value, CSROrdinal, uint8_t>(A, B, nbits, C, verb, statistic_verb );
    else if (uint16ForL)
        Multiplication<CSROrdinal, Value, CSROrdinal, uint16_t>(A, B, nbits, C, verb, statistic_verb );
    else
        Multiplication<CSROrdinal, Value, CSROrdinal, CSROrdinal>(A, B, nbits, C, verb, statistic_verb );
}
#endif
