#ifndef CSR_MATRIX_H_
#define CSR_MATRIX_H_

#include <cstddef>
#include <vector>
#include <string>

#include "pigo.hpp" // for matrix input and malloc/free for shared pointers
#include <omp.h>

class MtxHeader {
    private:
        bool hasValue_;
        bool symmetric_;
    public:
        MtxHeader(std::string & file_name) : hasValue_(false), symmetric_(false) {
            pigo::ROFile rf { file_name };
            std::string head { rf.fp(), std::min((size_t)64, rf.size()) };
            if (head.find("real") != std::string::npos) hasValue_ = true;
            if (head.find("complex") != std::string::npos) hasValue_ = true;
            if (head.find("symmetric") != std::string::npos) symmetric_ = true;
        }
        bool hasValue() { return hasValue_; }
        bool symmetric() { return symmetric_; }
};

template <class Ordinal, class Offset, class Value>
class Matrix {
  private:
    std::shared_ptr<Offset> ptr_;
    std::shared_ptr<Ordinal> colindices_;
    std::shared_ptr<Value> values_;

    Offset nrows_;
    Offset ncols_;
    Offset nnz_;
    Offset maxDegree_;
    bool has_value_;

    bool checkFileSuffix(std::string &file_name) {
        auto endswith = [](std::string & str, std::string & suffix) {
            if (str.length() < suffix.length()) return false;
            return (0 == str.compare(str.length() - suffix.length(), suffix.length(), suffix));
        };

        // Currently, we only process mtx format
        return endswith(file_name, ".mtx");
    }

  public:

    // Constructor
    Matrix() : nrows_(0), ncols_(0), nnz_(0), maxDegree_ (0), has_value_ (false) {}
    Matrix(
        std::string &file_name,
        bool storeValue = true);

    void init(Offset nrows, Offset ncols) { nrows_ = nrows; ncols_ = ncols; }

    void init(Matrix & that) {
        ptr_ = that.PtrSP();
        colindices_ = that.ColIndicesSP();
        values_ = that.ValuesSP();
        nrows_ = that.nrows();
        ncols_ = that.ncols();
        nnz_ = that.nnz();
        maxDegree_ = that.maxDegree();
        has_value_ = that.hasValue();
    }

    void init(std::string &file_name,
        bool storeValue = true);

    void AllocatePtr(
      Offset s
    );

    void AllocateColIndices(
      Offset s
    );

    void AllocateValues(
      Offset s
    );

    Offset * Ptr() { return ptr_.get(); }

    Ordinal * ColIndices() { return colindices_.get(); }

    Value * Values() { return values_.get(); }

    std::shared_ptr<Offset> PtrSP() { return ptr_; }

    std::shared_ptr<Ordinal> ColIndicesSP() { return colindices_; }

    std::shared_ptr<Value> ValuesSP() { return values_; }

    Offset maxDegree();

    Offset nnz();

    Offset nrows() { return nrows_; }

    Offset ncols() { return ncols_; }

    bool hasValue() { return has_value_; }

    void Serialize(std::string outname);

    void PrintStatis(size_t str_len, size_t num_len, const std::string name, const std::string filename = "")
    {
        std::cout << "-> Matrix " + name << (filename == "" ? "" : (" at " + filename) ) << " statistics" << std::endl;
        std::cout << std::setw(str_len) << "Num of rows" << " : " << std::setw(num_len) << nrows_ << std::endl;
        std::cout << std::setw(str_len) << "Num of cols" << " : " << std::setw(num_len) << ncols_ << std::endl;
        std::cout << std::setw(str_len) << "Num of non-zeros" << " : " << std::setw(num_len) << nnz_ << std::endl;
        std::cout << std::setw(str_len) << "Max degree" << " : " << std::setw(num_len) << maxDegree() << std::endl;
    }

    ~Matrix() {
        Free();
        resetSize();
    }

    void Free() {
        ptr_.reset();
        colindices_.reset();
        values_.reset();
        has_value_ = false;
    }

    void resetSize() {
        nrows_=0;
        ncols_=0;
        nnz_=0;
        maxDegree_=0;
    }

    Matrix(const Matrix &) = delete;
    Matrix(Matrix &&) = delete;
    Matrix& operator=(const Matrix &) = delete;
    Matrix& operator=(Matrix &&) = delete;

};

// Method implementations
#include <cstdint>
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <fstream>

template <class Ordinal, class Offset, class Value>
void Matrix<Ordinal, Offset, Value>::init(
    std::string &file_name,
    bool storeValue) {
    maxDegree_ = 0;
    has_value_ = storeValue;

    assert(checkFileSuffix(file_name));

    MtxHeader h(file_name);

    if (h.hasValue() && storeValue) {
        if (h.symmetric()) {
            pigo::COO< Ordinal, Ordinal, std::shared_ptr<Ordinal>, true, false, false, true, Value, std::shared_ptr<Value>> coo { file_name };
            pigo::CSR< Ordinal, Offset, std::shared_ptr<Ordinal>, std::shared_ptr<Ordinal>, true, Value, std::shared_ptr<Value> > mtx { coo };
            mtx.sort();

            nrows_ = mtx.nrows();
            ncols_ = mtx.ncols();
            nnz_ = mtx.m();

            ptr_ = mtx.offsets();
            colindices_ = mtx.endpoints();
            values_ = mtx.weights();
        } else {
            pigo::CSR< Ordinal, Offset, std::shared_ptr<Ordinal>, std::shared_ptr<Ordinal>, true, Value, std::shared_ptr<Value> > mtx { file_name };
            mtx.sort();
            nrows_ = mtx.nrows();
            ncols_ = mtx.ncols();
            nnz_ = mtx.m();

            ptr_ = mtx.offsets();
            colindices_ = mtx.endpoints();
            values_ = mtx.weights();
        }
    } else {
        if (h.symmetric()) {
            pigo::COO< Ordinal, Ordinal, std::shared_ptr<Ordinal>, true, false, false, false, Value, std::shared_ptr<Value>> coo { file_name };
            pigo::CSR< Ordinal, Offset, std::shared_ptr<Ordinal>, std::shared_ptr<Ordinal>, false, Value, std::shared_ptr<Value> > mtx { coo };
            mtx.sort();

            nrows_ = mtx.nrows();
            ncols_ = mtx.ncols();
            nnz_ = mtx.m();

            ptr_ = mtx.offsets();
            colindices_ = mtx.endpoints();
        } else {
            pigo::CSR< Ordinal, Offset, std::shared_ptr<Ordinal>, std::shared_ptr<Ordinal>, false, Value, std::shared_ptr<Value> > mtx { file_name };
            mtx.sort();
            nrows_ = mtx.nrows();
            ncols_ = mtx.ncols();
            nnz_ = mtx.m();

            ptr_ = mtx.offsets();
            colindices_ = mtx.endpoints();
        }
        if (storeValue) {
            pigo::detail::allocate_mem_(values_, nnz_);
#pragma omp parallel for schedule (static)
            for (size_t i=0; i<nnz_; i++) values_.get()[i] = 1;
        }
    }
}

template<class Ordinal, class Offset, class Value>
Matrix<Ordinal, Offset, Value>::Matrix(
    std::string &file_name,
    bool storeValue) {
    init(file_name, storeValue);
}

template <class Ordinal, class Offset, class Value>
void Matrix<Ordinal, Offset, Value>::AllocatePtr(Offset s) {
  pigo::detail::allocate_mem_(ptr_, s);
  nrows_ = s-1;
  if (ncols_ == 0) ncols_ = nrows_;
}

template <class Ordinal, class Offset, class Value>
void Matrix<Ordinal, Offset, Value>::AllocateColIndices(Offset s){
  pigo::detail::allocate_mem_(colindices_, s);
  nnz_ = s;
}

template <class Ordinal, class Offset, class Value>
void Matrix<Ordinal, Offset, Value>::AllocateValues(Offset s){
  pigo::detail::allocate_mem_(values_, s);
  has_value_ = true;
}

template <class Ordinal, class Offset, class Value>
Offset Matrix<Ordinal, Offset, Value>::nnz() {

    if (nnz_ == 0) {
        if (ptr_ && nrows_ && ptr_.use_count() != 0)
            nnz_ = Ptr()[nrows_];
    }
    return nnz_;
}

template <class Ordinal, class Offset, class Value>
Offset Matrix<Ordinal, Offset, Value>::maxDegree(){
  if(maxDegree_ > 0) {
    return maxDegree_;
  }
  Offset max = 0, t=0;
  for(Offset i=0; i<nrows_; i++) {
    t = ptr_.get()[i+1]-ptr_.get()[i];
    max = (max < t) ? t : max;
  }
  maxDegree_ = max;
  return max;
}

template <class Ordinal, class Offset, class Value>
void Matrix<Ordinal, Offset, Value>::Serialize(
    std::string outname) {

  // TODO: use PIGO for this
  std::ofstream o ( outname, std::ios::out | std::ios::binary );

  o.write ((char*)&nrows_, sizeof (Ordinal));
  o.write ((char*)&ncols_, sizeof (Ordinal));
  o.write ((char*)&nnz_, sizeof (Offset));

  for( size_t i=0 ; i<=nrows() ; i++ ) {
      o.write( (char*)&ptr_.get()[i], sizeof (Offset) );
  }
  o.write ((char*)ptr_.get(), sizeof(Offset) * (nrows_ + 1));
  o.write ((char*)colindices_.get(), sizeof(Ordinal) * nnz_);
  o.write ((char*)values_.get(), sizeof(Value) * nnz_);

  o.close();
}

#endif // CSR_MATRIX_H_
