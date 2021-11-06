/**
 * CSeg (Column-Segmented Sparse Matrix-Matrix Multiplication on Multicore CPUs)
 *
 * Copyright (c) 2021, GT-TDAlab (Umit V. Catalyurek)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#include <iostream>

#include "utils.h" // to use Config, and parseCommandLine
#include "spgemm.hpp" // to call SpGEMM, and enable the usage of class Matrix as well.

void run(Config c)
{
    // input reading
    std::cout << "-> Loading inputs." << std::endl;
    Matrix<CSROrdinal, CSROrdinal, Value> A( c.A );
    A.PrintStatis(STR_LEN, NUM_LEN, "A", c.A);
    Matrix<CSROrdinal, CSROrdinal, Value> B;
    if (c.B == c.A) B.init( A );
    else B.init( c.B );
    B.PrintStatis(STR_LEN, NUM_LEN, "B", c.B);

    Matrix<CSROrdinal, CSROrdinal, Value> C;

    // call SpGEMM based on the number of warm-up iterations
    for (CSROrdinal i=0; i<c.InterationsWarmUp; ++i) {
        if (i==0)
            std::cout << std::endl << "-> Warmup runs of SpGEMM" << std::endl;

        HooksRegionBegin("Total runtime (warmup)");
        SpGEMM(A, B, C, c.TrackIndividualTime, i==0);
        HooksRegionEnd("Total runtime (warmup)", STR_LEN, NUM_LEN, c.TrackIndividualTime);

        C.Free();
    }

    // call SpGEMM based on the number of execution iteration
    for (CSROrdinal i=0; i<c.InterationsExecution; ++i) {
        if (i==0)
            std::cout << std::endl << "-> Perf runs of SpGEMM" << std::endl;

        HooksRegionBegin("Total runtime");
        SpGEMM(A, B, C, c.TrackIndividualTime);
        HooksRegionEnd("Total runtime", STR_LEN, NUM_LEN, c.TrackIndividualTime);

        if (i<c.InterationsExecution-1) C.Free();
    }

    if (c.InterationsExecution != 0) {
        std::cout << std::endl;
        C.PrintStatis(STR_LEN, NUM_LEN, "C", c.C);
    }

    // output writing
    if (c.C == "")
        C.Free();
    else
        C.Serialize(c.C);
}

int main(int argc, char** argv) {

    Config c = parseCommandLine(argc, argv);

    run(c);

    return 0;
}
