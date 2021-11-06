# CSeg

CSeg (Column-Segmented Sparse Matrix-Matrix Multiplication on Multicore CPUs) is
a OpenMP based SpGEMM, Sparse Matrix-Matrix Multiplication,
implementation.

💻 **Source Code:** [http://github.com/GT-TDAlab/CSeg]

CSeg is developed by the members of [GT-TDAlab](http://tda.gatech.edu).

## License

CSeg is distributed under BSD License. For details, see [`LICENSE`](LICENSE.md).

The [PIGO](https://github.com/GT-TDAlab/PIGO) library is included as an external program. See [`PIGO LICENSE`](./include/PIGO_LICENSE.md) for details.

## Contributors

- [Xiaojing An](https://xiaojingan.com/)
- [Umit V. Catalyurek](http://cc.gatech.edu/~umit)

## Contact

For questions or supports [open an issue](../../issues) or contact contributors via <tdalab@cc.gatech.edu>.

## Citation
Citation for the CSeg (BibTeX):

```bibtex
    @inproceedings{An21-HiPC,
        title = {Column-Segmented Sparse Matrix-Matrix Multiplication},
        author = {Xiaojing An and \"Umit V. \c{C}ataly\"{u}rek},
        year = {2021},
        month = {Dec},
        booktitle = {HiPC21: 28th IEEE International Conference on High Performance Computing, Data, & Analytics},
        keywords = {SpGEMM}
    }
```

## Compilers tested

- icpc version 2021.2.0 (gcc version 9.2.0 compatibility)
- gcc version 9.2.0 (GCC)

## How to build

Run `make` from current folder:

    make

## How to run

You can use `main` executable to run CSeg from command line.
We are using PIGO for parallel input, and currently support `mtx` inputs only.

Use `main -h` to see all options. An usage example, multiplying the matrices in ./exampleData

```shell
./main -A ./exampleData/gre_185_left.mtx -B ./exampleData/gre_185_right.mtx -t -w 1 -e 2
-> Loading inputs.
-> Matrix A at ./exampleData/gre_185_left.mtx statistics
                             Num of rows :        186
                             Num of cols :        201
                        Num of non-zeros :       1005
                              Max degree :          7
-> Matrix B at ./exampleData/gre_185_right.mtx statistics
                             Num of rows :        201
                             Num of cols :        211
                        Num of non-zeros :       1005
                              Max degree :          7

-> Warmup runs of SpGEMM
                   Statistics for SpGEMM :
                    Number of bits in Bh :          8
                    Number of bits in Bl :          0
                    Compression Ratio is :   0.483582, and compression is applied.
                         S has non-zeros :       1005, direct merging on dense array will be applied
   Construct S and get upperbound (secs) :   0.014287
                           max_widthComp :         26
                              max_widthH :         47
                               max_width :         47
                         max_collisionsH :         40
                                max_degA :          7
                                   flops :       5609
                        MallocPtr (secs) :   0.000001
                         Symbolic (secs) :   0.006328
                 MallocColIndices (secs) :   0.000008
                          Numeric (secs) :   0.006927
           Total runtime (warmup) (secs) :   0.027652


-> Perf runs of SpGEMM
   Construct S and get upperbound (secs) :   0.012966
                        MallocPtr (secs) :   0.000000
                         Symbolic (secs) :   0.010833
                 MallocColIndices (secs) :   0.000000
                          Numeric (secs) :   0.008451
                    Total runtime (secs) :   0.032293

   Construct S and get upperbound (secs) :   0.010876
                        MallocPtr (secs) :   0.000000
                         Symbolic (secs) :   0.011075
                 MallocColIndices (secs) :   0.000001
                          Numeric (secs) :   0.006498
                    Total runtime (secs) :   0.028485


-> Matrix C statistics
                             Num of rows :        186
                             Num of cols :        211
                        Num of non-zeros :       3093
                              Max degree :         24
```
