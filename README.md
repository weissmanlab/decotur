For questions please email rohan.mehta@elmhurst.edu.

# DeCoTUR
## Detecting Coevolving Traits Using Relatives
This package computes the "coevolution score" presented in [this manuscript](https://link.springer.com/article/10.1186/s12859-023-05363-4).

## Updates
3/7/2023 - Removed a line of debug code that I missed earlier.

2/2/2023 - Updated all code in preparation for manuscript submission

6/1/2022 - Changed the significance test. It is not longer based on Fisher's exact test, but it is now based on a Poisson Binomial distribution. This change better reflects the fact that not every closely related pair has an equivalent chance of producing a signal.

6/1/2022 - Fixed an issue in an internal function which wasn't handling block sizes properly.

## Installation
This package is implemented in R. If you don't have `devtools`, then run
`install.packages('devtools')`. If you don't have an updated Rtools installation, install it from 
`https://cran.r-project.org/bin/windows/Rtools/` or `https://www.r-project.org/nosvn/winutf8/ucrt3/`.

To install DeCoTUR, run

`library(devtools)`

`install_github("weissmanlab/decotur")`

You may also have to install the library `PoissonBinomial` via `install.packages('PoissonBinomial')`.

To load DeCoTUR, run
`library(DeCoTUR)`
`library(PoissonBinomial)`

## Using DeCoTUR
The primary function you would want to use is

`get_scores`

Everything else is a helper function. The function `get_scores` has many arguments. The two most important are the presence-absence matrix and the distance matrix.

For the presence-absence matrix, each row is a trait, and each column is a sample. For example

| Trait | Sample1 | Sample2 | ... | SampleN |
| ----- | ------- | ------- | --- | -------- |
| Trait1 | 1 | 0 | ... | 1 |
| Trait2 | 1 | 1 | ... | 0 |
| Trait3 | 1 | 1 | ... | 1 |
| ... | ... | ... | ... | ... |
| TraitM | 0 | 0 | ... | 1 |

For the distance matrix, each row and column is a sample. For example

| | Sample1 | Sample2 | ... | SampleN |
| ----- | ------- | ------- | --- | -------- |
| Sample1 | 0 | 0.001 | ... | 1.01 |
| Sample2 | 0.001 | 0 | ... | 0.023 |
| Sample3 | 6.07 | 1.97 | ... | 1.45 |
| ... | ... | ... | ... | ... |
| SampleN | 1.01 | 0.023 | ... | 0 |

The other inputs to this function are: A) a method for determining the close-pair cutoff. This can be 
1. Auto. It looks at the pairwise distance distribution and tries to find a cutoff that encompasses the first peak. Not recommended for first time use.
2. Distance. You give it a distance cutoff, it uses it.
3. Fraction. You give it a fraction (i.e. I want the 5% closest pairs, so my fraction is 0.05).
4. Fixed Number. I.e. I want 500 close pairs.

Each method requires particular parameters detailed in the function descriptions. Put those parameters in a list, and that's another input into the function. For instance, if you want to use the Fixed Number method, that requires 1) a number of close pairs (let's say 500), a random seed (let's say 1), and whether or not you want to see the histogram of the distance distribution (which I would recommend), so it's TRUE.

So for the closepair_params part of the argument, this example would need 
`closepair_params = list(500, 1, TRUE)`

B) A block size for computation. If you divide up the traits you want to evaluate into blocks, the computation goes much faster. Provided that you had a lot of traits to being with. If not, set this number to the number of traits.

C) Whether or not you want to downweight bushes to avoid a single bush dominating the score.

D) Whether or not you want p values with your results. Default is TRUE, though my current implementation of this is fairly slow so you might want to set it to FALSE initially. I'm currently working on speeding it up using R packages that can handle very large matrices.

E) Whether or not you want text output for progress. I would recommend keeping this at the default of TRUE.

F) Whether or not you want to run the relatively "speed" optimized version or the relatively "memory" optimized version. In my experience, there is no reason to choose "memory" over "speed" so long as you have a small enough "blocksize", but the option is there in case you run into memory issues. I am currently working on using R packages that can handle very large matrices to obviate the need for the "memory" option.

## Troubleshooting

If you encounter the error: "PoissonBinomial.cpp:4:10: fatal error: fftw3.h: No such file or directory" upon installing the PoissonBinomial dependency on a Linux machine, this error can be resolved by first running 
`sudo apt-get install libfftw3-dev`.
