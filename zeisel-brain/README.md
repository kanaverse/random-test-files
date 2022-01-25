# Zeisel brain single-cell RNA-seq dataset

Pretty much as it says on the tin.
This takes the classic 3000-cell dataset from [Zeisel _et al._ (2015)](https://pubmed.ncbi.nlm.nih.gov/25700174/) and provides it in the HDF5 or MatrixMarket formats.
For HDF5, we use a dense matrix (mostly for testing), the 10X HDF5 matrix format, or the H5AD format.
We also generate some subsets of each format for quickly testing applications.
Run `generate.R` to construct the files.
