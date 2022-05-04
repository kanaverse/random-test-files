#!/bin/bash

curl -L https://cf.10xgenomics.com/samples/cell-vdj/3.0.0/vdj_v1_hs_pbmc2_5gex_protein/vdj_v1_hs_pbmc2_5gex_protein_filtered_feature_bc_matrix.h5 > immune_3.0.0-tenx.h5

curl -L https://cf.10xgenomics.com/samples/cell-vdj/3.0.0/vdj_v1_hs_pbmc2_5gex_protein/vdj_v1_hs_pbmc2_5gex_protein_filtered_feature_bc_matrix.tar.gz > immune-mtx.tar.gz
tar -xvf immune-mtx.tar.gz
rm immune-mtx.tar.gz

mv filtered_feature_bc_matrix/barcodes.tsv.gz immune_3.0.0-barcodes.tsv.gz
mv filtered_feature_bc_matrix/features.tsv.gz immune_3.0.0-features.tsv.gz
mv filtered_feature_bc_matrix/matrix.mtx.gz immune_3.0.0-matrix.mtx.gz
