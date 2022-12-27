#!/bin/bash

curl https://cf.10xgenomics.com/samples/cell-exp/6.0.0/SC3_v3_NextGem_DI_CRISPR_A549_5K_SC3_v3_NextGem_DI_CRISPR_A549_5K/SC3_v3_NextGem_DI_CRISPR_A549_5K_SC3_v3_NextGem_DI_CRISPR_A549_5K_count_sample_feature_bc_matrix.h5 > crispr_6.0.0-tenx.h5

curl https://cf.10xgenomics.com/samples/cell-exp/6.0.0/SC3_v3_NextGem_DI_CRISPR_A549_5K_SC3_v3_NextGem_DI_CRISPR_A549_5K/SC3_v3_NextGem_DI_CRISPR_A549_5K_SC3_v3_NextGem_DI_CRISPR_A549_5K_count_sample_feature_bc_matrix.tar.gz > crispr-mtx.tar.gz
tar -xvf crispr-mtx.tar.gz
rm crispr-mtx.tar.gz

mv sample_feature_bc_matrix/barcodes.tsv.gz crispr_6.0.0-barcodes.tsv.gz
mv sample_feature_bc_matrix/features.tsv.gz crispr_6.0.0-features.tsv.gz
mv sample_feature_bc_matrix/matrix.mtx.gz crispr_6.0.0-matrix.mtx.gz
