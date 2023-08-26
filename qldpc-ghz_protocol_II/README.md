# qldpc-ghz_protocol_II

Code base for the simpler protocol to distill GHZ states using CSS codes -- Protocol II of https://arxiv.org/abs/2210.14143

Copyright (C) 2022 Narayanan Rengaswamy

This project is licensed under the terms of the GNU Affero General Public License (AGPL) v3.0. See LICENSE.md for details.

If you use this code base, please add a citation to the paper as well as the archived version on Zenodo (https://zenodo.org/record/8284904)

**Scripts**:

*LP_MSA_dep.m*: Standard error correction simulation of Lifted Product QLDPC code family with vanilla min-sum algorithm (MSA) decoder; this can be reinterpreted as simulating the Bell pair distillation protocol in the above paper according to the noise setting

*LP_ghz_simple.m*: Generate performance curves for GHZ distillation with Lifted Product QLDPC code family using Protocol II in the paper; this uses the same vanilla MSA decoder which is different from the optimized decoder used for the results in the paper

**Functions**:

*concatenate_pauli.m*: Concatenate the binary representations of two Pauli matrices

*find_logical_paulis_by_ghz_msmt.m*: Compute logical X and Z generators for a stabilizer code using stabilizer measurements on GHZ states

*gf2matinv.m*: MATLAB in-built function for inverting binary matrices

*gf2rref.m*: Modified version of MATLAB in-built function for performing reduced row echelon form for binary matrices

*pauli_multiply.m*: Find the product of two Pauli matrices

*paulis_multiply.m*: Find the product of multiple Pauli matrices

*print_matrix_with_spaces.m*: Print the rows of matrix with spaces and horizontal lines between specific columns and rows, respectively

*stabilizer_formalism_msmt.m*: Simulate the measurement of Pauli operators on a stabilizer state using the stabilizer formalism (https://arxiv.org/abs/quant-ph/9807006)

*syndrome_MSA.m*: Vanilla syndrome-based min-sum algorithm (MSA) decoder

**Figures**:

*LP118_12_16_20_MSAseqvars80.fig*: MATLAB figure corresponding to Fig. 1 (left) in the paper

*LP118_12_16_20_MSAseqvars80_threshold.fig*: MATLAB figure corresponding to Fig. 1 (right) in the paper

*GHZsimple_LP118_12_16_20_MSAseqvars80_2.fig*: MATLAB figure corresponding to Fig. 4 (left) in the paper

*GHZsimple_LP118_12_16_20_MSAseqvars80_threshold_2.fig*: MATLAB figure corresponding to Fig. 4 (right) in the paper
