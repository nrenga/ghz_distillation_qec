# ghz_distillation_qec

Code base for the protocol to distill GHZ states using stabilizer codes

Copyright (C) 2021 Narayanan Rengaswamy

This project is licensed under the terms of the GNU Affero General Public License (AGPL) v3.0. See LICENSE.md for details.

**Scripts**:

*ghz_distillation_simulation.m*: Perform one round of the GHZ distillation protocol using a given stabilizer code

*ghz_distillation_test_performance.m*: Generate performance curves for GHZ distillation using a given stabilizer code

**Functions**:

*concatenate_pauli.m*: Concatenate the binary representations of two Pauli matrices

*find_binary_symmetric_matrix.m*: Find a binary symmetric matrix S satisfying A*S = B (mod 2)

*find_logical_paulis_by_ghz_msmt.m*: Compute logical X and Z generators for a stabilizer code using stabilizer measurements on GHZ states

*find_logical_paulis_by_bell_msmt.m*: Compute logical X and Z generators for a stabilizer code using stabilizer measurements on Bell states

*find_min_wt_error.m*: Find the minimum weight Pauli error matching the syndrome for a stabilizer code

*pauli_multiply.m*: Find the product of two Pauli matrices

*paulis_multiply.m*: Find the product of multiple Pauli matrices

*print_matrix_with_spaces.m*: Print the rows of matrix with spaces and horizontal lines between specific columns and rows, respectively

*stabilizer_code_on_pauli_channel.m*: Simulate a stabilizer code on a standard Pauli channel

*stabilizer_formalism_msmt.m*: Simulate the measurement of Pauli operators on a stabilizer state using the stabilizer formalism (https://arxiv.org/abs/quant-ph/9807006)

**Figures**:

*performance_5qubit_code.fig*: MATLAB figure showing the performance curves for the 5-qubit code (as in the paper)
