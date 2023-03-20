# SparseCoLO

This GitHub repository is a "modern" presentation of the original SparseCoLO hosted at http://www.opt.c.titech.ac.jp/kojima/SparseCoLO/SparseCoLO.htm

SparseCoLO is a Matlab package for implementing the four conversion methods [1], proposed by Kim, Kojima, Mevissen and Yamashita, via positive semidefinite matrix completion for an optimization problem with matrix inequalities satisfying a sparse chordal graph structure. It is based on quite a general description of optimization problem including both primal and dual form of linear, semidefinite, second-order cone programs with equality/inequality constraints. Among the four conversion methods,  two methods utilize the domain-space sparsity  of a semidefinite matrix  variable  and the other two methods the range-space sparsity of a  linear matrix inequality (LMI) constraint of the given problem. SparseCoLO can be used as a preprocessor to reduce the size of the given problem before applying semidefinite programming solvers.

## References

[1] Sunyoung Kim, Masakazu Kojima, Martin Mevissen and Makoto Yamashita, "Exploiting sparsity in linear and nonlinear matrix inequalities via positive semidefinite matrix completion," Mathematical Programming, 129(1), 33–68. https://doi.org/10.1007/s10107-010-0402-6
