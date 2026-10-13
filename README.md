# Exact verification of a 34-commutator in five variables

This repository evaluates one explicit standard polynomial of polynomial vector
fields by exact integer subset recurrences. Its main result is

```text
s34(D1,...,D5, E11,E12,...,E15,E21,...,E55, Z1,Z2,Z3,Z4)(x5)
  = -3239792640.
```

Here `Di = ∂/∂xi`, `Eab = xa ∂/∂xb`, and `Zi = xi^2 ∂/∂xi`.
The result is a constant polynomial, so it is also its value at the origin.

## Exact conventions and witness

Work in `Q[x1,x2,x3,x4,x5]`. The 34 inputs occur in three consecutive blocks:

1. `D1, D2, D3, D4, D5`;
2. all 25 `Eab`, with `a=1,...,5` outermost and `b=1,...,5` innermost;
3. `Z1, Z2, Z3, Z4`.

The test polynomial is `x5`. The quadratic coefficients are ordinary monomials
`xi^2`, without division by `2!`. Define

```text
sk(X1,...,Xk) = sum over permutations σ of sign(σ) Xσ(1)...Xσ(k).
```

Composition means `XY(f)=X(Y(f))`: the rightmost field acts first. There is no
normalizing factor `1/k!`. The input ordering fixes the sign of the displayed
integer. No expected answer or algebraic multiplier is used in the evaluator.

## Build and reproduce

Use a 64-bit platform and a C++17 compiler supporting signed `__int128` and the
GCC/Clang bit-counting builtins. GCC and Clang on supported 64-bit targets provide
these extensions. This is not portable ISO C++17 code for every compiler.

```sh
make CXX=clang++
make check
make verify-left
```

`make check` compares both recurrence directions with literal permutation sums
on small inputs. `make verify-left` runs the full 34-input calculation and fails
if its output differs from the value above. `make verify-right` performs the
other recurrence, and `make verify` runs both sequentially. Fresh output and
progress logs go in `runs/`; distributed results in `results/` are not overwritten.

The equivalent direct build and full run are:

```sh
clang++ -O3 -std=c++17 direct_dense.cpp -o direct_dense
./direct_dense 5 Y left
./direct_dense 5 Y right
```

The evaluator accepts ranks 2 through 5, family `Y` or `C`, and direction `left`
or `right`. The main claim here uses `5 Y`; the retained `C` option evaluates a
different family and is not needed for this verification. `-march=native` is
optional; the source includes portable bit extraction when BMI2 is unavailable.

The rank-five runs recorded here took approximately 227 seconds (left) and
190 seconds (right). The memo array and index allocate 2,261,227,524 bytes in
total, about 2.11 GiB; allow 3 GiB of available memory to leave headroom. These are
allocation figures, not a measured peak resident-memory bound. Runtime depends
on the compiler and machine. Run the two large calculations sequentially.

## Why the recurrences compute the claimed quantity

Write each input as `Xi=x^{αi} D_{di}` and put `δi=αi-e_{di}`. Then

```text
Xi(x^μ) = μ[di] x^(μ+δi).
```

For any fixed subset, every surviving composition has the same final exponent.
Its alternating coefficient can therefore be stored as one integer. Let `S`
be an ordered subset of the input indices and let `i` be its `j`th member in
increasing order. All signs below use this relative position, not the original
absolute index. The sum of the 34 exponent shifts is `-e5`, which proves the
final output is constant.

### Grouping by the leftmost input

Set `ν(S)=e5+sum_{i∈S} δi`. Let `L(S)` be the coefficient of `x^ν(S)` in the
alternator on `S` applied to `x5`. Moving `i` to the leftmost position gives

```text
L(empty) = 1,
L(S) = sum_{i∈S} (-1)^(j-1) ν(S\{i})[di] L(S\{i}).
```

A term with a negative exponent or zero derivative factor contributes zero.
The initial exponent is `ν(full)=0`. The implementation removes `δi` to obtain
the exponent of the child state.

### Grouping by the rightmost input

For a remaining subset `S`, set `μ(S)=-sum_{i∈S} δi`. Thus `μ(full)=e5`.
Let `R(S)` be the constant coefficient obtained by applying the alternator on
`S` to `x^μ(S)`. Moving `i` to the rightmost position gives

```text
R(empty) = 1,
R(S) = sum_{i∈S} (-1)^(|S|-j) μ(S)[di] R(S\{i}).
```

After this rightmost action, the exponent is
`μ(S)+δi=μ(S\{i})`. The same vanishing rules apply. These formulas are obtained
by partitioning the defining permutation sum; neither uses a conjectural
vanishing identity. Both methods have the same witness and share storage code,
so they are not wholly independent implementations.

## Storage and exact arithmetic

The five diagonal fields `Eii` have zero exponent shift. Their bits can be
omitted when deciding whether an exponent is nonnegative, but remain in every
memo key. The other 29 bits specify the five constants, twenty off-diagonal
linear fields, and four quadratics. A bit table marks valid choices; prefix
population counts assign each valid choice a unique rank. The memo index is
`(rank << 5) | diagonal_bits`. This is an injective encoding of valid subsets,
not a quotient identifying distinct inputs. The code additionally checks
100,000 reproducible sampled table entries against the original field list.

The input monomial has degree one and there are four quadratic fields. Any
surviving intermediate monomial has degree at most five, so each derivative
factor is at most five. The evaluator checks this bound. Child values are
checked signed 64-bit integers. Each sum has at most 34 terms and uses signed
128-bit arithmetic: even the bound `34*5*2^63 < 2^71` is safely within its range.
Every completed sum is range-checked before conversion to 64 bits.

The memo array stores 32-bit values where possible; larger values and reserved
sentinel values use a separate map keyed by the full subset. All coefficient
arithmetic is integer arithmetic. Floating point is used only for timing.

## Tests and recorded results

`test_dense_literal.cpp` enumerates permutations directly, computes permutation
parity, and applies monomial fields from right to left. Its 160 deterministic
random lists have arities 2 through 8, with 18 nonzero results. Both recurrences
agree with the literal sum; 320 repeated calls also agree. Those repeats test
repeatability, not hash collisions or independent algorithms. A fixed test
checks

```text
s4(D1,D2,x1 D1,x2 D1)(x1^2) = -4.
```

The distributed full results were regenerated on 13 October 2026 using the exact
`direct_dense.cpp` included here. See `results/Y5_left.json`,
`results/Y5_right.json`, and `results/environment.txt`. The left run evaluated
499,457,290 nonempty subset states; the right run evaluated 468,931,939. Both
returned `-3239792640`. The small-test result is in
`results/dense_literal_checks.json`. Small literal tests check the recurrence
implementation; they do not enumerate the `34!` terms of the full witness.

## Mathematical scope

With `N=n^2+2n-2`, rank five has `N=33` and `N+1=34`. This computation disproves
the unrestricted assertion `s_{N+1}=0`. It does not refute the different
conjecture that `s_N` is nonzero. Indeed, the standard recurrence shows that
nonvanishing of `s34` implies nonvanishing of `s33` in rank five.

The all-rank critical-vanishing assertion appears in the remark after Lemma 7.4
of [Dzhumadil'daev's 2002 preprint](https://arxiv.org/abs/math/0203036v2) and in
Section 1.5 of [Lebedev and Leites](https://arxiv.org/abs/1906.01380).
[Dzhumadil'daev's 2004 Theorem 2.1(ii)](https://doi.org/10.1007/s00014-004-0807-2)
gives the higher bound `s35=0` on `W5`. Combining that bound with this witness
shows that 35 is the least standard-identity arity on `W5`. These are literature
connections; the evaluator does not assume these results. Rank-three `s14=0`
is known from [the 2008 Theorem 1](https://doi.org/10.2991/jnmp.2008.15.1.7), so
this example does not imply critical nonvanishing for every odd rank.

## Citation and release

Author: Nurtas Shyntas. Title: *Exact verification of a 34-commutator in five
variables*. Citation metadata is supplied in `CITATION.cff`. Cite the actual
release or commit used for reproduction. No archival DOI or release has yet
been issued. The repository is
[nurtasshyntas/n-commutator-verifier](https://github.com/nurtasshyntas/n-commutator-verifier).
