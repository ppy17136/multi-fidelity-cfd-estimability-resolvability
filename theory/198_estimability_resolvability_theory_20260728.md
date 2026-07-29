# Formal two-gate estimability-resolvability framework

Date: 2026-07-28

## 1. Setting

Let the numerical observations satisfy the linear mean model

\[
\mathbf y = \mathbf X\boldsymbol\beta + \boldsymbol\varepsilon_{\rm num},
\]

where \(\mathbf X\in\mathbb R^{n\times p}\) is fixed by the acquired
multi-fidelity support, \(\boldsymbol\beta\) contains the model coefficients,
and \(\boldsymbol\varepsilon_{\rm num}\) denotes numerical perturbations rather
than experimental sampling noise.

The scientific claim is represented by a predeclared linear contrast

\[
\theta = \mathbf l^\mathsf T\boldsymbol\beta.
\]

The two gates answer different questions:

1. **Contrast estimability:** Is \(\theta\) uniquely determined by the acquired
   support?
2. **Numerical resolvability:** If estimable, is the magnitude of its estimate
   separated from a predeclared numerical evidence floor?

## 2. Gate 1 proposition

### Proposition 1: row-space criterion

For the linear mean model above, the contrast
\(\theta=\mathbf l^\mathsf T\boldsymbol\beta\) is estimable if and only if
\(\mathbf l^\mathsf T\) lies in the row space of \(\mathbf X\). Equivalently,
there exists a vector \(\mathbf a\in\mathbb R^n\) such that

\[
\mathbf a^\mathsf T\mathbf X=\mathbf l^\mathsf T.
\]

When this condition holds,
\(\widehat\theta=\mathbf a^\mathsf T\mathbf y\) has expectation
\(\theta\) under a zero-centred perturbation model, and its mean component is
independent of which generalized inverse is used to represent
\(\boldsymbol\beta\).

### Proof

If \(\mathbf l^\mathsf T=\mathbf a^\mathsf T\mathbf X\), then

\[
\mathbf a^\mathsf T E(\mathbf y)
=\mathbf a^\mathsf T\mathbf X\boldsymbol\beta
=\mathbf l^\mathsf T\boldsymbol\beta
=\theta,
\]

so the contrast is estimable.

Conversely, suppose \(\mathbf l^\mathsf T\) is not in the row space of
\(\mathbf X\). By the fundamental theorem of linear algebra, there exists
\(\mathbf v\in\operatorname{null}(\mathbf X)\) for which
\(\mathbf l^\mathsf T\mathbf v\ne0\). The coefficient vectors
\(\boldsymbol\beta\) and \(\boldsymbol\beta+\mathbf v\) then produce the same
mean response,

\[
\mathbf X(\boldsymbol\beta+\mathbf v)
=\mathbf X\boldsymbol\beta,
\]

but different contrast values,

\[
\mathbf l^\mathsf T(\boldsymbol\beta+\mathbf v)
\ne \mathbf l^\mathsf T\boldsymbol\beta.
\]

Therefore the observed support cannot determine the contrast uniquely, and the
contrast is not estimable. QED.

### Numerical implementation

Let \(\mathbf X^+\) be the Moore-Penrose pseudoinverse. The row-space
projection of \(\mathbf l\) is

\[
\mathbf l_{\rm proj}=\mathbf X^+\mathbf X\mathbf l.
\]

Gate 1 passes when

\[
\frac{\|\mathbf l-\mathbf l_{\rm proj}\|_2}
     {\max(\|\mathbf l\|_2,1)}
\le \tau_{\rm rank},
\]

with a declared numerical tolerance such as
\(\tau_{\rm rank}=10^{-10}\).

## 3. Two-by-two fidelity lattice

For binary fidelity factors \(A,B\in\{0,1\}\), use

\[
y_{AB}=\beta_0+\beta_A A+\beta_B B+\beta_{AB}AB.
\]

The complete design matrix, ordered as \(00,10,01,11\), is

\[
\mathbf X_{\rm full}=
\begin{bmatrix}
1&0&0&0\\
1&1&0&0\\
1&0&1&0\\
1&1&1&1
\end{bmatrix}.
\]

The interaction contrast is

\[
\Delta_{AB}
=y_{11}-y_{10}-y_{01}+y_{00}
=\beta_{AB},
\]

so \(\mathbf l=(0,0,0,1)^\mathsf T\) in coefficient space. The complete matrix
has rank four and the contrast is estimable.

If the \(11\) corner is absent,

\[
\mathbf X_{\rm miss}=
\begin{bmatrix}
1&0&0&0\\
1&1&0&0\\
1&0&1&0
\end{bmatrix}.
\]

Its fourth column is zero. The vector \(\mathbf l\) is not in its row space,
and \(\beta_{AB}\) is not estimable. An emulator may still output a value at
the missing corner, but that value is imposed by model structure or prior
assumptions rather than identified by the acquired factorial support.

## 4. Gate 2 definition

Let \(\widehat\Delta\) be an estimable target contrast and let
\(F_{\rm num}>0\) be a predeclared, response-aligned numerical evidence floor.
Define

\[
R=\frac{|\widehat\Delta|}{F_{\rm num}}.
\]

For a predeclared threshold \(\gamma\), gate 2 passes if

\[
R>\gamma.
\]

The current project uses \(\gamma=3\). This threshold is a materiality rule,
not a \(p\)-value, confidence level, or universal physical constant.

## 5. Conservative bounded-error guarantee

For a four-cell contrast, suppose every cell response has a deterministic
numerical perturbation bounded by

\[
|\varepsilon_{ab}|\le q.
\]

Then the contrast perturbation obeys

\[
|\varepsilon_{11}-\varepsilon_{10}-\varepsilon_{01}
 +\varepsilon_{00}|
\le 4q
\]

by the triangle inequality. Therefore \(F_{\rm num}=4q\) is a conservative
response-aligned floor.

If the true interaction is \(\Delta\), the observed contrast satisfies

\[
|\widehat\Delta-\Delta|\le F_{\rm num}.
\]

Consequently:

- if \(|\Delta|\le(\gamma-1)F_{\rm num}\), gate 2 cannot pass under the bounded
  perturbation model;
- if \(|\Delta|>(\gamma+1)F_{\rm num}\), gate 2 must pass;
- values between these bounds form a transition region in which the decision
  can depend on the realised numerical perturbation.

For \(\gamma=3\), a true interaction of \(2F_{\rm num}\) is a guaranteed
non-pass, while a true interaction greater than \(4F_{\rm num}\) is a
guaranteed pass.

## 6. Interpretation

Gate 1 is an algebraic statement about support. Gate 2 is an evidence statement
about magnitude relative to audited numerical variation. Passing gate 1 does
not imply passing gate 2.

Passing both gates justifies the bounded statement:

> The predeclared contrast is estimable from the acquired design and resolved
> relative to the audited numerical evidence floor.

It does not prove that the contrast has a unique physical cause, that all model
discrepancy has been separated, or that the response is free of unmodelled
numerical uncertainty.
