"""Exact endpoint semantics for finite stored binary floating-point inputs.

This arithmetic does not remove uncertainty in the scientific input values.
"""
from fractions import Fraction
import math


def rational(value):
    x = float(value)
    if not math.isfinite(x):
        raise ValueError("scalar inputs and thresholds must be finite")
    return Fraction.from_float(x)


def scalar_endpoints(z):
    rows = tuple(tuple(rational(v) for v in row) for row in z)
    coherent = tuple(sum((row[p] for row in rows), Fraction())
                     for p in range(len(rows[0])))
    return (min(coherent), max(coherent),
            sum((min(row) for row in rows), Fraction()),
            sum((max(row) for row in rows), Fraction()))


def exact_state(lower, upper, threshold):
    if lower > threshold:
        return "resolved_above"
    if upper < threshold:
        return "resolved_below"
    return "indeterminate"


def member_value(z, labels):
    return sum((rational(z[i, p]) for i, p in enumerate(labels)), Fraction())
