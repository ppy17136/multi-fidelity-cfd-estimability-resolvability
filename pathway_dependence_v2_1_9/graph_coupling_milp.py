"""MILP certificates for graph-constrained pathway-coupling relaxation.

SciPy's HiGHS backend is used as a transparent reference implementation.  The
formulations are binary Potts labellings: vertex variables choose pathway (or
ordered extremizer-pair) labels, and edge variables record disagreements.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix

from graph_constrained_coupling import coherent_scalar_range, interval_state
from scalar_semantics import rational, scalar_endpoints, exact_state, member_value


@dataclass(frozen=True)
class MilpCertificate:
    status: str
    objective: float | None
    dual_bound: float | None
    relative_gap: float | None
    node_count: int | None
    labels: tuple | None
    attained_value: float | None = None


def _edges(q: int, edges: Iterable[Sequence[float]]):
    checked = []
    seen = set()
    for raw in edges:
        if len(raw) != 3:
            raise ValueError("each edge must be (left, right, weight)")
        left, right, weight = int(raw[0]), int(raw[1]), float(raw[2])
        if not (0 <= left < q and 0 <= right < q) or left == right:
            raise ValueError("invalid edge endpoints")
        if weight < 0 or not np.isfinite(weight):
            raise ValueError("edge weights must be finite and nonnegative")
        key = tuple(sorted((left, right)))
        if key in seen:
            raise ValueError("parallel edges are not supported")
        seen.add(key)
        checked.append((key[0], key[1], weight))
    return tuple(checked)


def _require_connected(q, edges):
    adjacency = [[] for _ in range(q)]
    for left, right, _ in edges:
        adjacency[left].append(right)
        adjacency[right].append(left)
    reached = {0}
    stack = [0]
    while stack:
        node = stack.pop()
        for neighbour in adjacency[node]:
            if neighbour not in reached:
                reached.add(neighbour)
                stack.append(neighbour)
    if len(reached) != q:
        raise ValueError("decision-relaxation graphs must be connected")


def _solve_potts(
    allowed_labels,
    edges,
    *,
    contribution_values=None,
    threshold=None,
    threshold_sense=None,
    time_limit=None,
):
    allowed = tuple(tuple(labels) for labels in allowed_labels)
    q = len(allowed)
    if not allowed or any(not labels for labels in allowed):
        raise ValueError("every vertex must have at least one allowed label")
    labels = tuple(sorted(set().union(*map(set, allowed))))
    label_index = {label: j for j, label in enumerate(labels)}
    label_count = len(labels)
    checked_edges = _edges(q, edges)
    x_count = q * label_count
    variable_count = x_count + len(checked_edges)

    objective = np.zeros(variable_count)
    for edge_index, (_, _, weight) in enumerate(checked_edges):
        objective[x_count + edge_index] = weight
    lower_bounds = np.zeros(variable_count)
    upper_bounds = np.ones(variable_count)
    for i, choices in enumerate(allowed):
        permitted = set(choices)
        for label, j in label_index.items():
            if label not in permitted:
                upper_bounds[i * label_count + j] = 0.0

    row_indices = []
    column_indices = []
    values = []
    row_lower = []
    row_upper = []

    def append_row(entries, lower, upper):
        row = len(row_lower)
        for column, value in entries:
            row_indices.append(row)
            column_indices.append(column)
            values.append(value)
        row_lower.append(lower)
        row_upper.append(upper)

    for i in range(q):
        append_row(
            [(i * label_count + j, 1.0) for j in range(label_count)],
            1.0,
            1.0,
        )

    for edge_index, (left, right, _) in enumerate(checked_edges):
        y = x_count + edge_index
        for j in range(label_count):
            left_x = left * label_count + j
            right_x = right * label_count + j
            append_row(((y, 1.0), (left_x, -1.0), (right_x, 1.0)), 0.0, np.inf)
            append_row(((y, 1.0), (left_x, 1.0), (right_x, -1.0)), 0.0, np.inf)

    if contribution_values is not None:
        z = np.asarray(contribution_values, dtype=float)
        if z.shape != (q, label_count):
            raise ValueError("contribution_values must match vertices and labels")
        entries = [
            (i * label_count + j, float(z[i, j]))
            for i in range(q)
            for j in range(label_count)
        ]
        if threshold_sense == "le":
            append_row(entries, -np.inf, float(threshold))
        elif threshold_sense == "ge":
            append_row(entries, float(threshold), np.inf)
        else:
            raise ValueError("threshold_sense must be 'le' or 'ge'")

    matrix = coo_matrix(
        (values, (row_indices, column_indices)),
        shape=(len(row_lower), variable_count),
    ).tocsr()
    options = {"presolve": True, "mip_rel_gap": 0.0}
    if time_limit is not None:
        options["time_limit"] = float(time_limit)
    result = milp(
        objective,
        integrality=np.ones(variable_count),
        bounds=Bounds(lower_bounds, upper_bounds),
        constraints=LinearConstraint(matrix, row_lower, row_upper),
        options=options,
    )
    if result.x is None:
        return MilpCertificate(
            status=str(result.message),
            objective=None,
            dual_bound=float(result.mip_dual_bound) if getattr(result, "mip_dual_bound", None) is not None else None,
            relative_gap=float(result.mip_gap) if getattr(result, "mip_gap", None) is not None else None,
            node_count=int(result.mip_node_count) if getattr(result, "mip_node_count", None) is not None else None,
            labels=None,
        )
    # Solver success is not a semantic witness certificate. Reject invalid
    # assignments and threshold violations before exposing an incumbent.
    if not np.all(np.isfinite(result.x)) or np.any(np.abs(result.x-np.rint(result.x)) > 1e-6):
        raise RuntimeError("MILP returned a nonfinite or nonintegral witness")
    chosen_indices = tuple(
        int(np.argmax(result.x[i * label_count:(i + 1) * label_count]))
        for i in range(q)
    )
    chosen = tuple(labels[j] for j in chosen_indices)
    if any(label not in choices for label, choices in zip(chosen, allowed)):
        raise RuntimeError("MILP witness violates allowed label lists")
    assignments = np.rint(result.x[:x_count]).reshape(q, label_count)
    if np.any(assignments.sum(axis=1) != 1):
        raise RuntimeError("MILP witness violates one-label-per-vertex constraints")
    witness_cost = sum((rational(w) for i,j,w in checked_edges if chosen[i] != chosen[j]), rational(0))
    if not np.isclose(float(witness_cost), result.fun, rtol=1e-8, atol=1e-10):
        raise RuntimeError("MILP objective disagrees with recomputed cut cost")
    attained = None
    if contribution_values is not None:
        exact_attained = member_value(contribution_values, chosen_indices)
        exact_threshold = rational(threshold)
        valid = exact_attained <= exact_threshold if threshold_sense == "le" else exact_attained >= exact_threshold
        if not valid:
            raise RuntimeError("MILP witness fails the exact stored-input crossing inequality")
        attained = float(exact_attained)
    return MilpCertificate(
        status="optimal" if result.success else str(result.message),
        objective=float(result.fun),
        dual_bound=float(result.mip_dual_bound) if getattr(result, "mip_dual_bound", None) is not None else None,
        relative_gap=float(result.mip_gap) if getattr(result, "mip_gap", None) is not None else None,
        node_count=int(result.mip_node_count) if getattr(result, "mip_node_count", None) is not None else None,
        labels=chosen,
        attained_value=attained,
    )


def milp_directional_full_width(allowed_labels, edges, *, time_limit=None):
    """Minimum graph-relaxation cost attaining full directional width."""
    return _solve_potts(allowed_labels, edges, time_limit=time_limit)


def milp_decision_destroying_relaxation(
    contributions: np.ndarray,
    edges,
    threshold: float,
    *,
    time_limit=None,
):
    """Minimum cut cost that changes an initially decisive scalar state."""
    z = np.asarray(contributions, dtype=float)
    if z.ndim != 2 or z.shape[0] < 1 or z.shape[1] < 1 or not np.all(np.isfinite(z)):
        raise ValueError("contributions must be a finite nonempty (vertices, pathways) array")
    checked_edges = _edges(z.shape[0], edges)
    _require_connected(z.shape[0], checked_edges)
    lower, upper, independent_lower, independent_upper = scalar_endpoints(z)
    threshold_exact = rational(threshold)
    state = exact_state(lower, upper, threshold_exact)
    if state == "indeterminate":
        return state, MilpCertificate("initially_indeterminate", 0.0, 0.0, 0.0, 0, None)
    impossible = (independent_lower > threshold_exact if state == "resolved_above"
                  else independent_upper < threshold_exact)
    if impossible:
        return state, MilpCertificate("infeasible_exact_endpoint", None, None, None, 0, None)
    labels = tuple(range(z.shape[1]))
    allowed = tuple(labels for _ in range(z.shape[0]))
    sense = "le" if state == "resolved_above" else "ge"
    certificate = _solve_potts(
        allowed,
        checked_edges,
        contribution_values=z,
        threshold=threshold,
        threshold_sense=sense,
        time_limit=time_limit,
    )
    if certificate.labels is None and "infeasible" in certificate.status.lower():
        raise RuntimeError("MILP infeasibility contradicts the exact independent endpoint")
    return state, certificate


__all__ = [
    "MilpCertificate",
    "milp_decision_destroying_relaxation",
    "milp_directional_full_width",
]


