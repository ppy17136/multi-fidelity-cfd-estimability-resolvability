"""Graph-constrained numerical-pathway coupling and decision frontiers.

Vertices represent fidelity cells.  An edge with positive weight records the
cost of allowing its endpoint cells to use different numerical pathways.
Consequently, a pathway labelling induces a cut, and the connected components
after cutting are the independently selectable pathway blocks.

The routines here are exact reference implementations.  Brute force is used
for general small graphs, while directional full-width attainment on a forest
is solved by dynamic programming.
"""

from __future__ import annotations

from dataclasses import dataclass
import itertools
from typing import Iterable, Sequence

import numpy as np
import networkx as nx


Label = tuple[int, int]
Edge = tuple[int, int, float]


@dataclass(frozen=True)
class LabelingCertificate:
    cost: float
    labels: tuple


@dataclass(frozen=True)
class DecisionCertificate:
    coherent_state: str
    threshold: float
    minimum_relaxation_cost: float | None
    labels: tuple[int, ...] | None
    attained_value: float | None
    lower_frontier: tuple[tuple[float, float], ...]
    upper_frontier: tuple[tuple[float, float], ...]


def _normalise_edges(q: int, edges: Iterable[Sequence[float]]) -> tuple[Edge, ...]:
    result = []
    seen = set()
    for raw in edges:
        if len(raw) != 3:
            raise ValueError("each edge must be (left, right, nonnegative_weight)")
        left, right, weight = int(raw[0]), int(raw[1]), float(raw[2])
        if not (0 <= left < q and 0 <= right < q) or left == right:
            raise ValueError("edge endpoints must be distinct valid vertex indices")
        if weight < 0 or not np.isfinite(weight):
            raise ValueError("edge weights must be finite and nonnegative")
        key = tuple(sorted((left, right)))
        if key in seen:
            raise ValueError("parallel edges are not supported by this reference code")
        seen.add(key)
        result.append((key[0], key[1], weight))
    return tuple(result)


def _require_connected(q: int, edges: Sequence[Edge]) -> None:
    """Reject ambiguous zero-cost baselines for decision-frontier routines."""
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


def cut_cost(labels: Sequence, edges: Iterable[Sequence[float]]) -> float:
    labels = tuple(labels)
    checked = _normalise_edges(len(labels), edges)
    return float(sum(weight for left, right, weight in checked if labels[left] != labels[right]))


def pathway_pair_compatibility(
    projections: np.ndarray,
    *,
    relative_tolerance: float = 0.0,
    absolute_tolerance: float = 0.0,
) -> tuple[tuple[Label, ...], ...]:
    """Return ordered maximizing--minimizing pathway-pair labels.

    The exact graph-equivalence theorem requires exact extrema, so both
    tolerances default to zero. Positive tolerances deliberately enlarge the
    label lists to treat near-extrema as ties; a certificate obtained from such
    enlarged lists applies to that declared numerical-tie model and must not be
    called exact full-width attainment without recomputing the attained width.
    """
    z = np.asarray(projections, dtype=float)
    if z.ndim != 2 or z.shape[0] < 1 or z.shape[1] < 1:
        raise ValueError("projections must have nonempty shape (cells, pathways)")
    if not np.all(np.isfinite(z)):
        raise ValueError("projections must be finite")
    rtol = float(relative_tolerance)
    atol = float(absolute_tolerance)
    if rtol < 0 or atol < 0 or not np.isfinite(rtol) or not np.isfinite(atol):
        raise ValueError("extremizer tolerances must be finite and nonnegative")
    maxima = z.max(axis=1)
    minima = z.min(axis=1)
    allowed = []
    for i in range(z.shape[0]):
        maximizers = np.flatnonzero(
            np.isclose(z[i], maxima[i], rtol=rtol, atol=atol)
        )
        minimizers = np.flatnonzero(
            np.isclose(z[i], minima[i], rtol=rtol, atol=atol)
        )
        allowed.append(tuple((int(a), int(b)) for a in maximizers for b in minimizers))
    return tuple(allowed)


def exact_graph_relaxation_bruteforce(
    allowed_labels: Sequence[Sequence[Label]],
    edges: Iterable[Sequence[float]],
) -> LabelingCertificate:
    """Minimum cut cost whose components attain full independent width."""
    allowed = tuple(tuple(labels) for labels in allowed_labels)
    if not allowed or any(not labels for labels in allowed):
        raise ValueError("every vertex must have at least one allowed label")
    checked = _normalise_edges(len(allowed), edges)
    best_cost = np.inf
    best_labels = None
    for labels in itertools.product(*allowed):
        cost = cut_cost(labels, checked)
        if cost < best_cost:
            best_cost = cost
            best_labels = labels
    return LabelingCertificate(float(best_cost), tuple(best_labels))


def two_label_graph_relaxation(
    allowed_labels: Sequence[Sequence[Label]],
    edges: Iterable[Sequence[float]],
) -> LabelingCertificate:
    """Solve a one- or two-label general-graph instance by minimum cut.

    Vertices allowing only the first or second global label are attached to
    the corresponding terminal with a capacity larger than the total finite
    edge weight. Vertices allowing both labels remain free.
    """
    allowed = tuple(tuple(dict.fromkeys(labels)) for labels in allowed_labels)
    q = len(allowed)
    if not allowed or any(not labels for labels in allowed):
        raise ValueError("every vertex must have at least one allowed label")
    checked = _normalise_edges(q, edges)
    global_labels = tuple(sorted(set().union(*map(set, allowed))))
    if len(global_labels) > 2:
        raise ValueError("two_label_graph_relaxation supports at most two labels")
    if len(global_labels) == 1:
        return LabelingCertificate(0.0, (global_labels[0],) * q)

    first, second = global_labels
    forced_first = [i for i, choices in enumerate(allowed) if second not in choices]
    forced_second = [i for i, choices in enumerate(allowed) if first not in choices]
    if not forced_first or not forced_second:
        common = first if not forced_second else second
        return LabelingCertificate(0.0, (common,) * q)

    source = "__source__"
    sink = "__sink__"
    graph = nx.Graph()
    graph.add_nodes_from(range(q))
    total_weight = float(sum(weight for _, _, weight in checked))
    forcing_capacity = total_weight + 1.0
    for left, right, weight in checked:
        graph.add_edge(left, right, capacity=weight)
    for node in forced_first:
        graph.add_edge(source, node, capacity=forcing_capacity)
    for node in forced_second:
        graph.add_edge(node, sink, capacity=forcing_capacity)
    value, partition = nx.minimum_cut(graph, source, sink, capacity="capacity")
    source_side, sink_side = partition
    labels = tuple(first if i in source_side else second for i in range(q))
    return LabelingCertificate(float(value), labels)


def tree_graph_relaxation(
    allowed_labels: Sequence[Sequence[Label]],
    edges: Iterable[Sequence[float]],
) -> LabelingCertificate:
    """Exact forest dynamic program for minimum full-width relaxation cost.

    For supplied label lists, complexity is
    O(sum_i |L_i| + sum_{(i,j) in E} |L_i||L_j|), hence
    O(|V| L + |E| L^2) when each list has at most L labels.
    Disconnected forests, including isolated vertices, are handled componentwise.
    """
    allowed = tuple(tuple(labels) for labels in allowed_labels)
    q = len(allowed)
    if not allowed or any(not labels for labels in allowed):
        raise ValueError("every vertex must have at least one allowed label")
    checked = _normalise_edges(q, edges)
    adjacency = [[] for _ in range(q)]
    for left, right, weight in checked:
        adjacency[left].append((right, weight))
        adjacency[right].append((left, weight))

    # Build a rooted representation iteratively. Avoiding Python recursion is
    # important because a valid fidelity tree may contain thousands of cells.
    unseen = -2
    parent = [unseen] * q
    parent_weight = [0.0] * q
    roots = []
    order = []
    for root in range(q):
        if parent[root] != unseen:
            continue
        roots.append(root)
        parent[root] = -1
        stack = [root]
        while stack:
            node = stack.pop()
            order.append(node)
            for neighbour, weight in adjacency[node]:
                if neighbour == parent[node]:
                    continue
                if parent[neighbour] != unseen:
                    raise ValueError("tree_graph_relaxation requires a forest")
                parent[neighbour] = node
                parent_weight[neighbour] = weight
                stack.append(neighbour)

    children = [[] for _ in range(q)]
    for node in range(q):
        if parent[node] >= 0:
            children[parent[node]].append(node)

    dp = [None] * q
    decisions = [None] * q
    for node in reversed(order):
        node_dp = {}
        node_decisions = {}
        for label in allowed[node]:
            value = 0.0
            local = {}
            for child in children[node]:
                weight = parent_weight[child]
                child_label, child_value = min(
                    dp[child].items(),
                    key=lambda item: (
                        item[1] + (0.0 if item[0] == label else weight),
                        item[0],
                    ),
                )
                value += child_value + (0.0 if child_label == label else weight)
                local[child] = child_label
            node_dp[label] = value
            node_decisions[label] = local
        dp[node] = node_dp
        decisions[node] = node_decisions

    chosen = [None] * q
    total_cost = 0.0
    for root in roots:
        root_label, root_cost = min(dp[root].items(), key=lambda item: (item[1], item[0]))
        total_cost += root_cost
        stack = [(root, root_label)]
        while stack:
            node, label = stack.pop()
            chosen[node] = label
            stack.extend(
                (child, decisions[node][label][child])
                for child in children[node]
            )

    return LabelingCertificate(float(total_cost), tuple(chosen))


def coherent_scalar_range(contributions: np.ndarray) -> tuple[float, float]:
    z = np.asarray(contributions, dtype=float)
    if z.ndim != 2 or z.shape[0] < 1 or z.shape[1] < 1 or not np.all(np.isfinite(z)):
        raise ValueError("contributions must be a finite nonempty (cells, pathways) array")
    values = z.sum(axis=0)
    return float(values.min()), float(values.max())


def interval_state(lower: float, upper: float, threshold: float) -> str:
    if lower > threshold:
        return "resolved_above"
    if upper < threshold:
        return "resolved_below"
    return "indeterminate"


def _pareto_envelope(points, *, minimize_value: bool):
    # For equal costs retain the best attainable extremum, then retain only
    # strict improvements as relaxation cost increases.
    by_cost = {}
    for cost, value in points:
        if cost not in by_cost:
            by_cost[cost] = value
        elif minimize_value:
            by_cost[cost] = min(by_cost[cost], value)
        else:
            by_cost[cost] = max(by_cost[cost], value)
    frontier = []
    incumbent = np.inf if minimize_value else -np.inf
    for cost in sorted(by_cost):
        value = by_cost[cost]
        improves = value < incumbent if minimize_value else value > incumbent
        if improves:
            frontier.append((float(cost), float(value)))
            incumbent = value
    return tuple(frontier)


def exact_decision_relaxation_frontier(
    contributions: np.ndarray,
    edges: Iterable[Sequence[float]],
    threshold: float,
) -> DecisionCertificate:
    """Enumerate the exact cost/extremum frontier for a scalar decision.

    A common pathway is used at zero cut cost.  Different labels at the two
    ends of an edge pay its weight.  The certificate reports the least cost at
    which an initially decisive coherent state can become indeterminate.
    """
    z = np.asarray(contributions, dtype=float)
    if z.ndim != 2 or z.shape[0] < 1 or z.shape[1] < 1 or not np.all(np.isfinite(z)):
        raise ValueError("contributions must be a finite nonempty (cells, pathways) array")
    q, pathway_count = z.shape
    checked = _normalise_edges(q, edges)
    _require_connected(q, checked)
    coherent_lower, coherent_upper = coherent_scalar_range(z)
    state = interval_state(coherent_lower, coherent_upper, float(threshold))

    records = []
    for labels in itertools.product(range(pathway_count), repeat=q):
        value = float(sum(z[i, labels[i]] for i in range(q)))
        records.append((cut_cost(labels, checked), value, labels))

    lower_frontier = _pareto_envelope(
        ((cost, value) for cost, value, _ in records), minimize_value=True
    )
    upper_frontier = _pareto_envelope(
        ((cost, value) for cost, value, _ in records), minimize_value=False
    )

    feasible = []
    if state == "resolved_above":
        feasible = [record for record in records if record[1] <= threshold]
    elif state == "resolved_below":
        feasible = [record for record in records if record[1] >= threshold]

    if feasible:
        cost, value, labels = min(feasible, key=lambda item: (item[0], item[1], item[2]))
        return DecisionCertificate(
            state,
            float(threshold),
            float(cost),
            tuple(labels),
            float(value),
            lower_frontier,
            upper_frontier,
        )
    return DecisionCertificate(
        state,
        float(threshold),
        None,
        None,
        None,
        lower_frontier,
        upper_frontier,
    )


__all__ = [
    "DecisionCertificate",
    "LabelingCertificate",
    "coherent_scalar_range",
    "cut_cost",
    "exact_decision_relaxation_frontier",
    "exact_graph_relaxation_bruteforce",
    "interval_state",
    "pathway_pair_compatibility",
    "tree_graph_relaxation",
    "two_label_graph_relaxation",
]



