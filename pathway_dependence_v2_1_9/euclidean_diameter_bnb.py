"""Certified branch-and-bound for a finite Minkowski-sum diameter.

For each cell i, choose one vector d_i from a finite difference set D_i and
maximize ||sum_i d_i||_2.  This equals the Euclidean diameter of the fully
recombined pathway set when D_i contains all ordered pathway differences.

The upper bound is obtained by maximizing every remaining linear, quadratic,
and pairwise cross term independently.  It is conservative but rigorous.
"""

from __future__ import annotations

from dataclasses import dataclass
import itertools

import numpy as np


@dataclass(frozen=True)
class DiameterCertificate:
    diameter: float
    squared_diameter: float
    ordered_pair_indices: tuple[tuple[int, int], ...]
    nodes_visited: int
    leaves_evaluated: int
    nodes_pruned: int
    nominal_leaf_count: int
    intrinsic_rank: int


def intrinsic_pathway_coordinates(
    contributions: np.ndarray,
    relative_tolerance: float = 1.0e-12,
) -> tuple[np.ndarray, np.ndarray]:
    """Map q x P x ambient pathway contributions to exact Gram coordinates.

    Pathway zero is the affine reference in every cell.  Up to the supplied
    numerical rank tolerance, all within-cell differences are represented in
    at most q(P-1) dimensions.
    """
    contributions = np.asarray(contributions, dtype=float)
    if contributions.ndim != 3:
        raise ValueError("contributions must have shape (cells, pathways, dimension)")
    q, p_count, _ = contributions.shape
    generators = np.stack([
        contributions[i, p] - contributions[i, 0]
        for i in range(q)
        for p in range(1, p_count)
    ])
    gram = generators @ generators.T
    gram = (gram + gram.T) / 2.0
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    scale = max(float(eigenvalues[-1]), 0.0) if eigenvalues.size else 0.0
    keep = eigenvalues > relative_tolerance * scale
    kept_values = np.maximum(eigenvalues[keep], 0.0)
    generator_coordinates = eigenvectors[:, keep] * np.sqrt(kept_values)

    coordinates = np.zeros((q, p_count, int(np.count_nonzero(keep))), dtype=float)
    for i in range(q):
        for p in range(1, p_count):
            coordinates[i, p] = generator_coordinates[i * (p_count - 1) + p - 1]
    return coordinates, eigenvalues


def ordered_difference_options(
    pathway_coordinates: np.ndarray,
) -> tuple[list[np.ndarray], list[list[tuple[int, int]]]]:
    q, p_count, _ = pathway_coordinates.shape
    options = []
    labels = []
    for i in range(q):
        cell_labels = list(itertools.product(range(p_count), repeat=2))
        cell_options = np.stack([
            pathway_coordinates[i, left] - pathway_coordinates[i, right]
            for left, right in cell_labels
        ])
        options.append(cell_options)
        labels.append(cell_labels)
    return options, labels


def brute_force_diameter(options: list[np.ndarray]):
    best_squared = -np.inf
    best_choice = None
    for choice in itertools.product(*[range(len(cell)) for cell in options]):
        total = sum(options[i][selected] for i, selected in enumerate(choice))
        squared = float(total @ total)
        if squared > best_squared:
            best_squared = squared
            best_choice = choice
    return best_squared, best_choice


def _coordinate_ascent(options: list[np.ndarray], starts: int = 32):
    q = len(options)
    seed_choices = [tuple(0 for _ in range(q))]
    rng = np.random.default_rng(20260904)
    for _ in range(max(0, starts - 1)):
        seed_choices.append(tuple(rng.integers(len(options[i])) for i in range(q)))

    best_squared = -np.inf
    best_choice = None
    for initial in seed_choices:
        choice = list(initial)
        total = sum(options[i][choice[i]] for i in range(q))
        for _ in range(20):
            changed = False
            for i in range(q):
                without = total - options[i][choice[i]]
                scores = np.einsum("ij,ij->i", without + options[i], without + options[i])
                selected = int(np.argmax(scores))
                if selected != choice[i]:
                    total = without + options[i][selected]
                    choice[i] = selected
                    changed = True
            if not changed:
                break
        squared = float(total @ total)
        if squared > best_squared:
            best_squared = squared
            best_choice = tuple(choice)
    return best_squared, best_choice


def certified_diameter_bnb(
    options: list[np.ndarray],
    relative_tolerance: float = 1.0e-12,
) -> DiameterCertificate:
    if not options:
        return DiameterCertificate(0.0, 0.0, (), 1, 1, 0, 1, 0)
    intrinsic_rank = int(options[0].shape[1])
    q = len(options)

    # Put the largest/highest-interaction cells first to improve pruning.
    radii = np.array([
        np.sqrt(np.max(np.einsum("ij,ij->i", cell, cell))) for cell in options
    ])
    interaction = np.zeros((q, q), dtype=float)
    for i in range(q):
        for j in range(i + 1, q):
            interaction[i, j] = interaction[j, i] = float(
                np.max(options[i] @ options[j].T)
            )
    order = tuple(np.argsort(-(radii + interaction.sum(axis=1))))
    ordered = [options[i] for i in order]

    incumbent, incumbent_choice_ordered = _coordinate_ascent(ordered)
    nodes_visited = 0
    leaves_evaluated = 0
    nodes_pruned = 0

    norm2_max = [np.max(np.einsum("ij,ij->i", cell, cell)) for cell in ordered]
    cross_max = np.zeros((q, q), dtype=float)
    for i in range(q):
        for j in range(i + 1, q):
            cross_max[i, j] = float(np.max(ordered[i] @ ordered[j].T))

    selected = [0] * q

    def upper_bound_squared(depth: int, partial: np.ndarray) -> float:
        value = float(partial @ partial)
        for i in range(depth, q):
            value += 2.0 * float(np.max(ordered[i] @ partial))
            value += float(norm2_max[i])
        for i in range(depth, q):
            for j in range(i + 1, q):
                value += 2.0 * cross_max[i, j]
        return max(value, 0.0)

    def search(depth: int, partial: np.ndarray):
        nonlocal incumbent, incumbent_choice_ordered
        nonlocal nodes_visited, leaves_evaluated, nodes_pruned
        nodes_visited += 1
        upper = upper_bound_squared(depth, partial)
        if upper <= incumbent * (1.0 + relative_tolerance):
            nodes_pruned += 1
            return
        if depth == q:
            leaves_evaluated += 1
            squared = float(partial @ partial)
            if squared > incumbent:
                incumbent = squared
                incumbent_choice_ordered = tuple(selected)
            return

        candidates = ordered[depth]
        scores = np.einsum("ij,ij->i", partial + candidates, partial + candidates)
        for option_index in np.argsort(-scores):
            selected[depth] = int(option_index)
            search(depth + 1, partial + candidates[option_index])

    search(0, np.zeros(intrinsic_rank, dtype=float))

    original_choice = [0] * q
    for ordered_position, original_cell in enumerate(order):
        original_choice[original_cell] = incumbent_choice_ordered[ordered_position]

    # Option index a*P+b is decoded by the caller; retain indices generically.
    pair_count = int(round(np.sqrt(len(options[0]))))
    decoded = tuple(divmod(index, pair_count) for index in original_choice)
    return DiameterCertificate(
        diameter=float(np.sqrt(max(incumbent, 0.0))),
        squared_diameter=float(max(incumbent, 0.0)),
        ordered_pair_indices=decoded,
        nodes_visited=nodes_visited,
        leaves_evaluated=leaves_evaluated,
        nodes_pruned=nodes_pruned,
        nominal_leaf_count=int(np.prod([len(cell) for cell in options], dtype=object)),
        intrinsic_rank=intrinsic_rank,
    )

