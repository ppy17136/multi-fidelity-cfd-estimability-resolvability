"""Certify diameter coupling complexity through directional set cover.

The script is intentionally read-only: it prints a JSON certificate to stdout.
For the released four-cell audits, the independent diameter is still found by
complete endpoint enumeration.  The new part is the exact reduction of each
diameter direction to a set-cover problem over at most P**2 pathway pairs.
"""

from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
INPUT = HERE / "data" / "cfd_mapping"
RTOL = 1.0e-10


def exact_cover(candidate_masks: dict[tuple[int, int], int], full_mask: int):
    items = sorted(candidate_masks)
    for size in range(1, len(items) + 1):
        for selected in itertools.combinations(items, size):
            covered = 0
            for item in selected:
                covered |= candidate_masks[item]
            if covered == full_mask:
                return selected
    raise RuntimeError("No pathway-pair cover exists")


def greedy_cover(candidate_masks: dict[tuple[int, int], int], full_mask: int):
    uncovered = full_mask
    selected = []
    while uncovered:
        best = max(
            candidate_masks,
            key=lambda item: ((candidate_masks[item] & uncovered).bit_count(), item),
        )
        gain = candidate_masks[best] & uncovered
        if not gain:
            raise RuntimeError("Greedy cover stalled")
        selected.append(best)
        uncovered &= ~candidate_masks[best]
    return tuple(selected)


def partition_from_cover(cover, candidate_masks, q):
    remaining = (1 << q) - 1
    blocks = []
    for pair in cover:
        assigned = candidate_masks[pair] & remaining
        if assigned:
            blocks.append([i for i in range(q) if assigned & (1 << i)])
            remaining &= ~assigned
    if remaining:
        raise RuntimeError("Cover-to-partition conversion left cells unassigned")
    return blocks


def evaluate(path: Path):
    with np.load(path) as data:
        raw = np.asarray(data["weighted_cell_vectors"], dtype=np.float64)
        pathways = [str(x) for x in data["pathways"]]
        fidelities = [str(x) for x in data["fidelities"]]
        coefficients = np.asarray(data["coefficients"], dtype=np.float64)
        case_id = str(data["case_id"][0])
        support = str(data["support"][0])

    p_count = len(pathways)
    q = len(fidelities)
    contributions = np.stack([
        np.stack([
            coefficients[i] * raw[p * q + i]
            for p in range(p_count)
        ])
        for i in range(q)
    ])

    choices = list(itertools.product(range(p_count), repeat=q))
    vectors = np.stack([
        sum(contributions[i, choice[i]] for i in range(q))
        for choice in choices
    ])
    gram = vectors @ vectors.T
    norm2 = np.maximum(np.diag(gram), 0.0)
    distance2 = norm2[:, None] + norm2[None, :] - 2.0 * gram
    np.maximum(distance2, 0.0, out=distance2)
    max_distance2 = float(distance2.max())
    diameter = float(np.sqrt(max_distance2))
    diameter_pairs = [
        (int(i), int(j))
        for i, j in np.argwhere(
            np.isclose(distance2, max_distance2, rtol=RTOL, atol=0.0)
        )
        if i < j
    ]

    full_mask = (1 << q) - 1
    direction_results = []
    for left, right in diameter_pairs:
        delta = vectors[left] - vectors[right]
        direction = delta / np.linalg.norm(delta)
        projections = np.einsum("ipd,d->ip", contributions, direction)
        maxima = projections.max(axis=1)
        minima = projections.min(axis=1)
        scale = max(1.0, float(np.max(np.abs(projections))))
        atol = 1.0e-12 * scale

        masks = {}
        for a in range(p_count):
            for b in range(p_count):
                mask = 0
                for i in range(q):
                    is_max = np.isclose(projections[i, a], maxima[i], rtol=RTOL, atol=atol)
                    is_min = np.isclose(projections[i, b], minima[i], rtol=RTOL, atol=atol)
                    if is_max and is_min:
                        mask |= 1 << i
                if mask:
                    masks[(a, b)] = mask

        exact = exact_cover(masks, full_mask)
        greedy = greedy_cover(masks, full_mask)
        blocks = partition_from_cover(exact, masks, q)
        direction_results.append({
            "left_choice": [pathways[p] for p in choices[left]],
            "right_choice": [pathways[p] for p in choices[right]],
            "exact_cover_size": len(exact),
            "exact_pathway_pair_cover": [
                [pathways[a], pathways[b]] for a, b in exact
            ],
            "greedy_cover_size": len(greedy),
            "greedy_pathway_pair_cover": [
                [pathways[a], pathways[b]] for a, b in greedy
            ],
            "certified_partition": [
                [fidelities[i] for i in block] for block in blocks
            ],
        })

    best_size = min(item["exact_cover_size"] for item in direction_results)
    best = [item for item in direction_results if item["exact_cover_size"] == best_size]
    return {
        "case_id": case_id,
        "support": support,
        "pathway_count": p_count,
        "cell_count": q,
        "independent_member_count": p_count ** q,
        "independent_diameter": diameter,
        "diameter_pair_count_within_tolerance": len(diameter_pairs),
        "minimum_directional_cover_size": best_size,
        "compressed_member_count": p_count ** best_size,
        "compression_factor": (p_count ** q) / (p_count ** best_size),
        "greedy_is_exact_for_all_diameter_directions": all(
            item["greedy_cover_size"] == item["exact_cover_size"]
            for item in direction_results
        ),
        "best_certificates": best,
    }


def main():
    results = [evaluate(path) for path in sorted(INPUT.glob("*.npz"))]
    payload = {
        "status": "directional_set_cover_certificate",
        "relative_tolerance": RTOL,
        "input_count": len(results),
        "all_minimum_cover_sizes_equal_two": all(
            item["minimum_directional_cover_size"] == 2 for item in results
        ),
        "results": results,
        "scope": (
            "Independent Euclidean diameter endpoints are enumerated for the "
            "released four-cell data.  Conditional on each certified diameter "
            "direction, minimum coupling-block count is solved exactly as set cover."
        ),
    }
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
