"""Diagnose coupling partitions attaining independent-set extrema within tolerance."""

from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
INPUT = HERE / "data" / "cfd_mapping"
RTOL = 1.0e-10


def canonical(blocks):
    blocks = [tuple(sorted(block)) for block in blocks]
    return tuple(sorted(blocks, key=lambda b: (b[0], len(b), b)))


def label_partition(labels):
    groups = {}
    for i, label in enumerate(labels):
        groups.setdefault(label, []).append(i)
    return canonical(groups.values())


def pair_partition(left, right):
    return label_partition(tuple(zip(left, right)))


def common_refinement(left, right):
    blocks = []
    for a in left:
        for b in right:
            intersection = sorted(set(a).intersection(b))
            if intersection:
                blocks.append(intersection)
    return canonical(blocks)


def signature(partition, fidelities):
    return " | ".join("+".join(fidelities[i] for i in block) for block in partition)


def minimal_by_refinement(partitions):
    """Keep coarsest witnesses; finer versions add no admissibility information."""
    unique = set(partitions)
    kept = []
    for candidate in unique:
        redundant = False
        for other in unique:
            if candidate == other:
                continue
            # candidate refines other if each candidate block lies in an other block.
            if all(any(set(block) <= set(parent) for parent in other) for block in candidate):
                redundant = True
                break
        if not redundant:
            kept.append(candidate)
    return sorted(kept, key=lambda pi: (len(pi), pi))


def evaluate(path):
    with np.load(path) as data:
        raw = np.asarray(data["weighted_cell_vectors"], dtype=np.float64)
        pathways = [str(x) for x in data["pathways"]]
        fidelities = [str(x) for x in data["fidelities"]]
        coefficients = np.asarray(data["coefficients"], dtype=np.float64)
        case_id = str(data["case_id"][0])
        support = str(data["support"][0])

    cell = {
        (p, f): raw[p_i * len(fidelities) + f_i]
        for p_i, p in enumerate(pathways)
        for f_i, f in enumerate(fidelities)
    }
    choices = list(itertools.product(pathways, repeat=len(fidelities)))
    vectors = np.stack([
        sum(coefficients[i] * cell[(choice[i], fidelities[i])]
            for i in range(len(fidelities)))
        for choice in choices
    ])
    gram = vectors @ vectors.T
    norm2 = np.maximum(np.diag(gram), 0.0)
    distance2 = norm2[:, None] + norm2[None, :] - 2.0 * gram
    np.maximum(distance2, 0.0, out=distance2)
    min_value = float(norm2.min())
    max_value = float(distance2.max())

    min_indices = np.flatnonzero(np.isclose(norm2, min_value, rtol=RTOL, atol=0.0))
    max_pairs = np.argwhere(np.isclose(distance2, max_value, rtol=RTOL, atol=0.0))
    min_witnesses = minimal_by_refinement(label_partition(choices[i]) for i in min_indices)
    diameter_witnesses = minimal_by_refinement(
        pair_partition(choices[i], choices[j]) for i, j in max_pairs if i < j
    )
    joint = minimal_by_refinement(
        common_refinement(a, b) for a in min_witnesses for b in diameter_witnesses
    )

    return {
        "case_id": case_id,
        "support": support,
        "minimum_norm": float(np.sqrt(min_value)),
        "diameter": float(np.sqrt(max_value)),
        "minimum_member_count_within_tolerance": int(len(min_indices)),
        "diameter_pair_count_within_tolerance": int(sum(i < j for i, j in max_pairs)),
        "coarsest_minimum_norm_witness_partitions": [
            signature(pi, fidelities) for pi in min_witnesses
        ],
        "coarsest_diameter_witness_partitions": [
            signature(pi, fidelities) for pi in diameter_witnesses
        ],
        "coarsest_joint_endpoint_witness_partitions": [
            signature(pi, fidelities) for pi in joint
        ],
        "minimum_blocks_for_independent_minimum": min(map(len, min_witnesses)),
        "minimum_blocks_for_independent_diameter": min(map(len, diameter_witnesses)),
        "minimum_blocks_for_both_independent_extrema": min(map(len, joint)),
    }


def main():
    results = [evaluate(path) for path in sorted(INPUT.glob("*.npz"))]
    payload = {
        "status": "completed_extremal_coupling_certificate",
        "relative_tolerance": RTOL,
        "result_count": len(results),
        "results": results,
        "certificate_logic": (
            "A realization is admissible under a partition exactly when every "
            "block is constant in pathway label. A diameter pair is admissible "
            "exactly when every block is constant in the ordered pair of its two "
            "pathway labels. Therefore any refinement of a listed witness "
            "partition attains the corresponding fully independent extremum within the declared relative tolerance."
        ),
    }
    out = HERE / "06_extremal_coupling_certificate.json"
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"wrote {out.name}")
    for item in results:
        print(
            item["case_id"], item["support"],
            "min_blocks=", item["minimum_blocks_for_independent_minimum"],
            "diam_blocks=", item["minimum_blocks_for_independent_diameter"],
            "joint_blocks=", item["minimum_blocks_for_both_independent_extrema"],
        )


if __name__ == "__main__":
    main()
