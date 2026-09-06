"""Frozen scaling benchmark for exact forest coupling relaxation.

The generator uses four pathways and therefore at most 16 ordered
maximizer/minimizer labels per fidelity cell. Most cells admit all labels;
deterministically spaced anchor cells admit one seeded label. This preserves
the edge-dominated O(|E| L^2) workload while preventing the optimum from being identically
zero. Chain and random-tree topologies are tested with three frozen seeds.
"""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import platform
import statistics
import time
import tracemalloc

import numpy as np

from graph_constrained_coupling import cut_cost, tree_graph_relaxation


HERE = Path(__file__).resolve().parent
CSV_OUTPUT = HERE / "32_FOREST_DP_SCALING_BENCHMARK.csv"
JSON_OUTPUT = HERE / "32_FOREST_DP_SCALING_BENCHMARK.json"
SIZES = (256, 1024, 4096, 16384)
SEEDS = (0, 1, 2)
LABELS = tuple((left, right) for left in range(4) for right in range(4))
ANCHOR_PERIOD = 17


def allowed_labels(q, rng):
    allowed = []
    for cell in range(q):
        if cell % ANCHOR_PERIOD == 0:
            allowed.append((LABELS[int(rng.integers(0, len(LABELS)))],))
        else:
            allowed.append(LABELS)
    return tuple(allowed)


def chain_edges(q, rng):
    return tuple(
        (cell - 1, cell, float(rng.integers(1, 11)))
        for cell in range(1, q)
    )


def random_tree_edges(q, rng):
    return tuple(
        (int(rng.integers(0, cell)), cell, float(rng.integers(1, 11)))
        for cell in range(1, q)
    )


def benchmark(topology, q, seed):
    rng = np.random.default_rng(202609040000 + q * 100 + seed)
    allowed = allowed_labels(q, rng)
    edges = chain_edges(q, rng) if topology == "chain" else random_tree_edges(q, rng)
    tracemalloc.start()
    start = time.perf_counter()
    certificate = tree_graph_relaxation(allowed, edges)
    runtime = time.perf_counter() - start
    _, peak_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    recomputed = cut_cost(certificate.labels, edges)
    feasible = all(label in labels for label, labels in zip(certificate.labels, allowed))
    return {
        "topology": topology,
        "cells": q,
        "edges": len(edges),
        "maximum_labels_per_cell": max(map(len, allowed)),
        "anchor_period": ANCHOR_PERIOD,
        "seed": seed,
        "objective": certificate.cost,
        "objective_recomputed": recomputed,
        "objective_check_abs": abs(certificate.cost - recomputed),
        "labels_feasible": feasible,
        "runtime_s": runtime,
        "peak_tracemalloc_mib": peak_bytes / (1024.0 * 1024.0),
    }


def scaling_slope(rows, topology):
    medians = []
    for q in SIZES:
        values = [row["runtime_s"] for row in rows
                  if row["topology"] == topology and row["cells"] == q]
        medians.append(statistics.median(values))
    return float(np.polyfit(np.log(np.asarray(SIZES)), np.log(np.asarray(medians)), 1)[0])


def main():
    # Warm imports and bytecode paths outside timed records.
    tree_graph_relaxation((LABELS, LABELS), ((0, 1, 1.0),))
    rows = [
        benchmark(topology, q, seed)
        for topology in ("chain", "random_tree")
        for q in SIZES
        for seed in SEEDS
    ]
    with CSV_OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        topology: {
            "empirical_log_log_runtime_slope": scaling_slope(rows, topology),
            "largest_size_median_runtime_s": statistics.median(
                row["runtime_s"] for row in rows
                if row["topology"] == topology and row["cells"] == max(SIZES)
            ),
            "largest_size_max_peak_tracemalloc_mib": max(
                row["peak_tracemalloc_mib"] for row in rows
                if row["topology"] == topology and row["cells"] == max(SIZES)
            ),
        }
        for topology in ("chain", "random_tree")
    }
    payload = {
        "status": "frozen_forest_dp_scaling_benchmark",
        "generator": {
            "sizes": list(SIZES),
            "seeds": list(SEEDS),
            "pathways": 4,
            "maximum_ordered_extremizer_labels": len(LABELS),
            "anchor_period": ANCHOR_PERIOD,
            "edge_weights": "seeded integers from 1 through 10",
        },
        "machine": {
            "platform": platform.platform(),
            "python": platform.python_version(),
            "processor": platform.processor(),
            "logical_cpu_count": os.cpu_count(),
            "numpy": np.__version__,
        },
        "all_certificates_verified": all(
            row["labels_feasible"] and row["objective_check_abs"] == 0.0
            for row in rows
        ),
        "summary": summary,
        "rows": rows,
        "claim_boundary": (
            "The measurements validate this reference implementation on the "
            "frozen forest families; they do not establish performance on general graphs."
        ),
    }
    JSON_OUTPUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({
        "case_count": len(rows),
        "all_certificates_verified": payload["all_certificates_verified"],
        "summary": summary,
    }, indent=2))


if __name__ == "__main__":
    main()
