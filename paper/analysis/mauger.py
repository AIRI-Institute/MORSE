"""Reconstruct exact rank composites and fixed-score permutation statistics."""
from __future__ import annotations

import math
import numpy as np
from scipy.stats import rankdata

from ..io import DATA, read_json, rows, write_csv, write_json


def calculate(output=None, permutations=True):
    path = DATA / "SupplementaryFigure3_Mauger_U"
    observations = rows(path / "measurements_and_scores.csv")
    protocol = read_json(path / "analysis_protocol.json")
    display = read_json(path / "plot_data.json")
    components = protocol["individual_scorers_retained"]
    n = len(observations)
    if n != 39:
        raise ValueError("Mauger analysis requires exactly 39 constructs")
    y = np.asarray([float(row["U_bar_height_pt"]) for row in observations])
    ranks = {key: rankdata([float(row[key]) for row in observations], method="average") for key in components}
    for label, members in protocol["combinations"].items():
        products = [math.prod(int(round(2*ranks[key][i])) for key in members) for i in range(n)]
        # Python integers preserve exact products before tied ranks are formed.
        ranks[label] = np.asarray([1+sum(other < value for other in products)+(products.count(value)-1)/2 for value in products])
    keys = components + ["all translation", "All Relevant"]
    yr = rankdata(y, method="average")
    def unit(v):
        z = v-v.mean()
        norm=np.linalg.norm(z)
        if not np.isfinite(z).all() or not math.isfinite(float(norm)) or norm<=0:
            raise ValueError("Undefined correlation: nonfinite or constant rank input")
        return z/norm
    observed = {key: float(unit(yr) @ unit(ranks[key])) for key in keys}
    order = np.argsort(y, kind="stable")
    matrix = np.asarray([100*(yr[order]-1)/(n-1)] + [100*(ranks[key][order]-1)/(n-1) for key in keys])
    np.testing.assert_allclose(matrix, display["percentiles"], rtol=0, atol=1e-12)
    pvalues = {}
    if permutations:
        rng = np.random.default_rng(20261129)
        permuted = np.empty((100000, n))
        permuted[0] = unit(yr)
        for start in range(1, len(permuted), 4096):
            count = min(4096, len(permuted)-start)
            permuted[start:start+count] = rng.permuted(np.broadcast_to(permuted[0], (count,n)), axis=1)
        if not np.isfinite(permuted).all() or np.max(np.abs(permuted))>1:
            raise ValueError("Invalid bounded permutation rank vectors")
        for key in keys:
            # Small bounded 39-term reductions use an explicit non-BLAS path.
            # Some Apple BLAS builds leave spurious floating status flags on
            # matmul; finite/range checks still reject actual invalid statistics.
            null = np.einsum("ij,j->i",permuted,unit(ranks[key]),optimize=False)
            if not np.isfinite(null).all() or np.max(np.abs(null))>1+1e-12:
                raise ValueError(f"Invalid permutation correlations for {key}")
            pvalues[key] = float(np.mean(np.round(np.abs(null),12) >= round(abs(observed[key]),12)))
    expected = {row["scorer"]: row for row in display["statistics"]}
    result = []
    for key in keys:
        if not math.isclose(observed[key], expected[key]["spearman_rho"], abs_tol=1e-12):
            raise ValueError(f"Reconstructed Spearman statistic mismatch: {key}")
        if permutations and pvalues[key] != expected[key]["p_value"]:
            raise ValueError(f"Reconstructed permutation p-value mismatch: {key}")
        result.append({"scorer":key,"n":n,"spearman_rho":observed[key],"p_value":pvalues.get(key),
                       "random_permutations":99999 if permutations else 0,"permutation_seed":20261129})
    if output is not None:
        write_csv(output / "Figure2_recomputed_statistics.csv",result)
        write_json(output / "Figure2_recomputed_rank_heatmap.json",{"ids":[observations[i]["id"] for i in order],"rows":["experiment",*keys],"percentiles":matrix.tolist()})
    return {"n_constructs":n,"score_rows":len(keys),"rank_heatmap_exact":True,"all_correlations_match":True,
            "all_permutation_pvalues_match":True if permutations else None,"permutations":99999 if permutations else 0}
