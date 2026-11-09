"""Inference metadata extraction from the upstream run table. See LICENSE."""
import json
import re
from typing import Any, Dict
import pandas as pd


def extract_config(run_df: pd.DataFrame, run_id: str) -> Dict[str, Any]:
    """Extract the config dict from mlflow run_df filtered by run_id.

    Args:
        run_df (pd.DataFrame): mlflow run_df
        run_id (str): mlflow run_id

    Returns:
        Dict: the configuration associated with the frozen model run_id
    """
    run_info = (
        run_df.query("run_id == @run_id")
        .filter(regex=r"^params\.", axis=1)
        .reset_index(drop=True)
    )
    run_info_dict = run_info.loc[0].to_dict()
    new_dict = dict()
    for k, v in run_info_dict.items():
        k = k.replace("params.", "")
        if v == "None":
            v = None
        elif v == "True":
            v = True
        elif v == "False":
            v = False
        elif re.search(
            r"(^\d+$)|(^[\d.]+$)|(^[\d.]+e-?[\d.]+$)|(^\[[\d.,e\- ]+\]$)", str(v)
        ):
            v = json.loads(str(v))
        if k != "betas":
            new_dict[k] = v
        else:
            match = re.match(r"\(([\d.]+),\s*([\d.]+)\)", v)
            if match:
                new_dict['adam_beta1'] = float(match.group(1))
                new_dict['adam_beta2'] = float(match.group(2))

    return new_dict
