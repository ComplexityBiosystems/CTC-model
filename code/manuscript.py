import glob

import pandas as pd
import numpy as np
from scipy.stats import linregress

from typing import List
from typing import Dict

organs = """brain large_intestine left_lung pancreas right_lung stomach heart left_kidney liver right_kidney small_intestine prostate gallbladder urinary_bladder left_adrenal_gland right_adrenal_gland""".split(
    " ")
coefs = [0.28150614,  1.02239547,  0.38434443, -0.10483313,  0.47272915,  0.20556007,
         1.3477815,   0.3658177,   0.42736813,  0.59098282,  0.42546357, 0.30061422,
         0.32449906]


def analyse_one_cohort(data, model):
    data_hq = data.loc[data._num_studies >= 2, ["cancer", "metastatic_site",
                                                "fraction_of_patients_with_metastasis", "patients_with_some_metastasis", "Study"]]

    _metada_data = []
    for cancer, metastatic_site, fraction_data, num_patients_data, study in data_hq.values:
        _x = model.loc[(model.cancer == cancer) & (
            model["metastatic site"] == metastatic_site)]

        assert len(_x) in [0, 1]

        if len(_x) == 1:
            fraction_model = _x.iloc[0]["value"]
            _metada_data.append([
                fraction_data,
                num_patients_data,
                fraction_model,
                cancer,
                metastatic_site,
                study
            ])
    metadata = pd.DataFrame(
        data=_metada_data,
        columns=["fraction_data", "num_patients_data",
                 "fraction_model", "cancer", "metastatic_site", "Study"]
    )
    assert len(metadata) > 4
    fraction_data, num_patients_data, fraction_model = metadata[[
        "fraction_data", "num_patients_data", "fraction_model"]].values.T
    num_patients_data = num_patients_data.astype(int)

    result = compute_explained_variance(
        fraction_data=fraction_data,
        num_patients_data=num_patients_data,
        fraction_model=fraction_model
    )

    return result, metadata


def compute_explained_variance(
    fraction_data: List[float],
    num_patients_data: List[int],
    fraction_model: List[float],
) -> Dict[str, float]:
    """Compute the explained variance of the model wrt the data.

    Arguments:
        fraction_data {List[float]} -- Fraction of patients in autopsy data
        num_patients_data {List[int]} -- Number of patients in autopsy data
        fraction_model {List[float]} -- Fraction of cells in model

    Returns:
        Dict[str, float] -- The various bits of the variance decomposition.
    """
    # check correctness of input data
    assert len(fraction_model) == len(fraction_data) == len(num_patients_data)
    for x, y, z in zip(fraction_model, fraction_data, num_patients_data):
        assert isinstance(x, float)
        assert isinstance(y, float)
        assert isinstance(z, (int, np.int64))
        assert x >= 0 and x <= 1
        assert z > 0

    # remove zeros
    _fraction_data, _fraction_model, _num_patients_data = np.array([
        [x, y, z]
        for x, y, b in zip(
            fraction_data,
            fraction_model,
            fraction_data * fraction_model > 0
        )
        if b
    ]).T

    # take logarithms
    logfraction_data = np.log(_fraction_data)
    logfraction_model = np.log(_fraction_model)

    # do linear regression
    x = logfraction_model
    y = logfraction_data
    slope, intercept, _, pvalue, _ = linregress(x, y)

    f = intercept + slope * x
    E = np.std(f - y)
    S = np.std(y)
    R2 = (S**2 - E**2) / S**2

    # estimate measurement error
    # error propagation formula
    # f(x) = log(x) --> sigma_f = |sigma_x / x|
    _N = _num_patients_data
    _p = _fraction_data
    _sigma = np.sqrt(_p * (1 - _p) / _N)
    measurement_errors = np.abs(_sigma / _p)
    e = np.sqrt(np.mean(measurement_errors**2))

    # adjusted pearson correlation coefficient
    r2 = (S**2 - E**2) / (S**2 - e**2)

    result = {
        "R2": R2,
        "R2_adjusted": r2,
        "pvalue": pvalue,
        "S2": S**2,
        "E2": E**2,
        "e2": e**2,
        "flow": S**2 - E**2,
        "other": E**2 - e**2,
        "noise": e**2,
        "meansquare": np.mean(y)**2
    }
    return result


def get_switch_probabilities(tracers_df):
    df = tracers_df
    contact_points_dict = {
        (fpath.split("/")[-1].split("_solved_")[0]
         ): set(pd.read_csv(fpath, header=None).values.T[0])
        for fpath in glob.glob("../output/data/organs_contact_points/*_solved_main_network_threshold10.0_contacts.txt")
    }

    _to_given_from = {
        organ: df.final_node.apply(
            lambda x: x in contact_points_dict[organ]).mean()
        for organ in contact_points_dict.keys()
    }
    S = sum(_to_given_from.values())

    to_given_from = pd.Series({k: v/S for k, v in _to_given_from.items()})
    to_given_from.sort_values(ascending=False, inplace=True)

    return to_given_from


def load_trajectories_sims():
    """Load trajectory simulations and create a primary-organ-to-metastatic-site
    probability matrix.

    Keyword Arguments:
        chunks {int} -- number of chunks to load (default: {10})

    Returns:
        pd.DataFrame -- Primary-organ-to-metastatic-site probability matrix
    """
    dfs = []
    for i in range(1):
        df = pd.DataFrame({
            organ_name: get_switch_probabilities(pd.read_pickle(
                f"../output/data/tracers_notrajs/{organ_name}_10K_tracers.p"))
            for organ_name in organs
        })
        df.columns.name = "Site of origin"
        df.index.name = "Metastatic site"
        df = df.T
        dfs.append(df)
    df = sum(dfs) / len(dfs)
    return df
