#!/usr/bin/env python3
"""Derive metabolic, kidney, and cardiovascular baseline variables."""

import csv
from collections import Counter
import math
from pathlib import Path
import re
from tempfile import NamedTemporaryFile


BASELINE = Path("UKB_baseline.csv")

TG_CKM_CUTOFF = 135.0
BMI_CKM_CUTOFF = 25.0

# Apply the same waist cutoffs to all participants.
WAIST_FEMALE_CM = 88.0
WAIST_MALE_CM = 102.0

# Out-of-range PREVENT predictors remain missing.
PREVENT_RANGE_POLICY = "strict"

# Remove the BMI aliases created by the previous version.
DROP_COLUMNS = {
    "overweight_or_obesity",
    "obesity_bmi_ge30",
}

# Risk codes: 1=low, 2=moderate, 3=high, 4=very high.
# Each tuple corresponds to A1, A2, and A3.
KDIGO_RISK = {
    "G1": (1, 2, 3),
    "G2": (1, 2, 3),
    "G3a": (2, 3, 4),
    "G3b": (3, 4, 4),
    "G4": (4, 4, 4),
    "G5": (4, 4, 4),
}

REQUIRED = {
    "eid",
    "sex",
    "age_years",
    "bmi_kg_m2",
    "waist_cm",
    "hdl_mg_dl",
    "tc_mg_dl",
    "tg_mg_dl",
    "sbp_mmhg",
    "dbp_mmhg",
    "hba1c_pct",
    "serum_creatinine_mg_dl",
    "cystatin_c_mg_l",
    "urine_albumin_mg_l",
    "urine_creatinine_umol_l",
    "antihypertensive_medication_use",
    "lipid_lowering_medication_use",
    "diabetes_status",
    "current_smoking",
}

# Sex: 0=female, 1=male.
# PREVENT 10-year total-CVD base-model coefficients.
# Term order:
# intercept, A, N, H, SL, SH, D, C, GL, GH, B, T,
# B*SH, T*N, A*N, A*H, A*SH, A*D, A*C, A*GL.
PREVENT_COEFFICIENTS = {
    0: (
        -3.307728,
        0.7939329,
        0.0305239,
        -0.1606857,
        -0.2394003,
        0.360078,
        0.8667604,
        0.5360739,
        0.6045917,
        0.0433769,
        0.3151672,
        -0.1477655,
        -0.0663612,
        0.1197879,
        -0.0819715,
        0.0306769,
        -0.0946348,
        -0.27057,
        -0.078715,
        -0.1637806,
    ),
    1: (
        -3.031168,
        0.7688528,
        0.0736174,
        -0.0954431,
        -0.4347345,
        0.3362658,
        0.7692857,
        0.4386871,
        0.5378979,
        0.0164827,
        0.288879,
        -0.1337349,
        -0.0475924,
        0.150273,
        -0.0517874,
        0.0191169,
        -0.1049477,
        -0.2251948,
        -0.0895067,
        -0.1543702,
    ),
}


def number(raw):
    """Parse a finite number."""
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None

    return value if math.isfinite(value) else None


def positive(raw):
    """Parse a positive measurement without adding assay limits."""
    value = number(raw)
    return value if value is not None and value > 0 else None


def nonnegative(raw):
    """Parse a measurement that may legitimately equal zero."""
    value = number(raw)
    return value if value is not None and value >= 0 else None


def binary(raw):
    value = number(raw)
    return int(value) if value in (0, 1) else None


def threshold(value, cutoff):
    return None if value is None else int(value >= cutoff)


def any_positive(*flags):
    return 1 if 1 in flags else (None if None in flags else 0)


def classify_mets(
    waist_flag,
    sex,
    hdl,
    tg,
    sbp,
    dbp,
    bp_med,
    hba1c,
):
    """Classify the HbA1c-adapted MetS study phenotype."""
    bp = any_positive(
        threshold(sbp, 130),
        threshold(dbp, 80),
        bp_med,
    )

    low_hdl = (
        None
        if sex is None or hdl is None
        else int(hdl < (50 if sex == 0 else 40))
    )

    components = (
        waist_flag,
        low_hdl,
        threshold(tg, 150),
        bp,
        threshold(hba1c, 5.7),
    )

    positives = components.count(1)

    if positives >= 3:
        return 1

    if positives + components.count(None) < 3:
        return 0

    return None


def calculate_egfr(age, sex, scr, scys):
    """Calculate eGFR using the original equation hierarchy."""
    # Age: years; Scr: mg/dL; Scys: mg/L.
    # Output: mL/min/1.73 m2.
    # Legacy race multipliers are fixed at 1 for all participants.
    if age is None or age < 18 or sex not in (0, 1):
        return None, None

    scr = scr if scr is not None and scr > 0 else None
    scys = scys if scys is not None and scys > 0 else None

    female = sex == 0
    kappa = 0.7 if female else 0.9

    # Prefer the combined equation when both biomarkers are available.
    if scr is not None and scys is not None:
        x = scr / kappa
        y = scys / 0.8
        alpha = -0.248 if female else -0.207

        egfr = (
            135
            * min(x, 1) ** alpha
            * max(x, 1) ** -0.601
            * min(y, 1) ** -0.375
            * max(y, 1) ** -0.711
            * 0.995 ** age
            * (0.969 if female else 1)
        )

        return (
            egfr,
            "CKD-EPI_2012_cr_cys_race_multiplier_1",
        )

    if scr is not None:
        x = scr / kappa
        alpha = -0.329 if female else -0.411

        egfr = (
            141
            * min(x, 1) ** alpha
            * max(x, 1) ** -1.209
            * 0.993 ** age
            * (1.018 if female else 1)
        )

        return (
            egfr,
            "CKD-EPI_2009_cr_race_multiplier_1",
        )

    if scys is not None:
        y = scys / 0.8

        egfr = (
            133
            * min(y, 1) ** -0.499
            * max(y, 1) ** -1.328
            * 0.996 ** age
            * (0.932 if female else 1)
        )

        return egfr, "CKD-EPI_2012_cystatin_C"

    return None, None


def albumin_category(value):
    # KDIGO: A1 <30, A2 30-300, A3 >300 mg/g.
    if value < 30:
        return "A1"

    if value <= 300:
        return "A2"

    return "A3"


def calculate_uacr(row):
    """Use upstream substituted values directly as analysis point values."""
    # Albumin: mg/L.
    # Urine creatinine: umol/L.
    # UACR: mg/g.
    albumin = nonnegative(row.get("urine_albumin_mg_l"))
    creatinine = positive(row.get("urine_creatinine_umol_l"))

    ratio = None

    if albumin is not None and creatinine is not None:
        value = albumin / creatinine * 8840.0

        if math.isfinite(value):
            ratio = value

    return {
        "uacr_mg_g": ratio,
        "uacr_category": (
            None if ratio is None else albumin_category(ratio)
        ),
    }


def kidney_risk(egfr, a_category):
    """Assign KDIGO G category and risk code."""
    g = None if egfr is None else next(
        (
            name
            for lower, name in (
                (90, "G1"),
                (60, "G2"),
                (45, "G3a"),
                (30, "G3b"),
                (15, "G4"),
                (0, "G5"),
            )
            if egfr >= lower
        ),
        None,
    )

    if g is None:
        risk = None

    elif a_category is not None:
        risk = KDIGO_RISK[g][int(a_category[1]) - 1]

    else:
        # G4/G5 have very high risk regardless of the albumin category.
        risk = (
            KDIGO_RISK[g][0]
            if len(set(KDIGO_RISK[g])) == 1
            else None
        )

    return g, risk


def statin_input(row):
    """Use the lipid-lowering indicator as the protocol's statin proxy."""
    return (
        binary(row.get("lipid_lowering_medication_use")),
        "lipid_lowering_proxy",
    )


def prevent_risk(
    sex,
    age,
    tc,
    hdl,
    sbp,
    egfr,
    diabetes,
    smoking,
    bp_med,
    statin,
):
    """Calculate risk using the study-specified inputs and ranges."""
    # No exclusion based on baseline clinical CVD.
    predictors = {
        "sex": sex,
        "age": age,
        "tc": tc,
        "hdl": hdl,
        "sbp": sbp,
        "egfr": egfr,
        "diabetes": diabetes,
        "smoking": smoking,
        "bp_med": bp_med,
        "statin": statin,
    }

    for name, value in predictors.items():
        if value is None:
            return None, f"missing_{name}"

    if not 30 <= age <= 79:
        return None, "age_out_of_range"

    if not 15 <= egfr <= 150:
        return None, "egfr_out_of_range"

    limits = (
        (tc, 130, 320),
        (hdl, 20, 100),
        (sbp, 90, 200),
    )

    outside = any(
        not low <= value <= high
        for value, low, high in limits
    )

    if PREVENT_RANGE_POLICY == "strict":
        if outside:
            return None, "lipids_or_sbp_out_of_range"

    elif PREVENT_RANGE_POLICY == "winsorize":
        tc, hdl, sbp = (
            min(max(value, low), high)
            for value, low, high in limits
        )

    else:
        raise ValueError("Unsupported PREVENT_RANGE_POLICY")

    a = (age - 55) / 10

    # Convert cholesterol from mg/dL to mmol/L.
    tc = tc * 0.02586
    hdl = hdl * 0.02586

    n = tc - hdl - 3.5
    h = (hdl - 1.3) / 0.3

    sl = (min(sbp, 110) - 110) / 20
    sh = (max(sbp, 110) - 130) / 20

    gl = (min(egfr, 60) - 60) / -15
    gh = (max(egfr, 60) - 90) / -15

    b, t, d, c = bp_med, statin, diabetes, smoking

    terms = (
        1,
        a,
        n,
        h,
        sl,
        sh,
        d,
        c,
        gl,
        gh,
        b,
        t,
        b * sh,
        t * n,
        a * n,
        a * h,
        a * sh,
        a * d,
        a * c,
        a * gl,
    )

    log_odds = math.fsum(
        beta * term
        for beta, term in zip(
            PREVENT_COEFFICIENTS[sex],
            terms,
        )
    )

    # Numerically stable logistic transformation.
    if log_odds >= 0:
        probability = 1 / (1 + math.exp(-log_odds))
    else:
        exp_value = math.exp(log_odds)
        probability = exp_value / (1 + exp_value)

    status = "ok_winsorized" if outside else "ok"

    return 100 * probability, status


def derive(row):
    sex = binary(row.get("sex"))
    age = positive(row.get("age_years"))
    bmi = positive(row.get("bmi_kg_m2"))
    waist = positive(row.get("waist_cm"))

    hdl = positive(row.get("hdl_mg_dl"))
    tc = positive(row.get("tc_mg_dl"))
    tg = positive(row.get("tg_mg_dl"))

    # Use the existing first baseline BP reading.
    sbp = positive(row.get("sbp_mmhg"))
    dbp = positive(row.get("dbp_mmhg"))

    hba1c = positive(row.get("hba1c_pct"))

    bp_med = binary(
        row.get("antihypertensive_medication_use")
    )

    waist_flag = (
        None
        if sex is None
        else threshold(
            waist,
            WAIST_FEMALE_CM if sex == 0 else WAIST_MALE_CM,
        )
    )

    # Serum creatinine has already been converted to mg/dL upstream.
    egfr, equation = calculate_egfr(
        age,
        sex,
        positive(row.get("serum_creatinine_mg_dl")),
        positive(row.get("cystatin_c_mg_l")),
    )

    urine = calculate_uacr(row)
    g, risk = kidney_risk(egfr, urine["uacr_category"])

    statin, statin_source = statin_input(row)

    # Pass the selected eGFR directly into PREVENT.
    cvd_risk, prevent_status = prevent_risk(
        sex,
        age,
        tc,
        hdl,
        sbp,
        egfr,
        binary(row.get("diabetes_status")),
        binary(row.get("current_smoking")),
        bp_med,
        statin,
    )

    return {
        "mets": classify_mets(
            waist_flag,
            sex,
            hdl,
            tg,
            sbp,
            dbp,
            bp_med,
            hba1c,
        ),
        "metabolic_definition": (
            "CKM_2023_HbA1c_adapted_uniform_waist"
        ),
        "egfr_ml_min_1_73m2": egfr,
        "egfr_equation": equation,
        "egfr_primary_policy": "legacy_best_available",
        **urine,
        "egfr_category": g,
        "ckd_risk": risk,
        "albuminuria_ge30": (
            None
            if urine["uacr_category"] is None
            else int(urine["uacr_category"] != "A1")
        ),
        "hypertriglyceridemia": threshold(
            tg, TG_CKM_CUTOFF
        ),
        "obesity": threshold(bmi, BMI_CKM_CUTOFF),
        "abdominal_obesity": waist_flag,
        "prevent_statin_use": statin,
        "prevent_statin_source": statin_source,
        "prevent_cvd_10y_pct": cvd_risk,
        "prevent_cvd_10y_ge20": threshold(cvd_risk, 20),
        "prevent_status": prevent_status,
        "prevent_range_policy": PREVENT_RANGE_POLICY,
    }


def main():
    temporary = None
    seen = set()
    counts = Counter()

    try:
        with BASELINE.open(
            encoding="utf-8-sig", newline=""
        ) as source:
            reader = csv.DictReader(source, strict=True)

            if not reader.fieldnames:
                raise ValueError("Empty input or missing header")

            header = [
                name.strip()
                for name in reader.fieldnames
            ]

            if "" in header or len(set(header)) != len(header):
                raise ValueError("Empty or duplicate headers")

            reader.fieldnames = header

            absent = sorted(REQUIRED - set(header))

            if absent:
                raise ValueError(
                    f"Missing baseline columns: {absent}"
                )

            derived_columns = list(derive({}))

            output_header = [
                name
                for name in header
                if name not in DROP_COLUMNS
            ]

            columns = output_header + [
                name
                for name in derived_columns
                if name not in output_header
            ]

            with NamedTemporaryFile(
                mode="w",
                encoding="utf-8",
                newline="",
                dir=BASELINE.parent,
                delete=False,
            ) as target:
                temporary = Path(target.name)

                writer = csv.DictWriter(
                    target,
                    fieldnames=columns,
                )
                writer.writeheader()

                for row in reader:
                    if None in row or any(
                        value is None
                        for value in row.values()
                    ):
                        raise ValueError(
                            "Incorrect cell count on CSV line "
                            f"{reader.line_num}"
                        )

                    eid = row["eid"].strip()

                    if (
                        not re.fullmatch(r"[0-9]+", eid)
                        or eid in seen
                    ):
                        raise ValueError(
                            "Invalid or duplicate eid on CSV line "
                            f"{reader.line_num}"
                        )

                    seen.add(eid)

                    result = derive(row)

                    counts["rows"] += 1
                    counts[
                        f"prevent:{result['prevent_status']}"
                    ] += 1

                    row.update(result)

                    for name in DROP_COLUMNS:
                        row.pop(name, None)

                    writer.writerow({
                        name: "NA" if value is None else value
                        for name, value in row.items()
                    })

        temporary.replace(BASELINE)
        print(dict(counts))

    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


if __name__ == "__main__":
    main()