#!/usr/bin/env python3
"""Audit assigned CKM stages and retain verified records with all covariates."""

import csv
from collections import Counter
from contextlib import ExitStack
import json
import math
from pathlib import Path
import re
from tempfile import TemporaryDirectory


INPUT = Path("UKB_baseline.csv")

OUTPUT_NAMES = {
    "defined": "UKB_baseline defined CKM.csv",
    "passed": "UKB_baseline for CKM.csv",
    "audit": "UKB_CKM_audit.csv",
    "failed": "UKB_CKM_failed.csv",
    "excluded": "UKB_CKM_excluded.csv",
    "summary": "UKB_CKM_audit_summary.json",
}

MISSING = {"", "na", "n/a", "nan", "null", "none"}

CVD = ("chd", "stroke", "af", "pad", "hf")

SAVED_FLAGS = (
    "obesity",
    "abdominal_obesity",
    "prediabetes",
    "hypertriglyceridemia",
    "hypertension_status",
    "diabetes_status",
    "mets",
    "clinical_cvd",
    "prevent_cvd_10y_ge20",
)

REQUIRED = {
    "eid",
    "sex",
    "bmi_kg_m2",
    "waist_cm",
    "hba1c_pct",
    "tg_mg_dl",
    "hdl_mg_dl",
    "sbp_mmhg",
    "dbp_mmhg",
    "egfr_ml_min_1_73m2",
    "urine_albumin_mg_l",
    "urine_creatinine_umol_l",
    "uacr_mg_g",
    "ckd_risk",
    "egfr_category",
    "uacr_category",
    "prevent_cvd_10y_pct",
    "antihypertensive_medication_use",
    "glucose_lowering_medication_use",
    "hypertension_self_report",
    "diabetes_self_report",
    "hypertension_prebaseline",
    "diabetes_prebaseline",
    "chd_self_report",
    "stroke_self_report",
    *SAVED_FLAGS,
    *(f"{name}_prebaseline" for name in CVD),
}

# Columns correspond to A1, A2 and A3.
G_A_RISK = {
    "G1": (1, 2, 3),
    "G2": (1, 2, 3),
    "G3a": (2, 3, 4),
    "G3b": (3, 4, 4),
    "G4": (4, 4, 4),
    "G5": (4, 4, 4),
}

AUDIT_FIELDS = (
    "audit_input_line",
    "audit_assigned_stage",
    "audit_result",
    "audit_reasons",
    "audit_current_rule",
    "audit_higher_stages_excluded",
    "audit_obvious_violation",
    "audit_current_false",
    "audit_current_unknown",
    "audit_higher_true",
    "audit_higher_unknown",
    "audit_covariate_conflicts",
    "audit_unverified_flags",
    "audit_stage4",
    "audit_stage3",
    "audit_stage2",
    "audit_stage1",
    "audit_stage0",
    "audit_criteria_json",
)


def missing(value):
    """Recognize missing values without treating numeric zero as missing."""
    return value is None or str(value).strip().lower() in MISSING


def read_number(row, name, allowed=None, minimum=None, positive=False):
    """Read a finite numeric value and validate its coding."""
    if missing(row[name]):
        return None

    raw = str(row[name]).strip()

    try:
        value = float(raw)
    except ValueError as exc:
        raise ValueError(f"Non-numeric {name}: {raw!r}") from exc

    if not math.isfinite(value):
        raise ValueError(f"Non-finite {name}: {raw!r}")

    if allowed is not None and value not in allowed:
        raise ValueError(f"Invalid code in {name}: {raw!r}")

    if minimum is not None and value < minimum:
        raise ValueError(f"Invalid negative value in {name}: {raw!r}")

    if positive and value <= 0:
        raise ValueError(f"Non-positive {name}: {raw!r}")

    return value


def tri_all(values):
    """Three-valued AND: a definite failure overrides missing information."""
    values = tuple(values)
    return 0 if 0 in values else (None if None in values else 1)


def tri_any(values):
    """Three-valued OR: a definite positive overrides missing information."""
    values = tuple(values)
    return 1 if 1 in values else (None if None in values else 0)


def negate(value):
    return None if value is None else 1 - value


def ge(value, cutoff):
    return None if value is None else int(value >= cutoff)


def stage_value(raw):
    """Accept only stages 0, 1, 2, 3 and 4."""
    if missing(raw):
        return None

    try:
        value = float(raw)
    except ValueError as exc:
        raise ValueError("invalid_stage") from exc

    if not math.isfinite(value) or value not in (0, 1, 2, 3, 4):
        raise ValueError("invalid_stage")

    return int(value)


def kidney_categories(egfr, uacr):
    """Rebuild KDIGO G category, A category and risk code."""
    g = None if egfr is None else next(
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
    )

    a = (
        None
        if uacr is None
        else "A1" if uacr < 30
        else "A2" if uacr <= 300
        else "A3"
    )

    # G4/G5 imply very high risk even when the A category is missing.
    if g in ("G4", "G5"):
        risk = 4
    elif g is not None and a is not None:
        risk = G_A_RISK[g][int(a[1]) - 1]
    else:
        risk = None

    return g, a, risk


def audited_features(row):
    """Rebuild stage inputs independently from saved measurements and evidence."""
    flag = lambda name: read_number(row, name, allowed=(0, 1))
    measure = lambda name: read_number(row, name, positive=True)

    # UKB coding: 0 = female, 1 = male.
    sex = flag("sex")

    bmi = measure("bmi_kg_m2")
    waist = measure("waist_cm")
    a1c = measure("hba1c_pct")
    tg = measure("tg_mg_dl")
    hdl = measure("hdl_mg_dl")

    # Use the baseline BP values already prepared upstream.
    sbp = measure("sbp_mmhg")
    dbp = measure("dbp_mmhg")

    bp_med = flag("antihypertensive_medication_use")
    glucose_med = flag("glucose_lowering_medication_use")

    obesity = ge(bmi, 25)

    # Apply the same waist thresholds to every ethnicity.
    abdominal = (
        None
        if sex is None
        else ge(waist, 88 if sex == 0 else 102)
    )

    # Include prebaseline hospital evidence, self-report, medication and BP.
    hypertension = tri_any((
        flag("hypertension_prebaseline"),
        flag("hypertension_self_report"),
        bp_med,
        ge(sbp, 130),
        ge(dbp, 80),
    ))

    # Include prebaseline hospital evidence, self-report, medication and HbA1c.
    diabetes = tri_any((
        flag("diabetes_prebaseline"),
        flag("diabetes_self_report"),
        glucose_med,
        ge(a1c, 6.5),
    ))

    if diabetes == 1 or (
        a1c is not None and not 5.7 <= a1c < 6.5
    ):
        prediabetes = 0
    else:
        prediabetes = (
            1 if diabetes == 0 and a1c is not None else None
        )

    low_hdl = (
        None
        if sex is None or hdl is None
        else int(hdl < (50 if sex == 0 else 40))
    )

    # MetS BP criteria use measured BP or BP medication.
    components = (
        abdominal,
        low_hdl,
        ge(tg, 150),
        tri_any((ge(sbp, 130), ge(dbp, 80), bp_med)),
        ge(a1c, 5.7),
    )

    positives = components.count(1)

    mets = (
        1 if positives >= 3
        else 0 if positives + components.count(None) < 3
        else None
    )

    clinical = tri_any(
        tuple(flag(f"{name}_prebaseline") for name in CVD)
        + (
            flag("chd_self_report"),
            flag("stroke_self_report"),
        )
    )

    albumin = read_number(row, "urine_albumin_mg_l", minimum=0)
    urine_cr = measure("urine_creatinine_umol_l")

    # Preserve the agreed upstream detection-limit substitution policy.
    # UACR (mg/g) = albumin (mg/L) / creatinine (umol/L) * 8840.
    uacr = (
        None
        if albumin is None or urine_cr is None
        else albumin / urine_cr * 8840.0
    )

    if uacr is not None and not math.isfinite(uacr):
        raise ValueError("Non-finite calculated UACR")

    # Use saved eGFR and PREVENT values without recalculating either model.
    egfr = measure("egfr_ml_min_1_73m2")
    g, a, risk = kidney_categories(egfr, uacr)

    probability = read_number(
        row, "prevent_cvd_10y_pct", minimum=0
    )

    if probability is not None and probability > 100:
        raise ValueError("PREVENT percentage exceeds 100")

    expected = {
        "obesity": obesity,
        "abdominal_obesity": abdominal,
        "prediabetes": prediabetes,
        "hypertriglyceridemia": ge(tg, 135),
        "hypertension_status": hypertension,
        "diabetes_status": diabetes,
        "mets": mets,
        "clinical_cvd": clinical,
        "prevent_cvd_10y_ge20": ge(probability, 20),
        "ckd_risk": risk,
        "egfr_category": g,
        "uacr_category": a,
        "uacr_mg_g": uacr,
    }

    conflicts = []
    unverified = []

    # Compare saved variables with independently rebuilt values.
    for name, calculated in expected.items():
        if name in ("egfr_category", "uacr_category"):
            saved = None if missing(row[name]) else row[name].strip()

            valid_categories = (
                set(G_A_RISK)
                if name == "egfr_category"
                else {"A1", "A2", "A3"}
            )

            if saved is not None and saved not in valid_categories:
                raise ValueError(f"Invalid category in {name}")

        else:
            allowed = (
                (1, 2, 3, 4)
                if name == "ckd_risk"
                else (0, 1) if name in SAVED_FLAGS
                else None
            )

            saved = read_number(
                row, name, allowed=allowed, minimum=0
            )

        if saved is not None and calculated is not None:
            same = (
                math.isclose(
                    saved,
                    calculated,
                    rel_tol=1e-9,
                    abs_tol=1e-9,
                )
                if name == "uacr_mg_g"
                else saved == calculated
            )

            if not same:
                conflicts.append(
                    f"{name}:saved={saved},calculated={calculated}"
                )

        elif saved is not None and calculated is None:
            unverified.append(name)

    expected["normoglycemia"] = (
        None if a1c is None else int(a1c < 5.7)
    )

    return expected, conflicts, unverified


def criterion_sets(f):
    """List every necessary condition separately for each stage."""
    b = f["obesity"]
    w = f["abdominal_obesity"]
    p = f["prediabetes"]
    t = f["hypertriglyceridemia"]
    h = f["hypertension_status"]
    d = f["diabetes_status"]
    m = f["mets"]
    c = f["clinical_cvd"]
    r = f["prevent_cvd_10y_ge20"]
    k = f["ckd_risk"]

    low = None if k is None else int(k == 1)
    moderate = None if k is None else int(k in (2, 3))
    very_high = None if k is None else int(k == 4)

    # Stage 4 requires one of these nine risks; Stage 3 does not.
    background = tri_any((
        b, w, p, t, h, d, m, moderate, very_high
    ))

    return {
        4: {
            "clinical_cvd_present": c,
            "background_risk_present": background,
        },
        3: {
            "no_clinical_cvd": negate(c),
            "very_high_ckd_or_prevent_ge20": tri_any((
                very_high, r
            )),
        },
        2: {
            "no_clinical_cvd": negate(c),
            "no_very_high_ckd": negate(very_high),
            "prevent_lt20": negate(r),
            "metabolic_risk_or_moderate_high_ckd": tri_any((
                t, h, d, m, moderate
            )),
        },
        1: {
            "no_clinical_cvd": negate(c),
            "prevent_lt20": negate(r),
            "low_ckd": low,
            "no_hypertension": negate(h),
            "tg_lt135": negate(t),
            "no_diabetes": negate(d),
            "no_mets": negate(m),
            "adiposity_or_prediabetes": tri_any((b, w, p)),
        },
        0: {
            "no_clinical_cvd": negate(c),
            "prevent_lt20": negate(r),
            "low_ckd": low,
            "bmi_lt25": negate(b),
            "normal_waist": negate(w),
            "hba1c_lt5_7": f["normoglycemia"],
            "no_diabetes": negate(d),
            "no_hypertension": negate(h),
            "tg_lt135": negate(t),
        },
    }


def audit_patient(row, assigned, line_number):
    """Check the current stage, higher stages and covariate consistency."""
    details = {name: "" for name in AUDIT_FIELDS}
    details.update(
        audit_input_line=line_number,
        audit_assigned_stage=assigned,
    )

    try:
        features, conflicts, unverified = audited_features(row)
        criteria = criterion_sets(features)
    except ValueError as exc:
        details.update(
            audit_result="invalid_input",
            audit_reasons=str(exc),
        )
        return False, details, ["invalid_input"]

    # Evaluate every stage in descending order.
    rules = {
        stage: tri_all(criteria[stage].values())
        for stage in (4, 3, 2, 1, 0)
    }

    current = rules[assigned]
    higher = list(range(4, assigned, -1))

    higher_true = [
        stage for stage in higher if rules[stage] == 1
    ]

    higher_unknown = [
        stage for stage in higher if rules[stage] is None
    ]

    # Stage 4 has no higher stages, so higher_clear is automatically 1.
    higher_clear = (
        0 if higher_true
        else None if higher_unknown
        else 1
    )

    obvious = (
        current == 0
        or bool(higher_true)
        or bool(conflicts)
    )

    # Missing information never proves that a higher stage is absent.
    passed = (
        current == 1
        and higher_clear == 1
        and not conflicts
    )

    reasons = []

    if current == 0:
        reasons.append("current_stage_failed")
    elif current is None:
        reasons.append("current_stage_uncertain")

    if higher_true:
        reasons.append("higher_stage_satisfied")

    if higher_unknown:
        reasons.append("higher_stage_uncertain")

    if conflicts:
        reasons.append("covariate_conflict")

    details.update(
        audit_result=(
            "passed" if passed
            else "failed" if obvious
            else "uncertain"
        ),
        audit_reasons=";".join(reasons),
        audit_current_rule=current,
        audit_higher_stages_excluded=higher_clear,
        audit_obvious_violation=int(obvious),
        audit_current_false=";".join(
            name
            for name, value in criteria[assigned].items()
            if value == 0
        ),
        audit_current_unknown=";".join(
            name
            for name, value in criteria[assigned].items()
            if value is None
        ),
        audit_higher_true=";".join(map(str, higher_true)),
        audit_higher_unknown=";".join(map(str, higher_unknown)),
        audit_covariate_conflicts=";".join(conflicts),
        audit_unverified_flags=";".join(unverified),
        audit_criteria_json=json.dumps(
            criteria,
            ensure_ascii=True,
            sort_keys=True,
        ),
    )

    details.update({
        f"audit_stage{stage}": rules[stage]
        for stage in (4, 3, 2, 1, 0)
    })

    return passed, details, reasons


def csv_values(row):
    """Write unknown audit results as NA."""
    return {
        name: "NA" if value is None else value
        for name, value in row.items()
    }


def main():
    counts = Counter()
    reasons = Counter()
    exclusions = Counter()

    by_stage = {
        str(stage): {
            "defined": 0,
            "passed": 0,
            "not_passed": 0,
        }
        for stage in range(5)
    }

    with INPUT.open(encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source, strict=True)
        header = reader.fieldnames

        if (
            not header
            or "" in header
            or len(header) != len(set(header))
        ):
            raise ValueError(
                "Missing, empty, or duplicate column names"
            )

        # Accept either CKM_stage or ckm_stage without renaming it.
        stage_columns = [
            name
            for name in header
            if name.lower() == "ckm_stage"
        ]

        if not stage_columns:
            raise ValueError("Missing ckm_stage or CKM_stage")

        absent = sorted(REQUIRED - set(header))

        if absent:
            raise ValueError(
                "Missing covariates needed for independent checking: "
                f"{absent}"
            )

        if set(header) & set(AUDIT_FIELDS):
            raise ValueError(
                "Input contains reserved audit column names"
            )

        # Stage outputs in temporary files until processing is complete.
        with TemporaryDirectory(
            prefix=".ckm_audit_",
            dir=INPUT.parent,
        ) as temporary:
            temporary = Path(temporary)

            with ExitStack() as stack:
                writers = {}

                for kind in (
                    "defined", "passed", "audit", "failed", "excluded"
                ):
                    path = temporary / OUTPUT_NAMES[kind]

                    stream = stack.enter_context(
                        path.open(
                            "w",
                            encoding="utf-8-sig",
                            newline="",
                        )
                    )

                    if kind in ("defined", "passed"):
                        columns = header
                    elif kind == "audit":
                        columns = ["eid", *AUDIT_FIELDS]
                    else:
                        columns = [*header, *AUDIT_FIELDS]

                    writers[kind] = csv.DictWriter(
                        stream, fieldnames=columns
                    )
                    writers[kind].writeheader()

                seen = set()

                for row in reader:
                    counts["input_rows"] += 1

                    if None in row or any(
                        value is None for value in row.values()
                    ):
                        raise ValueError(
                            "Incorrect cell count on CSV line "
                            f"{reader.line_num}"
                        )

                    eid = row["eid"].strip()

                    if not re.fullmatch(r"[0-9]+", eid) or eid in seen:
                        raise ValueError(
                            "Invalid or duplicate eid on CSV line "
                            f"{reader.line_num}"
                        )

                    seen.add(eid)

                    try:
                        stages = [
                            stage_value(row[name])
                            for name in stage_columns
                        ]

                        if len(set(stages)) > 1:
                            raise ValueError(
                                "conflicting_stage_columns"
                            )

                        assigned = stages[0]
                        exclusion = (
                            "missing_stage"
                            if assigned is None
                            else None
                        )

                    except ValueError as exc:
                        assigned = None
                        exclusion = str(exc)

                    if exclusion:
                        exclusions[exclusion] += 1

                        diagnostic = {
                            name: "" for name in AUDIT_FIELDS
                        }

                        diagnostic.update(
                            audit_input_line=reader.line_num,
                            audit_result="excluded",
                            audit_reasons=exclusion,
                        )

                        writers["excluded"].writerow({
                            **row, **diagnostic
                        })
                        continue

                    counts["defined_stage_rows"] += 1
                    by_stage[str(assigned)]["defined"] += 1

                    # Retain all clearly staged records before auditing.
                    writers["defined"].writerow(row)

                    passed, detail, patient_reasons = audit_patient(
                        row,
                        assigned,
                        reader.line_num,
                    )

                    writers["audit"].writerow(
                        csv_values({"eid": row["eid"], **detail})
                    )

                    reasons.update(patient_reasons)
                    counts[f"audit_{detail['audit_result']}"] += 1

                    if passed:
                        # Preserve all original columns and values.
                        writers["passed"].writerow(row)
                        by_stage[str(assigned)]["passed"] += 1
                    else:
                        writers["failed"].writerow(
                            csv_values({**row, **detail})
                        )
                        by_stage[str(assigned)]["not_passed"] += 1

            report = {
                "input_file": str(INPUT),
                "input_rows": counts["input_rows"],
                "defined_stage_rows": counts["defined_stage_rows"],
                "excluded_rows": sum(exclusions.values()),
                "passed_rows": counts["audit_passed"],
                "not_passed_rows": (
                    counts["defined_stage_rows"]
                    - counts["audit_passed"]
                ),
                "audit_outcomes": {
                    name: counts[f"audit_{name}"]
                    for name in (
                        "passed",
                        "failed",
                        "uncertain",
                        "invalid_input",
                    )
                },
                "excluded_reasons": dict(exclusions),
                "failure_reasons": dict(reasons),
                "by_assigned_stage": by_stage,
                "definitions": {
                    "stage3_background_required": False,
                    "stage4_background_required": True,
                    "hospital_evidence_in_hypertension_and_diabetes": True,
                },
                "validation_scope": (
                    "Stage rules and saved covariates; existing eGFR "
                    "and PREVENT values are not recalculated."
                ),
                "outputs": OUTPUT_NAMES,
            }

            summary_path = temporary / OUTPUT_NAMES["summary"]
            summary_path.write_text(
                json.dumps(report, indent=2, ensure_ascii=True),
                encoding="utf-8",
            )

            # Replace outputs only after the entire input is processed.
            for name in OUTPUT_NAMES.values():
                (temporary / name).replace(INPUT.with_name(name))

    print(json.dumps(report, indent=2, ensure_ascii=True))


if __name__ == "__main__":
    main()