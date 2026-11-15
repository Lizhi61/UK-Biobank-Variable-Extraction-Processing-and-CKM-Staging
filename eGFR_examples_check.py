#!/usr/bin/env python3
"""Verify the three supplied eGFR examples against independent reference values."""

import csv
import math
from pathlib import Path

OUTPUT_PATH = Path("egfr_examples_check.csv")


def calculate_egfr(age, sex, scr, scys):
    """Use sex 0=female, 1=male; Scr in mg/dL; Scys in mg/L."""
    if age is None or age < 18 or sex not in (0, 1):
        return None, None

    scr = scr if scr is not None and scr > 0 else None
    scys = scys if scys is not None and scys > 0 else None

    female = sex == 0
    kappa = 0.7 if female else 0.9

    if scr is not None and scys is not None:
        x = scr / kappa
        y = scys / 0.8
        alpha = -0.248 if female else -0.207

        # CKD-EPI 2012 combined equation:
        # eGFR = 135
        #   * min(Scr/kappa, 1)**alpha
        #   * max(Scr/kappa, 1)**(-0.601)
        #   * min(Scys/0.8, 1)**(-0.375)
        #   * max(Scys/0.8, 1)**(-0.711)
        #   * 0.995**age
        #   * (0.969 if female else 1)
        #   * race_multiplier.
        #
        # kappa = 0.7 for females, 0.9 for males.
        # alpha = -0.248 for females, -0.207 for males.
        # race_multiplier is fixed at 1 under this study protocol.

        egfr = (
            135
            * min(x, 1) ** alpha
            * max(x, 1) ** -0.601
            * min(y, 1) ** -0.375
            * max(y, 1) ** -0.711
            * 0.995 ** age
            * (0.969 if female else 1)
        )

        return egfr, "CKD-EPI_2012_cr_cys_race_multiplier_1"

    if scr is not None:
        x = scr / kappa
        alpha = -0.329 if female else -0.411

        # CKD-EPI 2009 creatinine equation:
        # eGFR = 141
        #   * min(Scr/kappa, 1)**alpha
        #   * max(Scr/kappa, 1)**(-1.209)
        #   * 0.993**age
        #   * (1.018 if female else 1)
        #   * race_multiplier.
        #
        # kappa = 0.7 for females, 0.9 for males.
        # alpha = -0.329 for females, -0.411 for males.
        # race_multiplier is fixed at 1 under this study protocol.

        egfr = (
            141
            * min(x, 1) ** alpha
            * max(x, 1) ** -1.209
            * 0.993 ** age
            * (1.018 if female else 1)
        )

        return egfr, "CKD-EPI_2009_cr_race_multiplier_1"

    if scys is not None:
        y = scys / 0.8

        # CKD-EPI 2012 cystatin C equation:
        # eGFR = 133
        #   * min(Scys/0.8, 1)**(-0.499)
        #   * max(Scys/0.8, 1)**(-1.328)
        #   * 0.996**age
        #   * (0.932 if female else 1).
        #
        # This equation has no race coefficient.

        egfr = (
            133
            * min(y, 1) ** -0.499
            * max(y, 1) ** -1.328
            * 0.996 ** age
            * (0.932 if female else 1)
        )

        return egfr, "CKD-EPI_2012_cystatin_C"

    return None, None


def main():
    # Independent reference values were evaluated at higher precision.
    # Each tuple contains:
    # name, age, sex, Scr (umol/L), Scys (mg/L),
    # expected equation, expected eGFR.
    examples = (
        (
            "Example 1",
            60,
            1,
            88.40,
            None,
            "CKD-EPI_2009_cr_race_multiplier_1",
            81.44289751798162,
        ),
        (
            "Example 2",
            60,
            0,
            None,
            1.00,
            "CKD-EPI_2012_cystatin_C",
            72.46551140565252,
        ),
        (
            "Example 3",
            60,
            1,
            88.40,
            1.00,
            "CKD-EPI_2012_cr_cys_race_multiplier_1",
            80.04152272869973,
        ),
    )

    results = []

    for (
        name,
        age,
        sex,
        scr_umol,
        scys,
        expected_equation,
        expected,
    ) in examples:
        # Convert serum creatinine exactly once before calculation.
        scr_mg = None if scr_umol is None else scr_umol / 88.4

        calculated, equation = calculate_egfr(
            age, sex, scr_mg, scys
        )

        passed = (
            equation == expected_equation
            and calculated is not None
            and math.isclose(
                calculated,
                expected,
                rel_tol=1e-12,
                abs_tol=1e-10,
            )
        )

        if not passed:
            raise AssertionError(
                f"{name} failed: equation={equation}, "
                f"calculated={calculated}, expected={expected}"
            )

        results.append({
            "example": name,
            "age_years": age,
            "sex": "Female" if sex == 0 else "Male",
            "serum_creatinine_umol_l": scr_umol,
            "serum_creatinine_mg_dl": scr_mg,
            "cystatin_c_mg_l": scys,
            "equation": equation,
            "egfr_calculated": calculated,
            "egfr_expected": expected,
            "absolute_difference": abs(calculated - expected),
            "check": "PASS",
        })

        print(
            f"{name}: eGFR={calculated:.8f} "
            "mL/min/1.73 m2 | PASS"
        )
        print(f"  Equation: {equation}")

    with OUTPUT_PATH.open(
        "w", encoding="utf-8-sig", newline=""
    ) as target:
        writer = csv.DictWriter(
            target,
            fieldnames=list(results[0]),
        )
        writer.writeheader()

        for row in results:
            writer.writerow({
                key: "NA" if value is None else value
                for key, value in row.items()
            })

    print("All three examples passed.")
    print("Results:", OUTPUT_PATH.resolve())


if __name__ == "__main__":
    main()