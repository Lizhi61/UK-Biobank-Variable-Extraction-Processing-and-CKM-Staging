#!/usr/bin/env python3
"""Verify HbA1c conversion using the study-specified coefficients."""

import csv
from decimal import Decimal, InvalidOperation, ROUND_HALF_UP
from pathlib import Path

OUTPUT_PATH = Path("hba1c_examples_check.csv")

SLOPE = Decimal("0.0915")
INTERCEPT = Decimal("2.15")


def convert_hba1c(value):
    """Convert HbA1c from IFCC mmol/mol to NGSP percent."""
    if value is None:
        return None

    try:
        value = Decimal(str(value))
    except InvalidOperation as exc:
        raise ValueError("HbA1c must be numeric") from exc

    if not value.is_finite() or value <= 0:
        raise ValueError("HbA1c must be a finite positive number")

    # Complete conversion equation:
    # HbA1c (NGSP %) = 0.0915 * HbA1c (IFCC mmol/mol) + 2.15.
    # Preserve the unrounded value for subsequent threshold comparisons.
    return SLOPE * value + INTERCEPT


def main():
    # Example label, HbA1c in mmol/mol, expected percentage to two decimals.
    examples = (
        (
            "Normal range",
            30,
            "4.90",
        ),
        (
            "Near the prediabetes threshold, below 5.7%",
            38,
            "5.63",
        ),
        (
            "Near the prediabetes threshold, above 5.7%",
            39,
            "5.72",
        ),
        (
            "Prediabetes range",
            42,
            "5.99",
        ),
        (
            "Near the diabetes threshold, below 6.5%",
            47,
            "6.45",
        ),
        (
            "Near the diabetes threshold, above 6.5%",
            48,
            "6.54",
        ),
        (
            "Diabetes range",
            53,
            "7.00",
        ),
        (
            "Markedly elevated",
            64,
            "8.01",
        ),
    )

    results = []

    for index, (label, mmol_mol, expected) in enumerate(
        examples, start=1
    ):
        calculated = convert_hba1c(mmol_mol)

        # Use exact decimal arithmetic and conventional half-up rounding.
        displayed = str(
            calculated.quantize(
                Decimal("0.01"),
                rounding=ROUND_HALF_UP,
            )
        )

        # The expected values are rounded to two decimals.
        # Compare displayed values with the supplied table.
        if displayed != expected:
            raise AssertionError(
                f"Example {index} failed: raw={calculated}, "
                f"displayed={displayed}, expected={expected}"
            )

        results.append({
            "example": f"Example {index}",
            "example_range": label,
            "hba1c_mmol_mol": mmol_mol,
            "hba1c_pct": calculated,
            "hba1c_pct_display": displayed,
            "expected_pct_display": expected,
            "check": "PASS",
        })

        print(
            f"Example {index}: {mmol_mol} mmol/mol -> "
            f"{displayed}% | PASS"
        )

    with OUTPUT_PATH.open(
        "w",
        encoding="utf-8-sig",
        newline="",
    ) as target:
        writer = csv.DictWriter(
            target,
            fieldnames=list(results[0]),
        )
        writer.writeheader()
        writer.writerows(results)

    print("All eight examples passed.")
    print("Results:", OUTPUT_PATH.resolve())


if __name__ == "__main__":
    main()