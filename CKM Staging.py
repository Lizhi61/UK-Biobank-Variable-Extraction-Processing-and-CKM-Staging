#!/usr/bin/env python3
"""Assign CKM stages after the four baseline-processing scripts."""

import csv
import math
import re
from collections import Counter
from pathlib import Path
from tempfile import NamedTemporaryFile

BASELINE = Path("UKB_baseline.csv")
MISSING = {"", "na", "n/a", "nan", "null", "none"}
FLAGS = (
    "obesity", "abdominal_obesity", "prediabetes",
    "hypertriglyceridemia", "hypertension_status", "diabetes_status",
    "mets", "clinical_cvd", "prevent_cvd_10y_ge20",
)
REQUIRED = {"eid", "ckd_risk", "hba1c_pct", *FLAGS}
OUTPUTS = ("ckm_stage", "ckm_stage_status")


def read_value(row, name, allowed=None):
    raw = row[name].strip()
    if raw.lower() in MISSING:
        return None
    try:
        value = float(raw)
    except ValueError as exc:
        raise ValueError(f"Non-numeric value in {name}: {raw!r}") from exc
    if not math.isfinite(value) or (allowed is not None and value not in allowed):
        raise ValueError(f"Invalid value in {name}: {raw!r}")
    return value


def tri_all(*values):
    return 0 if 0 in values else (None if None in values else 1)


def tri_any(*values):
    return 1 if 1 in values else (None if None in values else 0)


def tri_not(value):
    return None if value is None else 1 - value


def stage_rules(row):
    b, w, p, t, h, d, m, c, r = (
        read_value(row, name, (0, 1)) for name in FLAGS
    )
    k = read_value(row, "ckd_risk", (1, 2, 3, 4))
    a1c = read_value(row, "hba1c_pct")
    if a1c is not None and a1c <= 0:
        raise ValueError("HbA1c must be positive or missing")
    low = None if k is None else int(k == 1)
    moderate_high = None if k is None else int(k in (2, 3))
    very_high = None if k is None else int(k == 4)
    normoglycemia = None if a1c is None else int(a1c < 5.7)
    background = tri_any(b, w, p, t, h, d, m, moderate_high, very_high)
    rules = {
        # Stage 4 requires clinical CVD and at least one background risk.
        4: tri_all(c, background),
        # Stage 3 does not require the optional background-risk criteria.
        3: tri_all(tri_not(c), tri_any(very_high, r)),
        2: tri_all(
            tri_not(c), tri_not(very_high), tri_not(r),
            tri_any(t, h, d, m, moderate_high),
        ),
        1: tri_all(
            tri_not(c), tri_not(r), low, tri_not(h), tri_not(t),
            tri_not(d), tri_not(m), tri_any(b, w, p),
        ),
        0: tri_all(
            tri_not(c), tri_not(r), low, tri_not(b), tri_not(w),
            normoglycemia, tri_not(d), tri_not(h), tri_not(t),
        ),
    }
    return rules, c, background


def classify(row):
    rules, cvd, background = stage_rules(row)
    if cvd == 1 and background == 0:
        return None, "clinical_cvd_without_background_risk"
    higher_unknown = False
    # Assign only when every higher stage has been ruled out.
    for stage in (4, 3, 2, 1, 0):
        if rules[stage] == 1:
            if higher_unknown:
                return None, "higher_stage_uncertain"
            return stage, "classified"
        if rules[stage] is None:
            higher_unknown = True
    return None, "insufficient_data" if higher_unknown else "criteria_not_met"


def main():
    temporary = None
    seen = set()
    counts = Counter({f"stage_{stage}": 0 for stage in range(5)})
    try:
        with BASELINE.open(encoding="utf-8-sig", newline="") as source:
            reader = csv.DictReader(source, strict=True)
            if not reader.fieldnames:
                raise ValueError("Empty input or missing header")
            header = [name.strip() for name in reader.fieldnames]
            if "" in header or len(header) != len(set(header)):
                raise ValueError("Empty or duplicate headers")
            reader.fieldnames = header
            absent = sorted(REQUIRED - set(header))
            if absent:
                raise ValueError(f"Missing CKM input columns: {absent}")
            columns = header + [name for name in OUTPUTS if name not in header]
            with NamedTemporaryFile(
                mode="w", encoding="utf-8", newline="",
                dir=BASELINE.parent, delete=False,
            ) as target:
                temporary = Path(target.name)
                writer = csv.DictWriter(target, fieldnames=columns)
                writer.writeheader()
                for row in reader:
                    if None in row or any(value is None for value in row.values()):
                        raise ValueError(f"Incorrect cell count on CSV line {reader.line_num}")
                    eid = row["eid"].strip()
                    if not re.fullmatch(r"[0-9]+", eid) or eid in seen:
                        raise ValueError(f"Invalid or duplicate eid on CSV line {reader.line_num}")
                    seen.add(eid)
                    try:
                        stage, status = classify(row)
                    except ValueError as exc:
                        raise ValueError(f"CSV line {reader.line_num}: {exc}") from exc
                    row["ckm_stage"] = "NA" if stage is None else stage
                    row["ckm_stage_status"] = status
                    writer.writerow(row)
                    counts["total_rows"] += 1
                    counts["stage_NA" if stage is None else f"stage_{stage}"] += 1
                    counts[f"status:{status}"] += 1
        temporary.replace(BASELINE)
        print(dict(counts))
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


if __name__ == "__main__":
    main()