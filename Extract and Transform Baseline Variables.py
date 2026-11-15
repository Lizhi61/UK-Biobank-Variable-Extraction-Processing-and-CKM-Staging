#!/usr/bin/env python3
"""Extract baseline CKM input variables from a wide UKB CSV."""

import csv
from datetime import date
import math
import re
from collections import Counter
from pathlib import Path
from tempfile import NamedTemporaryFile

INPUT_PATH = Path("visit_data.csv")
OUTPUT_PATH = Path("UKB_baseline.csv")

MISSING = {"", "na", "n/a", "nan", "null", "none"}

# Field ID, output name, optional range, optional reportability/urine assay flag field.
MEASURES = (
    (21003, "age_years", None, None),
    (21001, "bmi_kg_m2", None, None),
    (48, "waist_cm", None, None),
    (30750, "hba1c_mmol_mol", (15, 184), 30756),
    (30870, "tg_mmol_l", (0.10, 11.3), 30876),
    (30760, "hdl_mmol_l", (0.05, 4.65), 30766),
    (30690, "tc_mmol_l", (0.50, 18), 30696),
    (30700, "serum_creatinine_umol_l", (4, 4420), 30706),
    (30720, "cystatin_c_mg_l", (0.10, 8.99), 30726),
    (30500, "urine_albumin_mg_l", None, 30505),
    (30510, "urine_creatinine_umol_l", None, 30515),
)

URINE_DETECTION_LIMITS = {
    30500: ("<6.7", 6.7),
    30510: ("<88", 88.0),
}

ETHNICITY = dict.fromkeys((1, 1001, 1002, 1003), "White")
ETHNICITY.update(dict.fromkeys(
    (
        2, 3, 4, 5, 6,
        2001, 2002, 2003, 2004,
        3001, 3002, 3003, 3004,
        4001, 4002, 4003,
    ),
    "Non-White",
))

BP_FIELDS = {"sbp": (4080, 93), "dbp": (4079, 94)}

CONVERSIONS = (
    ("hba1c_pct", "hba1c_mmol_mol", 0.09148, 2.15),
    ("tg_mg_dl", "tg_mmol_l", 88.57, 0),
    ("hdl_mg_dl", "hdl_mmol_l", 38.67, 0),
    ("tc_mg_dl", "tc_mmol_l", 38.67, 0),
    ("serum_creatinine_mg_dl", "serum_creatinine_umol_l", 1 / 88.4, 0),
    ("urine_creatinine_g_l", "urine_creatinine_umol_l", 1 / 8840, 0),
)


def text(value):
    """Safely convert a cell value to stripped text."""
    return "" if value is None else str(value).strip()


def normalize_column(column):
    """Normalize supported column names to fieldID-instance.array."""
    column = text(column)

    if column in {"eid", "f.eid"}:
        return "eid"

    # The RAP sex field is not repeated across visits.
    if column == "p31":
        return "31-0.0"

    patterns = (
        r"(\d+)-(\d+)\.(\d+)",
        r"f\.(\d+)\.(\d+)\.(\d+)",
        r"p(\d+)_i(\d+)(?:_a(\d+))?",
    )

    for pattern in patterns:
        match = re.fullmatch(pattern, column)
        if match:
            field, instance, array = match.groups()
            return f"{int(field)}-{int(instance)}.{int(array or 0)}"

    return column


def measurement(raw, limits=None, reportability=""):
    raw = text(raw)
    reportability = text(reportability)

    if raw.lower() in MISSING:
        return None

    if raw.startswith(("<", ">", "≤", "≥")):
        return None

    # Specific urine flags are handled earlier in extract().
    # Accept only code 1 when blood reportability is provided.
    # Missing flags permit numeric extraction but do not establish assay quality.
    if reportability.lower() not in MISSING:
        try:
            if float(reportability) != 1:
                return None
        except ValueError:
            return None

    try:
        value = float(raw)
    except ValueError:
        return None

    if not math.isfinite(value) or value <= 0:
        return None

    if limits is not None and not limits[0] <= value <= limits[1]:
        return None

    return value


def category(raw, mapping):
    try:
        return mapping.get(float(text(raw)))
    except ValueError:
        return None


def extract(row, array_indices, eid, baseline_date):
    def raw(field, index=0):
        return text(row.get(f"{field}-0.{index}"))

    out = {"eid": eid, "baseline_date": baseline_date}

    for field, name, limits, flag in MEASURES:
        reportability = raw(flag) if flag is not None else ""

        if flag is not None:
            out[f"{name}_flag_raw"] = reportability or None

        detection_limit = URINE_DETECTION_LIMITS.get(field)

        if detection_limit:
            below_lod = (
                "".join(reportability.split()) == detection_limit[0]
            )
            out[f"{name}_raw"] = raw(field) or None
            out[f"{name}_lod_substituted"] = int(below_lod)

            # Preserve the original policy: substitute the detection limit and record the substitution.
            if below_lod:
                out[name] = detection_limit[1]
                continue

        out[name] = measurement(raw(field), limits, reportability)

    # For SBP and DBP separately, prefer automated readings; use manual readings only if no valid automated reading exists.
    for name, (auto_field, manual_field) in BP_FIELDS.items():
        candidates = {}

        for source, field in (
            ("automated", auto_field),
            ("manual", manual_field),
        ):
            candidates[source] = []

            for index in array_indices[field]:
                value = measurement(raw(field, index))
                out[f"{name}_{source}_a{index}_mmhg"] = value

                if value is not None:
                    candidates[source].append(value)

        out[f"{name}_mmhg"] = None
        out[f"{name}_source"] = None

        for source, values in candidates.items():
            if values:
                out[f"{name}_mmhg"] = sum(values) / len(values)
                out[f"{name}_source"] = source
                break

    for target, source, factor, offset in CONVERSIONS:
        value = out[source]
        out[target] = None if value is None else value * factor + offset

    out["sex"] = category(raw(31), {0: 0, 1: 1})
    out["ethnicity"] = category(raw(21000), ETHNICITY)

    smoking = category(
        raw(20116),
        {0: "never", 1: "former", 2: "current"},
    )
    out["smoking"] = smoking
    out["current_smoking"] = (
        None if smoking is None else int(smoking == "current")
    )

    return out


def write_baseline():
    temporary = None
    counts = Counter()
    seen_ids = set()

    try:
        with INPUT_PATH.open(
            encoding="utf-8-sig", newline=""
        ) as source:
            reader = csv.DictReader(source, strict=True)

            if not reader.fieldnames:
                raise ValueError("Input file is empty or has no header")

            header = [normalize_column(c) for c in reader.fieldnames]
            duplicates = [
                c for c, n in Counter(header).items() if n > 1
            ]

            if "" in header or duplicates:
                raise ValueError(
                    f"Empty column names or duplicate columns after normalization: {duplicates}"
                )

            reader.fieldnames = header

            required = {
                "eid", "53-0.0", "31-0.0",
                "21000-0.0", "20116-0.0",
            }
            required.update(
                f"{field}-0.0" for field, _, _, _ in MEASURES
            )

            absent = sorted(required - set(header))
            if absent:
                raise ValueError(
                    f"Missing required columns or unsupported column name format: {absent}"
                )

            absent_flags = [
                f"{flag}-0.0"
                for _, _, _, flag in MEASURES
                if flag is not None and f"{flag}-0.0" not in header
            ]
            if absent_flags:
                print(
                    "Note: the following flag columns are unavailable; full "
                    f"reportability/detection-limit verification is not possible: {absent_flags}"
                )

            array_indices = {}
            columns = ["eid", "baseline_date"]
            columns.extend(name for _, name, _, _ in MEASURES)

            for name, (auto_field, manual_field) in BP_FIELDS.items():
                for method, field in (
                    ("automated", auto_field),
                    ("manual", manual_field),
                ):
                    prefix = f"{field}-0."
                    indices = sorted(
                        int(column[len(prefix):])
                        for column in header
                        if column.startswith(prefix)
                        and column[len(prefix):].isdigit()
                    )
                    array_indices[field] = indices
                    columns.extend(
                        f"{name}_{method}_a{index}_mmhg"
                        for index in indices
                    )

                if not (
                    array_indices[auto_field]
                    or array_indices[manual_field]
                ):
                    raise ValueError(
                        f"All baseline {name} blood pressure columns are missing"
                    )

                columns.extend(
                    (f"{name}_mmhg", f"{name}_source")
                )

            columns.extend(
                target for target, _, _, _ in CONVERSIONS
            )
            columns.extend(
                ("sex", "ethnicity", "smoking", "current_smoking")
            )

            for field, name, _, flag in MEASURES:
                if flag is not None:
                    columns.append(f"{name}_flag_raw")
                if field in URINE_DETECTION_LIMITS:
                    columns.extend(
                        (f"{name}_raw", f"{name}_lod_substituted")
                    )

            with NamedTemporaryFile(
                mode="w",
                encoding="utf-8",
                newline="",
                dir=OUTPUT_PATH.parent,
                delete=False,
            ) as target:
                temporary = Path(target.name)
                writer = csv.DictWriter(target, fieldnames=columns)
                writer.writeheader()

                for row in reader:
                    counts["input_rows"] += 1

                    if None in row or any(
                        value is None for value in row.values()
                    ):
                        raise ValueError(
                            f"Column count on CSV line {reader.line_num} "
                            "does not match the header"
                        )

                    eid = text(row["eid"])

                    if not re.fullmatch(r"[0-9]+", eid):
                        raise ValueError(
                            f"CSV line {reader.line_num}: "
                            "eid is missing or invalid"
                        )

                    if eid in seen_ids:
                        raise ValueError(
                            f"CSV line {reader.line_num}: duplicate eid"
                        )
                    seen_ids.add(eid)

                    visit_date = text(row["53-0.0"])

                    if visit_date.lower() in MISSING:
                        counts["missing_baseline_date"] += 1
                        continue

                    try:
                        if not re.fullmatch(
                            r"[0-9]{4}-[0-9]{2}-[0-9]{2}",
                            visit_date,
                        ):
                            raise ValueError("Date must use YYYY-MM-DD format")

                        baseline_date = date.fromisoformat(
                            visit_date
                        ).isoformat()

                    except ValueError:
                        counts["invalid_baseline_date"] += 1
                        continue

                    result = extract(
                        row, array_indices, eid, baseline_date
                    )

                    counts["written_rows"] += 1
                    counts["urine_albumin_lod_substitutions"] += (
                        result["urine_albumin_mg_l_lod_substituted"]
                    )
                    counts["urine_creatinine_lod_substitutions"] += (
                        result[
                            "urine_creatinine_umol_l_lod_substituted"
                        ]
                    )

                    writer.writerow({
                        name: "NA" if value is None else value
                        for name, value in result.items()
                    })

        temporary.replace(OUTPUT_PATH)

        print({
            key: counts[key]
            for key in (
                "input_rows",
                "written_rows",
                "missing_baseline_date",
                "invalid_baseline_date",
                "urine_albumin_lod_substitutions",
                "urine_creatinine_lod_substitutions",
            )
        })

    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


if __name__ == "__main__":
    write_baseline()