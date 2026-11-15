#!/usr/bin/env python3
"""Add baseline disease phenotypes to an existing UKB baseline CSV."""

import csv
from collections import Counter
from datetime import date
import math
from pathlib import Path
import re
from tempfile import NamedTemporaryFile


VISITS = Path("visit_data.csv")
BASELINE = Path("UKB_baseline.csv")

# Explicit study definitions.
INCLUDE_BASELINE_DAY = True
INCLUDE_SELF_REPORTED_CVD = True
EXCLUDE_CHD_ANEURYSM = True

SBP_CUTOFF = 130.0
DBP_CUTOFF = 80.0
HBA1C_DIABETES_PCT = 6.5
HBA1C_PREDIABETES_PCT = 5.7

CVD = ("chd", "stroke", "af", "pad", "hf")
DISEASES = CVD + ("hypertension", "diabetes")

# Diagnosis field, corresponding date field, ICD version.
FIELD_PAIRS = (
    (41270, 41280, 10),
    (41271, 41281, 9),
)

MISSING = {"", "na", "n/a", "nan", "null", "none"}

# Always retain all PAD codes supplied in the original study definition.
PAD_ICD10_REQUIRED = {
    "I700", "I7000", "I7001",
    "I702", "I7020", "I7021",
    "I708", "I7080",
    "I709", "I7090",
    "I738", "I739",
}

PAD_ICD9_REQUIRED = {
    "4400", "4402", "4438", "4439",
}

# Always retain all additional HF codes supplied originally.
HF_ICD10_REQUIRED = {
    "I110", "I130", "I132",
    "I420", "I425", "I428", "I429",
}

REQUIRED_BASELINE = {
    "eid",
    "baseline_date",
    "sex",
    "sbp_mmhg",
    "dbp_mmhg",
    "hba1c_pct",
    "antihypertensive_medication_use",
    "glucose_lowering_medication_use",
}


def text(value):
    return "" if value is None else str(value).strip()


def missing(value):
    return text(value).lower() in MISSING


def normalize_column(value):
    value = text(value)

    if value in {"eid", "f.eid"}:
        return "eid"

    for pattern in (
        r"(\d+)-(\d+)\.(\d+)",
        r"f\.(\d+)\.(\d+)\.(\d+)",
        r"p(\d+)_i(\d+)(?:_a(\d+))?",
    ):
        match = re.fullmatch(pattern, value)

        if match:
            field, instance, array = match.groups()
            return f"{int(field)}-{int(instance)}.{int(array or 0)}"

    # Support flattened arrays from singular RAP fields.
    match = re.fullmatch(r"p(\d+)_a(\d+)", value)

    if match:
        field, array = match.groups()
        return f"{int(field)}-0.{int(array)}"

    return value


def reader_for(source):
    reader = csv.DictReader(source, strict=True)

    if not reader.fieldnames:
        raise ValueError("Empty input or missing header")

    header = [normalize_column(c) for c in reader.fieldnames]

    duplicates = [
        column
        for column, count in Counter(header).items()
        if count > 1
    ]

    if "" in header or duplicates:
        raise ValueError(
            "Empty or duplicate/alias-conflicting headers: "
            f"{duplicates}"
        )

    if "eid" not in header:
        raise ValueError("Missing eid or f.eid column")

    reader.fieldnames = header
    return reader, header


def checked_eid(row, line_number, seen):
    if None in row or any(value is None for value in row.values()):
        raise ValueError(
            f"Incorrect cell count on CSV line {line_number}"
        )

    eid = text(row["eid"])

    if not re.fullmatch(r"[0-9]+", eid) or eid in seen:
        raise ValueError(
            f"Invalid or duplicate eid on CSV line {line_number}"
        )

    seen.add(eid)
    return eid


def parse_date(value):
    value = text(value)

    if not re.fullmatch(r"[0-9]{4}-[0-9]{2}-[0-9]{2}", value):
        return None

    try:
        return date.fromisoformat(value)
    except ValueError:
        return None


def number(value):
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None

    return value if math.isfinite(value) else None


def binary(value):
    value = number(value)
    return int(value) if value in (0, 1) else None


def positive_number(value):
    value = number(value)
    return value if value is not None and value > 0 else None


def any_positive(*flags):
    # Missing evidence remains unknown unless another criterion is positive.
    return 1 if 1 in flags else (None if None in flags else 0)


def threshold(value, cutoff):
    return None if value is None else int(value >= cutoff)


def normalize_icd(value, version):
    value = text(value).upper()

    if version == 10:
        pattern = r"[A-Z][0-9]{2}(?:\.?[A-Z0-9]{1,4})?"
    else:
        pattern = (
            r"(?:[0-9]{3}(?:\.?[0-9]{1,2})?"
            r"|[EV][0-9]{2,4}(?:\.[0-9])?)"
        )

    if re.fullmatch(pattern, value):
        return value.replace(".", "")

    return None


def outcomes_for(code, version):
    group = code[:3]
    outcomes = set()

    if version == 10:
        if group in {"I20", "I21", "I22", "I23", "I24", "I25"}:
            if (
                not EXCLUDE_CHD_ANEURYSM
                or not code.startswith("I253")
            ):
                outcomes.add("chd")

        if group in {"I60", "I61", "I63", "I64"}:
            outcomes.add("stroke")

        if group == "I48":
            outcomes.add("af")

        # Retain the original exact codes and existing prefix matching.
        if (
            code in PAD_ICD10_REQUIRED
            or code.startswith(("I702", "I738", "I739"))
        ):
            outcomes.add("pad")

        if group == "I50" or code in HF_ICD10_REQUIRED:
            outcomes.add("hf")

        # Keep the revised hypertension definition.
        if group in {"I10", "I11", "I12", "I13", "I15"}:
            outcomes.add("hypertension")

        # Keep the revised non-gestational diabetes definition.
        if group in {"E10", "E11", "E12", "E13", "E14"}:
            outcomes.add("diabetes")

    else:
        if group in {"410", "411", "412", "413", "414"}:
            if (
                not EXCLUDE_CHD_ANEURYSM
                or not code.startswith("4141")
            ):
                outcomes.add("chd")

        # Preserve every stroke group supplied originally.
        if group in {
            "430", "431", "432", "433", "434", "435", "436",
        }:
            outcomes.add("stroke")

        if code.startswith("4273"):
            outcomes.add("af")

        # Retain the original exact codes and existing prefix matching.
        if (
            code in PAD_ICD9_REQUIRED
            or code.startswith(("4402", "4438", "4439"))
        ):
            outcomes.add("pad")

        if group == "428":
            outcomes.add("hf")

        # Keep the revised hypertension definition.
        if group in {"401", "402", "403", "404", "405"}:
            outcomes.add("hypertension")

        # Keep the revised non-gestational diabetes definition.
        if group == "250":
            outcomes.add("diabetes")

    return outcomes


def array_columns(header, field):
    pattern = re.compile(rf"{field}-0\.([0-9]+)")

    return {
        int(match.group(1)): column
        for column in header
        if (match := pattern.fullmatch(column))
    }


def paired_columns(header):
    pairs = []

    for diagnosis, dates, version in FIELD_PAIRS:
        codes = array_columns(header, diagnosis)
        times = array_columns(header, dates)

        if not codes and not times:
            print(
                f"Warning: ICD-{version} hospital fields are absent"
            )
            continue

        if not codes or codes.keys() != times.keys():
            raise ValueError(
                f"Unpaired hospital columns: {diagnosis}/{dates}"
            )

        if set(codes) != set(range(max(codes) + 1)):
            raise ValueError(
                f"Non-contiguous hospital array indices: {diagnosis}"
            )

        pairs.extend(
            (codes[index], times[index], version)
            for index in sorted(codes)
        )

    if not pairs:
        raise ValueError(
            "No paired hospital diagnosis/date columns were provided"
        )

    # Pairing cannot establish whether trailing slots were omitted.
    return pairs


def questionnaire_flag(row, columns, targets, complete):
    codes = {
        number(row[column])
        for column in columns
        if not missing(row[column])
    }

    if codes & targets:
        return 1

    valid = (
        codes == {-7}
        or (bool(codes) and codes <= {1, 2, 3, 4})
    )

    return 0 if complete and valid else None


def diabetes_self_report(doctor, pregnancy_only, sex):
    # Conflicting pregnancy information must not create a negative.
    if pregnancy_only == 1 and (sex == 1 or doctor == 0):
        return None

    if doctor != 1:
        return doctor

    # Missing pregnancy information is expected in men.
    if sex == 1 or pregnancy_only == 0:
        return 1

    if pregnancy_only == 1:
        return 0

    return None


def read_diagnoses():
    records = {}
    seen = set()

    with VISITS.open(encoding="utf-8-sig", newline="") as source:
        reader, header = reader_for(source)
        pairs = paired_columns(header)

        q_columns = array_columns(header, 6150)
        q_complete = set(range(4)) <= set(q_columns)

        if not q_complete:
            print(
                "Warning: incomplete field 6150; "
                "negative self-reports remain unknown"
            )

        for field in (2443, 4041):
            if f"{field}-0.0" not in header:
                print(f"Warning: baseline field {field} is absent")

        for row in reader:
            eid = checked_eid(row, reader.line_num, seen)

            earliest = {}
            undated = set()
            data_issue = False

            for code_col, date_col, version in pairs:
                if missing(row[code_col]):
                    data_issue |= not missing(row[date_col])
                    continue

                code = normalize_icd(row[code_col], version)

                if code is None:
                    data_issue = True
                    continue

                outcomes = outcomes_for(code, version)

                if not outcomes:
                    continue

                recorded = parse_date(row[date_col])

                for outcome in outcomes:
                    if recorded is None:
                        undated.add(outcome)

                    elif (
                        outcome not in earliest
                        or recorded < earliest[outcome]
                    ):
                        earliest[outcome] = recorded

            q = [
                q_columns[index]
                for index in sorted(q_columns)
            ]

            reports = tuple(
                questionnaire_flag(
                    row, q, targets, q_complete
                )
                for targets in ({1, 2}, {3}, {4})
            )

            records[eid] = (
                earliest,
                undated,
                data_issue,
                reports,
                binary(row.get("2443-0.0")),
                binary(row.get("4041-0.0")),
            )

    return records


def phenotype_variables(record, baseline_date, row):
    if record is None:
        earliest = {}
        undated = set()
        issue = False
        chd_report = None
        stroke_report = None
        hypertension_report = None
        doctor = None
        pregnancy_only = None

    else:
        (
            earliest,
            undated,
            issue,
            reports,
            doctor,
            pregnancy_only,
        ) = record

        chd_report, stroke_report, hypertension_report = reports

    result = {
        "hospital_source_matched": int(record is not None)
    }

    result["hospital_data_issue"] = (
        int(issue) if record is not None else None
    )

    for outcome in DISEASES:
        first = earliest.get(outcome)

        before = (
            first is not None
            and baseline_date is not None
            and (
                first <= baseline_date
                if INCLUDE_BASELINE_DAY
                else first < baseline_date
            )
        )

        if before:
            flag = 1

        elif (
            record is None
            or baseline_date is None
            or outcome in undated
            or issue
        ):
            flag = None

        else:
            flag = 0

        result[f"{outcome}_first_date"] = (
            first.isoformat() if first else None
        )

        result[f"{outcome}_prebaseline"] = flag

        result[f"{outcome}_hospital_undated"] = (
            int(outcome in undated)
            if record is not None
            else None
        )

    result["chd_self_report"] = chd_report
    result["stroke_self_report"] = stroke_report
    result["hypertension_self_report"] = hypertension_report

    result["clinical_cvd_hospital"] = any_positive(
        *(result[f"{d}_prebaseline"] for d in CVD)
    )

    result["clinical_cvd"] = (
        any_positive(
            result["clinical_cvd_hospital"],
            chd_report,
            stroke_report,
        )
        if INCLUDE_SELF_REPORTED_CVD
        else result["clinical_cvd_hospital"]
    )

    sex = binary(row.get("sex"))

    result["diabetes_self_report_conflict"] = (
        int(
            pregnancy_only == 1
            and (sex == 1 or doctor == 0)
        )
        if record is not None
        else None
    )

    diabetes_report = diabetes_self_report(
        doctor, pregnancy_only, sex
    )
    result["diabetes_self_report"] = diabetes_report

    # Use existing baseline measurements without recalculating BP.
    sbp = positive_number(row.get("sbp_mmhg"))
    dbp = positive_number(row.get("dbp_mmhg"))
    hba1c = positive_number(row.get("hba1c_pct"))

    result["bp_elevated"] = any_positive(
        threshold(sbp, SBP_CUTOFF),
        threshold(dbp, DBP_CUTOFF),
    )

    result["hba1c_diabetes_range"] = threshold(
        hba1c, HBA1C_DIABETES_PCT
    )

    result["hypertension_status"] = any_positive(
        result["hypertension_prebaseline"],
        hypertension_report,
        binary(row.get("antihypertensive_medication_use")),
        result["bp_elevated"],
    )

    diabetes = any_positive(
        result["diabetes_prebaseline"],
        diabetes_report,
        binary(row.get("glucose_lowering_medication_use")),
        result["hba1c_diabetes_range"],
    )
    result["diabetes_status"] = diabetes

    if diabetes == 1 or (
        hba1c is not None
        and not (
            HBA1C_PREDIABETES_PCT
            <= hba1c
            < HBA1C_DIABETES_PCT
        )
    ):
        result["prediabetes"] = 0

    elif diabetes == 0 and hba1c is not None:
        result["prediabetes"] = 1

    else:
        result["prediabetes"] = None

    return result


def main():
    records = read_diagnoses()
    new_columns = list(phenotype_variables(None, None, {}))

    counts = Counter()
    seen = set()
    temporary = None

    try:
        with BASELINE.open(
            encoding="utf-8-sig", newline=""
        ) as source:
            reader, header = reader_for(source)

            absent = sorted(REQUIRED_BASELINE - set(header))

            if absent:
                raise ValueError(
                    "Missing baseline columns; "
                    f"run prior stages first: {absent}"
                )

            columns = header + [
                column
                for column in new_columns
                if column not in header
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
                    target, fieldnames=columns
                )
                writer.writeheader()

                for row in reader:
                    eid = checked_eid(
                        row, reader.line_num, seen
                    )

                    baseline_date = parse_date(
                        row["baseline_date"]
                    )

                    if baseline_date is None:
                        raise ValueError(
                            "Invalid baseline_date on CSV line "
                            f"{reader.line_num}"
                        )

                    record = records.get(eid)

                    counts["baseline_rows"] += 1

                    match_status = (
                        "matched_rows"
                        if record is not None
                        else "unmatched_rows"
                    )
                    counts[match_status] += 1

                    derived = phenotype_variables(
                        record, baseline_date, row
                    )

                    for name in (
                        "clinical_cvd",
                        "hypertension_status",
                        "diabetes_status",
                        "prediabetes",
                    ):
                        value = derived[name]
                        label = "NA" if value is None else value
                        counts[f"{name}:{label}"] += 1

                    row.update({
                        name: "NA" if value is None else value
                        for name, value in derived.items()
                    })

                    writer.writerow(row)

        temporary.replace(BASELINE)
        print(dict(counts))

    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


if __name__ == "__main__":
    main()