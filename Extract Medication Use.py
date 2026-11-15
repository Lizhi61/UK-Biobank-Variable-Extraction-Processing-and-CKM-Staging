#!/usr/bin/env python3
"""Add baseline medication indicators to UKB_baseline.csv."""

import csv
from pathlib import Path
from tempfile import NamedTemporaryFile
from collections import Counter
from decimal import Decimal, InvalidOperation
import re

VISITS = Path("visit_data.csv")
BASELINE = Path("UKB_baseline.csv")
MISSING = {"", "na", "n/a", "nan", "null", "none"}

# Field 20003 medication codes supplied for this study.
STATINS = {
    1141146234, 1140910632, 1140888594, 1140888648, 1141192410,
    1140861958, 1141192414, 1140864592, 1141146138, 1140861970,
    1141188146, 1140881748, 1141200040,
}
BETA_BLOCKERS = {
    1140866724, 1141164276, 1140866738, 1140879830, 1140879760,
    1140860292, 1140909368, 1140879842, 1140879762, 1140879854,
    1140879824, 1140866704, 1140879818, 1140866764, 1140860192,
    1140851556, 1141182904, 1141182968, 1140860362, 1140860250,
    1140866756, 1140860222, 1140866758, 1140860418, 1140860332,
    1140860304, 1141168498, 1140866782, 1140866802, 1140860266,
    1141152076, 1141171152, 1141156754, 1140866712, 1140866800,
    1140860498, 1140866804, 1140860194, 1140860274, 1140860492,
    1140860278, 1140860294, 1141146124, 1141146126, 1141146128,
    1141194810, 1140860426, 1141180778, 1140864950, 1140866692,
    1140866726, 1140860434, 1140860220, 1141164280, 1140860308,
    1140851484, 1140860404, 1140916868, 1140860312, 1140866766,
    1141194804, 1141187048,
}
CALCIUM_CHANNEL_BLOCKERS = {
    1140879802, 1140879810, 1140879806, 1140861088, 1140888646,
    1140872568, 1140861190, 1140888510, 1140861276, 1140861090,
    1141153026, 1140881702, 1141188920, 1140872472, 1140861282,
    1141201814, 1141200782, 1140928212, 1141157140, 1140866466,
    1141150538, 1141150500, 1141169730, 1140911698, 1140923572,
    1140861176, 1140861138, 1141199858, 1140926780, 1140861120,
    1140917428, 1140866554, 1141187094, 1140861166, 1141173766,
    1141157136, 1140927940, 1141150926, 1140861128, 1141169710,
    1140926188, 1141153032, 1140881692, 1140926778, 1141190160,
    1141165470, 1141187774, 1140866460, 1141188576, 1141188936,
    1141188152, 1140861202, 1141145870,
}
RENIN_ANGIOTENSIN_CODES = {
    1141156836, 1140923712, 1140860750, 1141193282, 1140860882,
    1140888560, 1140888552, 1140860728, 1141171336, 1140860806,
    1140888556, 1141166006, 1141164148, 1140860904, 1141145660,
    1141152998, 1140881706, 1140860696, 1141156846, 1140916356,
    1141153006, 1140860758, 1141193346, 1140864910, 1140881712,
    1140860706, 1140860878, 1141179974, 1141164154, 1140916362,
    1141171344, 1141145668, 1141188408, 1140860912, 1140860892,
    1140860776, 1140860714, 1141150560, 1140860764, 1141199940,
    1141151018, 1141172492, 1140860790, 1140864952, 1141153328,
    1141151016, 1141201038, 1141180592, 1141145658, 1141153316,
    1141150898, 1141187788,
}
# Indapamide is a diuretic. Keep it in the supplied antihypertensive phenotype.
ANTIHYPERTENSIVE_CODES = (
    BETA_BLOCKERS | CALCIUM_CHANNEL_BLOCKERS | RENIN_ANGIOTENSIN_CODES | {1140866078}
)
GLUCOSE_LOWERING_CODES = {
    1140883066, 1140874690, 1140884600, 1140874706, 1140874646,
    1140874712, 1140874650, 1140874716, 1140874658, 1140874718,
    1140874660, 1140874724, 1141173882, 1140868902, 1140874664,
    1140874726, 1140874666, 1140874728, 1140874674, 1140874736,
    1140874678, 1140874744, 1140874680, 1140874746, 1140874686,
    1141177600, 1140874646, 1141189090, 1141157284, 1141171646,
    1141153254, 1141152590,
}

# Each rule maps an output column to its touchscreen option and interview codes.
MEDICATION_RULES = (
    ("antihypertensive_medication_use", 2, ANTIHYPERTENSIVE_CODES),
    ("lipid_lowering_medication_use", 1, STATINS),
    ("glucose_lowering_medication_use", 3, GLUCOSE_LOWERING_CODES),
)


def integer_code(raw):
    """Parse an exact finite integer code without floating-point rounding."""
    raw = "" if raw is None else str(raw).strip()

    if raw.lower() in MISSING:
        return None

    try:
        value = Decimal(raw)
    except InvalidOperation:
        return None

    if not value.is_finite() or value != value.to_integral_value():
        return None

    return int(value)


def normalize_column(column):
    column = column.strip()

    if column in {"eid", "f.eid"}:
        return "eid"

    for pattern in (
        r"(\d+)-(\d+)\.(\d+)",
        r"f\.(\d+)\.(\d+)\.(\d+)",
        r"p(\d+)_i(\d+)(?:_a(\d+))?",
    ):
        match = re.fullmatch(pattern, column)

        if match:
            field, instance, array = match.groups()
            return f"{int(field)}-{int(instance)}.{int(array or 0)}"

    return column


def validated_reader(source):
    reader = csv.DictReader(source, strict=True)

    if not reader.fieldnames:
        raise ValueError("The input file is empty or has no header")

    header = [normalize_column(c) for c in reader.fieldnames]
    duplicates = [c for c, n in Counter(header).items() if n > 1]

    if "" in header or duplicates:
        raise ValueError(
            f"Empty or duplicate/alias-conflicting headers: {duplicates}"
        )

    if "eid" not in header:
        raise ValueError("Missing participant ID column: eid or f.eid")

    reader.fieldnames = header
    return reader, header


def checked_eid(row, line_number, seen):
    if None in row or any(value is None for value in row.values()):
        raise ValueError(
            f"CSV line {line_number} has an incorrect number of cells"
        )

    eid = row["eid"].strip()

    if not re.fullmatch(r"[0-9]+", eid):
        raise ValueError(f"CSV line {line_number} has an invalid eid")

    if eid in seen:
        raise ValueError(f"CSV line {line_number} has a duplicate eid")

    seen.add(eid)
    return eid


def recorded_codes(row, columns):
    codes = set()

    for column in columns:
        value = row[column].strip()

        if value.lower() not in MISSING:
            # Preserve unparseable entries rather than treating them as empty.
            codes.add(integer_code(value))

    return codes


def medication_flags(questionnaire, interview, medication_count=None):
    """Use conservative three-state rules with positive evidence first."""
    coded = {
        code
        for code in interview
        if code is not None and code > 0
    }

    valid_questionnaire = bool(questionnaire) and (
        questionnaire == {-7}
        or questionnaire <= {1, 2, 3, 4, 5}
    )

    # This is a recorded-zero rule, not proof of complete ascertainment.
    recorded_zero = medication_count == 0 and not interview

    flags = {}

    for name, option, drug_codes in MEDICATION_RULES:
        if option in questionnaire or coded & drug_codes:
            flags[name] = 1

        elif recorded_zero:
            flags[name] = 0

        elif option != 3 and valid_questionnaire:
            flags[name] = 0

        else:
            # No insulin response does not exclude non-insulin drugs.
            # Unmatched interview entries cannot establish absence.
            flags[name] = "NA"

    return flags


def load_medications():
    with VISITS.open(encoding="utf-8-sig", newline="") as source:
        reader, header = validated_reader(source)
        fields = {}

        # Require complete baseline multi-select slots for this input contract.
        for field, size in ((6153, 4), (6177, 3), (20003, 48)):
            fields[field] = [
                f"{field}-0.{index}" for index in range(size)
            ]

            missing = sorted(set(fields[field]) - set(header))

            if missing:
                raise ValueError(
                    f"Missing baseline array columns: {missing}"
                )

        if "137-0.0" not in header:
            print(
                "Warning: field 137 is absent; "
                "blank interviews remain unknown"
            )

        results = {}
        seen = set()
        touchscreen = fields[6153] + fields[6177]

        for row in reader:
            eid = checked_eid(row, reader.line_num, seen)

            results[eid] = medication_flags(
                recorded_codes(row, touchscreen),
                recorded_codes(row, fields[20003]),
                integer_code(row.get("137-0.0")),
            )

    return results


def main():
    medications = load_medications()
    missing_flags = {
        name: "NA" for name, _, _ in MEDICATION_RULES
    }

    counts = Counter()
    seen = set()
    temporary = None

    try:
        with BASELINE.open(encoding="utf-8-sig", newline="") as source:
            reader, header = validated_reader(source)

            columns = header + [
                name for name in missing_flags if name not in header
            ]

            with NamedTemporaryFile(
                mode="w",
                encoding="utf-8",
                newline="",
                dir=BASELINE.parent,
                delete=False,
            ) as target:
                temporary = Path(target.name)
                writer = csv.DictWriter(target, fieldnames=columns)
                writer.writeheader()

                for row in reader:
                    eid = checked_eid(row, reader.line_num, seen)
                    counts["baseline_rows"] += 1

                    match_status = (
                        "matched_rows"
                        if eid in medications
                        else "unmatched_rows"
                    )
                    counts[match_status] += 1

                    flags = medications.get(eid, missing_flags)
                    row.update(flags)
                    writer.writerow(row)

                    for name, value in flags.items():
                        counts[f"{name}:{value}"] += 1

        temporary.replace(BASELINE)
        print(dict(counts))

    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


if __name__ == "__main__":
    main()