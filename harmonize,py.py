import pandas as pd

DARWIN_CORE_FIELDS = {
    "scientificName": ["species", "scientific_name", "taxon", "name"],
    "decimalLatitude": ["latitude", "lat", "y"],
    "decimalLongitude": ["longitude", "lon", "lng", "x"],
    "eventDate": ["date", "observed_on", "recorded_date"],
    "individualCount": ["abundance", "count", "n", "individuals"],
    "locality": ["site_name", "location", "place"],
}

def map_to_darwin_core(df):
    rename_map = {}
    lower_cols = {c.lower().strip(): c for c in df.columns}
    for dwc_field, aliases in DARWIN_CORE_FIELDS.items():
        for alias in aliases:
            if alias in lower_cols:
                rename_map[lower_cols[alias]] = dwc_field
                break
    return df.rename(columns=rename_map)

def flag_quality_issues(df):
    flags = []
    for _, row in df.iterrows():
        issues = []
        lat, lon = row.get("decimalLatitude"), row.get("decimalLongitude")
        if pd.isna(lat) or pd.isna(lon):
            issues.append("missing coordinates")
        elif not (-90 <= lat <= 90 and -180 <= lon <= 180):
            issues.append("implausible coordinates")
        if pd.isna(row.get("scientificName")):
            issues.append("missing species name")
        flags.append("; ".join(issues) if issues else "OK")
    df["quality_flag"] = flags
    return df

def harmonization_report(df):
    total = len(df)
    flagged = (df["quality_flag"] != "OK").sum()
    return {
        "total_records": total,
        "flagged_records": int(flagged),
        "clean_records": total - int(flagged),
    }

if __name__ == "__main__":
    df = pd.read_csv("gbif_test.csv", sep="\t", low_memory=False)
    df = map_to_darwin_core(df)
    df = flag_quality_issues(df)
    print(harmonization_report(df))