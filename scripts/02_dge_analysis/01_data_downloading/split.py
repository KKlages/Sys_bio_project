import csv
import os
import re

with open('download_links.csv', newline='') as infile:
    reader = csv.DictReader(infile)
    batches = {}
    for row in reader:
        abbr = row['Abbrev'].strip()
        # Clean up abbreviation for a safe filename
        safe_abbr = re.sub(r'\W+', '_', abbr)
        if safe_abbr not in batches:
            batches[safe_abbr] = []
        batches[safe_abbr].append(row)

    for abbr, rows in batches.items():
        fname = f'batch_{abbr}.csv'
        with open(fname, 'w', newline='') as outfile:
            writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames)
            writer.writeheader()
            writer.writerows(rows)

