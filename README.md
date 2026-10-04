# py3plink2treemix

Convert PLINK allele-frequency output stratified by population (`.frq.strat`) into the input format required by [TreeMix](https://bitbucket.org/nygcresearch/treemix).

Python 3 port and refactor of the original `plink2treemix.py` (Python 2).

## Requirements

- Python 3.6+ (standard library only)
- PLINK 1.07 or 1.9 to generate the input file

## Workflow: VCF to TreeMix

1. **Prepare a population file** (`clusters.txt`) with three whitespace-separated columns (spaces or tabs): `FID` `IID` `POPULATION`.

```
   0 sample_A 1
   0 sample_B 1
   0 sample_C 2
```

2. **Compute per-population allele counts with PLINK:**

```bash
   plink --vcf data.vcf --allow-extra-chr --const-fid 0 \
         --within clusters.txt --freq gz --out countfreq
```

   - `--const-fid 0` makes each sample `FID=0`, `IID=<full sample name>`. If your `clusters.txt` has `FID = IID`, use `--double-id` instead.
   - `--set-missing-var-ids @:#` gives each SNP a unique ID (`chr:pos`) when the VCF ID column is `.`. Without unique IDs, SNPs get merged.
   - `--freq counts` writes the `MAC` column (recommended). Without `counts`, PLINK writes `MAF` and the script derives counts from it.
   - If PLINK warns `No samples named in --within file remain in the current analysis`, the FID/IID in `clusters.txt` do not match the VCF sample names, and no `.frq.strat.gz` file is produced.

3. **Run the conversion:**

```bash
   python3 py3plink2treemix.py -i freq1.frq.strat.gz -o treemix_in.gz
```

4. **Run TreeMix** (see notes below):

```bash
   treemix -i treemix_in.gz -m 1 -k 500 -noss -o tm_m1
```

## Usage

```
py3plink2treemix.py -i INPUT -o OUTPUT [-e {fail,impute}]
```

| Option | Description |
|---|---|
| `-i`, `--input` | Input `.frq.strat` file (gzip or plain text). Required. |
| `-o`, `--output` | Output file for TreeMix. Gzipped if the name ends in `.gz`. Required. |
| `-e`, `--on-missing` | What to do when a SNP has no data in some population: `fail` (default) stops with an error; `impute` fills with `0,0` and prints a warning to stderr. |

## Input

A PLINK `.frq.strat` file. Columns are detected by header name, so both PLINK versions work. Required columns: `SNP`, `CLST`, `NCHROBS`, and either `MAC` or `MAF`.

- If `MAC` is present, it is used as is.
- If only `MAF` is present, the count is computed as `round(MAF * NCHROBS)` and a notice is printed to stderr.
- Values of `NA`, or `NCHROBS = 0`, are treated as missing data for that SNP/population.
- The counted allele is `A1`, as defined by PLINK.

## Output

A TreeMix input file:

- Line 1: population names, space-separated, in order of first appearance in the input.
- One line per SNP, in input order, with one `count_A1,count_other` pair per population.

```
1 2 3
3,1 5,11 2,10
0,4 8,8 6,6
```

## Notes

- **SNP order matters.** TreeMix builds blocks (`-k`) from consecutive SNPs in file order, so the VCF must be sorted by position.
- **Use polymorphic sites only.** Invariant sites add no information.
- **Choosing `-k`.** Pick a block size that covers the scale of linkage disequilibrium while leaving at least ~100 blocks (e.g. `-k 500` for ~150,000 SNPs). With only a few hundred SNPs (proof of concept), use a small value such as `-k 5`.
- **Small populations.** Populations with very few samples give noisy allele frequencies; consider `-noss` in TreeMix and merging or enlarging such groups when possible.
- **Number of migrations (`-m`).** Identifiable migration edges are limited by the number of populations. With 3 populations, `-m` above 1-2 is not meaningful.
- **No outgroup.** `-root` is optional in TreeMix; without it the tree is rooted arbitrarily, so migration directions should not be over-interpreted.

## Differences from the original script

- Python 3 (print and `gzip` text mode, no `has_key`).
- `-i` / `-o` options via `argparse` instead of positional arguments.
- Columns located by header name instead of fixed positions.
- Supports `MAC` or `MAF`, and gzipped or plain input/output.
- Explicit handling of missing SNP/population combinations (`-e`).
- Output no longer has a trailing space at the end of each line.
