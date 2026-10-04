#!/usr/bin/env python3
"""Convierte la salida de plink --freq [counts] --within (.frq.strat[.gz]) al formato de TreeMix."""

import argparse
import gzip
import sys
from collections import defaultdict

NA_VALUES = {"NA", "nan", "NaN", "."}


def open_text(path, mode):
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t")
    return open(path, mode)


def read_counts(path):
    """Devuelve (snps, pop2rs): snps en orden de aparición y {pop: {snp: (mc, total)}}."""
    snps = []
    seen = set()
    pop2rs = defaultdict(dict)

    with open_text(path, "r") as infile:
        header = next(infile, "").split()
        col = {name: i for i, name in enumerate(header)}

        missing = [c for c in ("SNP", "CLST", "NCHROBS") if c not in col]
        if missing:
            sys.exit(f"ERROR: faltan columnas {missing} en la cabecera: {header}")
        if "MAC" not in col and "MAF" not in col:
            sys.exit("ERROR: el archivo no tiene columna MAC ni MAF. "
                     f"Cabecera encontrada: {header}")

        use_mac = "MAC" in col
        if not use_mac:
            print("AVISO: no hay columna MAC; se calcula MAC = round(MAF x NCHROBS).",
                  file=sys.stderr)

        i_snp, i_pop, i_tot = col["SNP"], col["CLST"], col["NCHROBS"]
        i_val = col["MAC"] if use_mac else col["MAF"]
        ncols = len(header)

        for n, line in enumerate(infile, start=2):
            fields = line.split()
            if not fields:
                continue
            if len(fields) < ncols:
                sys.exit(f"ERROR: la línea {n} tiene {len(fields)} columnas, se esperaban {ncols}.")

            rs, pop = fields[i_snp], fields[i_pop]
            if rs not in seen:
                seen.add(rs)
                snps.append(rs)
            pops_entry = pop2rs[pop]  # asegura que la población exista aunque falten datos

            val, tot = fields[i_val], fields[i_tot]
            if val in NA_VALUES or tot in NA_VALUES or int(float(tot)) == 0:
                continue  # sin datos: se trata como faltante

            total = int(float(tot))
            mc = int(float(val)) if use_mac else int(round(float(val) * total))
            # se conserva la primera aparición de cada (pop, snp)
            pops_entry.setdefault(rs, (mc, total))

    return snps, pop2rs


def write_treemix(path, snps, pop2rs, on_missing="fail"):
    pops = list(pop2rs)  # orden de primera aparición
    n_imputed = 0

    with open_text(path, "w") as outfile:
        outfile.write(" ".join(pops) + "\n")
        for rs in snps:
            row = []
            for pop in pops:
                try:
                    mc, total = pop2rs[pop][rs]
                except KeyError:
                    if on_missing == "fail":
                        sys.exit(f"ERROR: el SNP '{rs}' no tiene datos en la población '{pop}'. "
                                 "Usa -e impute para rellenar con 0,0.")
                    mc, total = 0, 0
                    n_imputed += 1
                row.append(f"{mc},{total - mc}")
            outfile.write(" ".join(row) + "\n")

    if n_imputed:
        print(f"ADVERTENCIA: {n_imputed} valores faltantes imputados como 0,0", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Convierte la salida de PLINK (.frq.strat[.gz]) a formato TreeMix."
    )
    parser.add_argument("-i", "--input", required=True, help="archivo de entrada (gzip o texto)")
    parser.add_argument("-o", "--output", required=True,
                        help="archivo de salida (gzip si termina en .gz)")
    parser.add_argument(
        "-e", "--on-missing", choices=["fail", "impute"], default="fail",
        help="qué hacer si un SNP falta o no tiene datos en alguna población: "
             "fail = error y salir (por defecto), impute = rellenar con 0,0",
    )
    args = parser.parse_args()

    snps, pop2rs = read_counts(args.input)
    write_treemix(args.output, snps, pop2rs, args.on_missing)


if __name__ == "__main__":
    main()
