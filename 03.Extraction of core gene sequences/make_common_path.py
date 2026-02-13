#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
make_common_path.py

Generate a single-sample "common path" VCF from a multi-sample pangenome VCF.
Selects the major allele (present in most samples) when frequency > MAJ_FRAC
and |SVLEN| >= MIN_LEN.

Rules:
- Allele is considered "present" if any GT > 0
- For multi-allelic sites: choose ALT with the highest number of carriers
- Only apply ALT if presence frequency > MAJ_FRAC and |SVLEN| >= MIN_LEN
- For symbolic ALTs (<INS>, <DEL>, ...), try to recover sequence from INFO/SEQ
- Outputs provenance TSV with origin and carrier information

Usage:
    ./make_common_path.py input.vcf.gz ref.fa output.common.vcf.gz provenance.tsv [MIN_LEN] [MAJ_FRAC]

Example:
    ./make_common_path.py Os.norm.sorted.vcf.gz Os_MSU.fa Os.common.vcf.gz Os.provenance.tsv 500 0.5
"""

import sys
import os
from collections import defaultdict
import pysam

try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False


def alt_is_sequence(alt_obj) -> bool:
    s = str(alt_obj)
    return not (s.startswith("<") and s.endswith(">"))


def infer_svlen(rec: pysam.VariantRecord, alt_idx: int) -> int:
    fallback = 0
    try:
        ref_len = len(rec.ref or "")
        alt_len = len(str(rec.alts[alt_idx]) or "")
        fallback = alt_len - ref_len
    except Exception:
        pass

    try:
        svlen = rec.info.get("SVLEN")
        if isinstance(svlen, (list, tuple)):
            val = svlen[alt_idx]
            if isinstance(val, int) and val != 0:
                return val
        elif isinstance(svlen, int) and svlen != 0:
            return svlen
    except Exception:
        pass
    return fallback


def sample_has_alt(gt_tuple):
    if gt_tuple is None:
        return False, True
    if any(a is None for a in gt_tuple):
        return False, True
    if any(a is not None and a > 0 for a in gt_tuple):
        return True, False
    return False, False


def choose_major_alt(rec: pysam.VariantRecord, sample_names):
    if not rec.alts:
        return None, set(), 0

    per_alt_carriers = {i: set() for i in range(len(rec.alts))}
    called = 0

    for s in sample_names:
        sm = rec.samples[s]
        gt = sm.get("GT", None)
        has_alt, missing = sample_has_alt(gt)
        if missing:
            continue
        called += 1
        if not has_alt:
            continue
        alts_present = {a for a in gt if a and a > 0}
        for a in alts_present:
            idx = a - 1
            if 0 <= idx < len(rec.alts):
                per_alt_carriers[idx].add(s)

    best_idx = max(per_alt_carriers, key=lambda i: len(per_alt_carriers[i]), default=None)
    if best_idx is None:
        return None, set(), called

    return best_idx, per_alt_carriers[best_idx], called


def header_has_info(hdr, key: str) -> bool:
    try:
        return key in hdr.info
    except Exception:
        return False


def header_has_format(hdr, key: str) -> bool:
    try:
        return key in hdr.formats
    except Exception:
        return False


def rebuild_single_sample_header(src_header: pysam.VariantHeader, sample_name="COMMON") -> pysam.VariantHeader:
    new_header = pysam.VariantHeader()

    # Copy contigs
    for ctg in src_header.contigs.values():
        new_header.contigs.add(ctg.name, length=ctg.length)

    # Copy INFO fields
    for key, info in src_header.info.items():
        try:
            new_header.info.add(key, number=info.number, type=info.type,
                                description=info.description or "")
        except Exception:
            pass

    # Copy FORMAT fields
    for key, fmt in src_header.formats.items():
        try:
            new_header.formats.add(key, number=fmt.number, type=fmt.type,
                                   description=fmt.description or "")
        except Exception:
            pass

    # Ensure GT exists
    if "GT" not in new_header.formats:
        new_header.formats.add("GT", number=1, type="String", description="Genotype")

    new_header.add_sample(sample_name)
    return new_header


def main(in_vcf, ref_fa, out_common_vcf, out_prov_tsv, min_abs_len=500, majority=0.5):
    invcf = pysam.VariantFile(in_vcf)
    sample_names = list(invcf.header.samples)
    if not sample_names:
        raise ValueError("Input VCF has no sample columns. Please provide a multi-sample VCF.")

    header = rebuild_single_sample_header(invcf.header, sample_name="COMMON")
    outvcf = pysam.VariantFile(out_common_vcf, "wz", header=header)

    with open(out_prov_tsv, "w") as prov:
        prov.write("\t".join([
            "CHROM", "POS", "END", "SVTYPE", "SVLEN", "AF_presence",
            "N_carriers", "N_called", "ALT_repr", "SAMPLES_WITH_ALT", "SOURCE_VCF"
        ]) + "\n")

        iterator = invcf
        if HAS_TQDM:
            iterator = tqdm(iterator, desc="Processing variants", unit="var")

        for rec in iterator:
            chrom = rec.chrom
            pos = rec.pos
            stop = rec.stop if hasattr(rec, "stop") else (rec.start + len(rec.ref or "") - 1)
            ref = rec.ref
            alts = rec.alts or ()

            if not alts:
                new = header.new_record(contig=chrom, start=rec.start, stop=stop,
                                        id=rec.id, qual=rec.qual, filter=rec.filter.keys())
                new.ref = ref
                new.alts = alts
                new.samples["COMMON"]["GT"] = (0, 0)
                outvcf.write(new)
                continue

            chosen_idx, carriers, called = choose_major_alt(rec, sample_names)
            if chosen_idx is None or called == 0:
                new = header.new_record(contig=chrom, start=rec.start, stop=stop,
                                        id=rec.id, qual=rec.qual, filter=rec.filter.keys())
                new.ref = ref
                new.alts = alts
                new.samples["COMMON"]["GT"] = (0, 0)
                outvcf.write(new)
                continue

            af_pres = len(carriers) / called if called > 0 else 0.0
            svlen = infer_svlen(rec, chosen_idx)
            is_large = abs(svlen) >= min_abs_len

            svtype = rec.info.get("SVTYPE", None)
            if not svtype:
                svtype = "INS" if svlen > 0 else ("DEL" if svlen < 0 else "VAR")

            take_alt = (af_pres > majority) and is_large

            chosen_alt = alts[chosen_idx]
            alt_seq = str(chosen_alt)
            has_seq = alt_is_sequence(chosen_alt)

            if take_alt and not has_seq:
                info_seq = rec.info.get("SEQ", None)
                if isinstance(info_seq, (bytes, bytearray)):
                    info_seq = info_seq.decode("utf-8", "ignore")
                if isinstance(info_seq, str) and info_seq:
                    alt_seq = info_seq
                    has_seq = True
                else:
                    take_alt = False

            new = header.new_record(contig=chrom, start=rec.start, stop=stop,
                                    id=rec.id, qual=rec.qual, filter=rec.filter.keys())
            new.ref = ref

            if take_alt:
                new.alts = (alt_seq,)
                if header_has_info(header, "SVTYPE"):
                    new.info["SVTYPE"] = svtype
                if header_has_info(header, "SVLEN"):
                    new.info["SVLEN"] = len(alt_seq) - len(ref or "")
                new.samples["COMMON"]["GT"] = (1, 1)

                alt_repr = alt_seq if has_seq and len(alt_seq) <= 50 else ("<SEQ>" if has_seq else "<SYMBOLIC>")
                carriers_str = ",".join(sorted(carriers)) if carriers else "."
                prov.write("\t".join(map(str, [
                    chrom, pos, stop, svtype, svlen,
                    f"{af_pres:.6f}", len(carriers), called,
                    alt_repr, carriers_str, os.path.abspath(in_vcf)
                ])) + "\n")
            else:
                new.alts = alts
                new.samples["COMMON"]["GT"] = (0, 0)

            outvcf.write(new)

    outvcf.close()


if __name__ == "__main__":
    if len(sys.argv) < 5:
        print(f"Usage: {sys.argv[0]} in.vcf.gz ref.fa out_common.vcf.gz provenance.tsv [MIN_LEN=500] [MAJ_FRAC=0.5]")
        sys.exit(1)

    in_vcf = sys.argv[1]
    ref_fa = sys.argv[2]
    out_vcf = sys.argv[3]
    out_tsv = sys.argv[4]
    min_len = int(sys.argv[5]) if len(sys.argv) > 5 else 500
    maj_frac = float(sys.argv[6]) if len(sys.argv) > 6 else 0.5

    # Optional: ensure index exists
    try:
        pysam.faidx(ref_fa)
    except Exception:
        pass

    main(in_vcf, ref_fa, out_vcf, out_tsv, min_len, maj_frac)
