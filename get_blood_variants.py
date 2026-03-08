from get_tumor_variants import patients

from cyvcf2 import VCF

patient_to_blood_vcf_path = {p: f"/path/to/data_analysis/numbat/bcftools/{p}/{p}_variants.vcf.gz" for p in patients}

patient_to_blood_variants = {
    p: {
        (rec.CHROM, rec.POS, rec.REF, alt)
        for rec in VCF(path)
        if rec.CHROM != "chrM" and (rec.FILTER is None or str(rec.FILTER).lower() == "pass")
        for alt in rec.ALT
    }
    for p, path in patient_to_blood_vcf_path.items()
}

print({p: len(s) for p, s in patient_to_blood_variants.items()})