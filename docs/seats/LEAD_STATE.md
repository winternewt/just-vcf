**just-vcf: where it stands** (updated 2026-11-23)

just-vcf reads 23andMe-style array exports, other genotype tables and VCF/gVCF into one set of tables, and writes VCF from them. Version 0.1.0 is ready; releasing it is the person's call. Most open plan items wait on answers on the person's desk.

**Headline results**
- ✅ **More 23andMe insertion/deletion calls can now be read.** About 3,900 of 23,000 such rows on the v3 chip (17%) and 2,200 of 144,000 on v4 (1.5%). This was learned from 36 public PGP volunteers whose array exports could be checked against their whole genomes. The rest stay unread, each with its reason.
- ⚠️ **How often two calling methods disagree depends on the person as well as the method.** Nine PGP volunteers were each sequenced twice by Complete Genomics, with the same two pipeline versions. Their disagreement on single-letter variants in easy regions ranges from about 0.1% to 1.8%. One pooled rate per method would hide that spread.
- ✅ **The error atlas has all its planned pairs.** The last three PGP volunteers are in: DRAGEN agrees with Complete Genomics on 99.93% of single-letter variants in easy regions, while an Illumina array disagrees with the same genome on 2% of them, unexplained so far. A fault that made array and SNV-only files look wrong on every insertion was found on the way and fixed, in the comparison and in the combined call alike.
- ✅ **Runs on polars 2.** It now uses polars-bio 0.37 and polars 2.1, and it still works with polars 1.44 for projects that haven't moved yet. The full test suite, 2,889 tests, passes on both versions.

**Next**
1. ⏳ Release 0.1.0 — the person.
2. 🚧 The same fix for ancestry and relatedness (RM302, building now), then arrays read probe by probe (RM303) — lead, unattended. Plain-text array exports are covered since RM301.
3. 🚧 Last rules for imputed data (RM287) — lead, waits on desk #212, #213, #214.
4. ⏳ The remaining v4 insertion/deletion rows (RM281) — lead, looked at again 2026-12-07.
