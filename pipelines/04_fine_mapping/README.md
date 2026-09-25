# 04 Fine-Mapping

Identification of independent signals at the genome-wide significant loci and
construction of credible sets.

> **Status:** scripts to be added.

## Expected contents

- Locus definition: window size around each sentinel and the rule for merging
  overlapping loci.
- Conditional and joint analysis (for example GCTA-COJO) and the LD reference
  used.
- Credible-set construction (for example Wakefield approximate Bayes factors,
  SuSiE or FINEMAP), including the prior and the credible-set coverage (for
  example 95%).

## Output used downstream

- Sentinel variants: `sentinel_file` in the V2G config.
- Credible-set variants with posterior probabilities: `credible_set_file` in
  the V2G config.
