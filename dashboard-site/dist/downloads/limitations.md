# Limitations and interpretation guardrails

1. **This V1 is single-instrument MR.** It does not have a complete,
   ancestry-matched regional LD matrix for every gene and therefore does not
   claim MR-link-2, generalized IVW, SuSiE-MR, heterogeneity, or pleiotropy-
   robust estimation genome-wide.

   Genome-wide gene-to-CAD estimates also lack regional colocalization. A
   single-SNP Wald ratio cannot distinguish expression mediation from linkage
   or horizontal pleiotropy, so this layer is ratio evidence rather than
   confirmed causality.

2. **Whole blood is a discovery context, not the whole biology of CAD.** It is
   scalable and well powered, but LDL metabolism has central liver biology and
   CAD also involves arterial and immune contexts. Only the SORT1 region has a
   liver audit in this release.

3. **The gene-to-gene layer has sample overlap.** Source and target eQTLGen
   effects were estimated in overlapping cohorts. Their covariance is not
   released and is set to zero, so the delta standard errors and spike-and-slab
   edge probabilities can be overconfident or otherwise miscalibrated.

4. **A trans association does not identify mediation.** A variant can affect
   the nominated source and target independently. Shared upstream regulation,
   linkage, horizontal pleiotropy, or an unmeasured nearby source gene can all
   mimic a source-expression → target-expression relationship.

5. **The trans-tested SNP universe is selected.** eQTLGen Phase I tested about
   10,000 known trait-associated variants in trans, not every common variant.
   The network is conditional on that universe and cannot be compared directly
   with a genome-wide exhaustive trans scan. Trait-associated selection also
   makes broad horizontal pleiotropy especially plausible.

6. **Public eQTLGen effect sizes are reconstructed approximately.** Z/sqrt(N)
   and 1/sqrt(N) are genotype-standardized approximations, not native
   expression beta/SE values. Real CAD allele frequencies are used to put its
   per-allele beta/SE on the same genotype-SD scale, but the approximation can
   still affect calibration.

7. **One top instrument invites winner's curse.** The significant cis table is
   used to select instruments. The gene-to-CAD outcome is independent, but
   exposure magnitudes can still be overestimated.

8. **Physical source-locus collapse is not LD inference.** The 1 Mb rule avoids
   counting nearby genes repeatedly but does not prove source loci are
   independent. Within-locus probabilities use a conservative maximum because
   the edges are correlated; this is not a joint locus posterior. Across-locus
   independence in the Poisson-binomial model also remains unverified.

9. **Candidate ascertainment is restrictive.** Results apply to genes surfaced
   at genome-wide significant LDL loci by Open Targets L2G ≥0.40 plus ten
   BayesV2G seeds. A gene absent from the table was not shown to be non-core.
   The ranking is therefore “core-like among LDL-locus, whole-blood-testable
   candidates,” not a genome-wide list of general CAD core genes.

10. **Core probabilities are conditional model outputs.** They depend on the
    point-normal priors, independence assumptions, tested-pair coverage, and
    operational thresholds. They do not quantify the chance that a
    philosophical or universally accepted “core gene” label is true.
    With many source loci, the edge prior alone gives nonzero probability of
    convergence; the report therefore shows the prior-only baseline.

11. **Evidence channels are not multiplied indiscriminately.** The liver
    colocalization result is kept as a separate locus audit because MR and
    colocalization reuse the same association statistics. Treating them as
    independent Bayes factors would double-count evidence.

12. **Classic coloc is single-signal.** Multiple causal variants or inaccurate
    regional priors can change PP(H4). A production extension should use
    ancestry-matched LD, fine mapping, and SuSiE-coloc.

13. **The annotation audit has no true negative class.** PanelApp is a monogenic
    diagnostic proxy. Retrieval/enrichment can support face validity but cannot
    calibrate `P(Core)` or establish that non-panel genes are false positives.

14. **Lifelong genetic effects are not interventions.** An MR effect of
    genetically predicted expression is not automatically the effect of an
    acute drug, a particular dose, or perturbation in a specific cell state.

This project should be described as a complete, auditable MVP and feasibility
analysis—not a clinically validated gene-discovery system.
