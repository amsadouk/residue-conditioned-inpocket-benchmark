# Reference list: structure-based generative design, and the benchmarks we compare against

Companion to `octara_references.md` (methods and software) and `octara_intro_and_theory.md`.

Every number below was taken from the original paper and the table it appears in. Where a value
could not be verified it is marked UNVERIFIED rather than filled in. Two conventions in this
literature are incompatible and are flagged at each use.

---

## The benchmark everyone reports on

1. Francoeur PG, Masuda T, Sunseri J, Jia A, Iovanisci RB, Snyder I, Koes DR. Three-Dimensional
   Convolutional Neural Networks and a Cross-Docked Data Set for Structure-Based Drug Design.
   *J Chem Inf Model* 2020;60(9):4200-4215.
   *CrossDocked2020. Approximately 18,450 bound crystal structures, clustered by pocket
   similarity with ProBiS, each ligand redocked to its cognate receptor and cross-docked to every
   receptor with a similar pocket. The standard test split is 100 pockets, obtained by filtering
   to RMSD > 1 A and clustering at 30% sequence identity with mmseqs2.*

**What the "reference ligand" in that benchmark actually is.** It is the ligand that happened to
be co-crystallised with the pocket. There is no potency selection of any kind. It may be a
fragment, a tool compound or a weak binder. The reference set averages **Vina -6.36 kcal/mol**
(verified in Luo Table 1, TargetDiff Table 1 and MolCRAFT Table 2, which agree). Consequently
"beats the reference" is a far weaker statement than it sounds, and it is not comparable to a
comparison against a target's most potent known actives.

---

## The generators

Vina values in kcal/mol. "HA %" is the reported high-affinity fraction: the percentage of
generated molecules scoring better than the reference ligand.

2. Luo S, Guan J, Ma J, Peng J. A 3D Generative Model for Structure-Based Drug Design.
   *NeurIPS* 2021. *(AR / 3D-SBDD.)* Table 1: Vina -6.344 mean / -6.200 median, HA 29.09 / 18.50,
   QED 0.525, SA 0.657, diversity 0.720.
3. Peng X, Luo S, Guan J, Xie Q, Peng J, Ma J. Pocket2Mol: Efficient Molecular Sampling Based on
   3D Protein Pockets. *ICML* 2022. Table 1, scored with QVina: -7.288 +/- 2.53, QED 0.563,
   SA 0.765, diversity 0.688.
4. Guan J, Qian WW, Peng X, Su Y, Peng J, Ma J. 3D Equivariant Diffusion for Target-Aware
   Molecule Generation and Affinity Prediction. *ICLR* 2023. *(TargetDiff.)* Table 1:
   Vina Score -5.47 / -6.30, Vina Min -6.64 / -6.83, Vina Dock -7.80 / -7.91, HA 58.1 / 59.1,
   QED 0.48, SA 0.58, diversity 0.72, 24.2 heavy atoms.
5. Guan J, Zhou X, Yang Y, et al. DecompDiff: Diffusion Models with Decomposed Priors for
   Structure-Based Drug Design. *ICML* 2023. Table 3: Vina Score -5.67 / -6.04,
   Vina Dock -8.39 / -8.43, HA 64.4, QED 0.45, SA 0.61, 29.4 heavy atoms.
6. Qu Y, Qiu K, Song Y, et al. MolCRAFT: Structure-Based Drug Design in Continuous Parameter
   Space. *ICML* 2024. Table 2: Vina Score -6.59 / -7.04, Vina Min -7.27 / -7.26,
   Vina Dock -7.92 / -8.01, QED 0.50, SA 0.69, diversity 0.72, 22.7 heavy atoms.
7. Ragoza M, Masuda T, Koes DR. Generating 3D molecules conditional on receptor binding sites
   with deep generative models. *Chem Sci* 2022. *(liGAN.)* **Evaluated on 10 CrossDocked
   targets, not the 100-pocket split**; 98.5% validity (posterior), 100% novelty, HA 30.8%
   (posterior) / 17.3% (prior). Any "liGAN on 100 pockets" figure is a third-party
   re-evaluation, not this paper.
8. Schneuing A, Du Y, Harris C, et al. Structure-based drug design with equivariant diffusion
   models. *(DiffSBDD.)* **Reports distributional comparisons, not absolute Vina/QED/SA tables**
   for CrossDocked. A widely circulated "DiffSBDD Table 1" of QED 0.106 / Vina 0.596 contains
   divergence values, not scores, and must not be quoted as performance.
9. Zhang O, Zhang J, Jin J, et al. ResGen: a pocket-aware 3D molecular generation model.
   *Nat Mach Intell* 2023;5:1020-1030. Full text was not accessible. The only figure obtained is
   second-hand via GenBench3D: top-10 mean Vina -9.08 after redocking, against Pocket2Mol -8.73.
   That is a post-redocking number and is **not** comparable to in-place Vina Score.
   All other ResGen metrics: UNVERIFIED.

**A trap in this table.** The liGAN and AR rows usually quoted in later papers are TargetDiff's
re-evaluation, not the original results. Pocket2Mol's own paper reports -7.288 while TargetDiff's
re-evaluation of the same model on the same split reports -5.14, a gap of roughly 2.1 kcal/mol,
because Pocket2Mol scored with QVina under a different protocol. Rows from the two sources must
never be mixed.

**Two incompatible synthetic-accessibility conventions.** All eight papers above report the
NORMALISED score (10 - SA_Ertl) / 9, where higher is better and the reference set sits at 0.73.
GenBench3D and PoseBusters report the RAW Ertl-Schuffenhauer score on 1 to 10, where lower is
better. A value of 0.73 and a value of 3.4 describe similar molecules. Converting is mandatory
before any comparison.

---

## The critiques, which are the reason this comparison matters

10. Buttenschoen M, Morris GM, Deane CM. PoseBusters: AI-based docking methods fail to generate
    physically valid poses or generalise to novel sequences. *Chem Sci* 2024;15(9):3130-3139.
    doi:10.1039/D3SC04185A
    *Concerns docking methods rather than generators. Finding: deep-learning methods produce
    physically implausible poses despite sub-2 A RMSD, and no deep-learning method yet
    outperforms classical docking once physical plausibility is required.*

11. Baillif B, Cole J, McCabe P, Bender A. Applying atomistic neural networks to bias conformer
    ensembles towards bioactive-like conformations. *(GenBench3D.)* arXiv:2407.04424.
    **This is the critique that covers the generators above, and it is severe.** Validity3D, the
    fraction of generated molecules with a valid 3D conformation: DiffSBDD 0.4%, liGAN 2%,
    Pocket2Mol 9%, 3D-SBDD 10%, TargetDiff 11%, ResGen 11%. The abstract states the range
    directly: "only between 0% and 11% of generated molecules have valid conformations", and
    local relaxation in the pocket raises every model by at least 40%.

    Median MMFF94s strain energy, **verified against the source 2026-09-21**: CrossDocked test
    set **102.5**, Pocket2Mol **194.9** (best of the generators), 3D-SBDD **592.2**, TargetDiff
    **1241.7**, DiffSBDD **1243.1** kcal/mol.

    CORRECTION. An earlier revision of this file recorded "133 kcal/mol (Pocket2Mol) to 2006
    (DiffSBDD), against 56 for the training set". **None of those three numbers appears in
    GenBench3D.** They were carried into a figure caption and a console summary, where "2006
    against our 12" was a comparison resting on a value the cited paper does not contain. The
    verified numbers above replace them. The provenance of the withdrawn figures is unknown and
    they must not be reinstated from memory.

    NOTE ON COMPARING OUR OWN STRAIN TO THESE. We do not. GenBench3D computes MMFF94s medians
    over its own generated sets under its own relaxation protocol; a single axis carrying both
    would be a cross-method ranking this work declines to make. Our value is reported alone.

    On raw in-place poses ResGen's median Vina collapses to -0.47, and its Vina correlates with
    clash count at Spearman 0.85.

12. Gao B, Ren M, Ni Y, et al. Rethinking Specificity in SBDD: Leveraging Delta Score and Energy
    -based Models. *ICML* 2024. arXiv:2403.12987.
    *States that "specificity has been a largely overlooked factor in benchmarking SBDD
    previously" and that generated molecules "bind to almost every protein pocket with high
    affinity". Of 20 top-scoring generated molecules, only 2 of 20 scored best against their
    intended target under cross-docking. Supplies the control that matters: **random molecules
    drawn from the training set beat the reference ligand 30.83% of the time**, so a reported
    high-affinity rate near 38% is barely distinguishable from random sampling.*

13. Gao B, et al. From Theory to Therapy: Reviewing the Current State of AI Methods in Drug
    Design. arXiv:2406.08980.
    *"Vina scores can be significantly inflated simply by increasing the number of atoms in a
    molecule, revealing a susceptibility to overfitting." MolCRAFT concedes the same point in its
    own paper, noting that Vina Dock "can potentially be hacked by generating larger molecules",
    and demonstrates it: DecompDiff attains the second-best Vina Dock at 29.4 heavy atoms against
    the reference's 22.8. Any comparison that is not size-controlled is confounded.*

---

## Selectivity: what is and is not prior art

Verified: **AR, Pocket2Mol, TargetDiff, DecompDiff, liGAN and MolCRAFT perform no off-target or
paralogue evaluation at all.**

One exception must be credited rather than overlooked:

14. Schneuing et al., DiffSBDD (ref 8) performs explicit negative design, optimising against the
    on-target kinase BIKE (PDB 4W9W) while simultaneously designing against the structurally
    similar off-target MPSK1 (PDB 2BUJ). Across five iterations the on-target score improves from
    -7.2 to -13.9 while the off-target degrades from -10.8 to -8.7. This is a single kinase-pair
    case study using docking scores as the selectivity readout, not a benchmark-wide evaluation,
    but it is a genuine paralogue-selectivity experiment and predates ours.

ResGen's full text was inaccessible, so whether it evaluates selectivity is UNVERIFIED. A later
publication from the same group designs selective LTK ligands; that is a different paper and must
not be attributed to ResGen.

---

## How to position a comparison honestly

Four things follow from the above and should govern how our own numbers are reported.

- **Report size-controlled scores.** MolCRAFT's stratification is the model to follow, because
  Vina inflates with heavy-atom count (refs 12, 13).
- **Report in-place Vina Score alongside Vina Dock.** They differ by 2 kcal/mol or more, and
  quoting only the redocked value flatters the method (refs 4, 5, 11).
- **Report physical validity and strain**, not only affinity, since 89% to 99.6% of molecules
  from the published models fail a conformer-validity check (ref 11).
- **Benchmark the high-affinity rate against the 30.83% random-molecule floor**, not against
  zero (ref 12).
