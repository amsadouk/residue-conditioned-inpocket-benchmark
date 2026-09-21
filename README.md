# Residue-Conditioned In-Pocket Design for Paralogue Selectivity

Supplementary code and data of the study .

Ahmed Amine SADOUK , Independent Researcher

This repository exist so that every number reported in the paper can be recomputed by anyone ,
without asking anything to the author and without any access to the design engine .

---

## What is in here and what is not

In here : the campaign records that every reported statistic is computed from , the analysis code
that compute them , the figure code , and the resolved structure of the 44 reference compounds .

Not in here : the generative design engine . It is proprietary and it is not distributed .

This separation is deliberate and it cost the reader nothing . The engine decide which molecule
to build next , the code of this repository only measure what was already built .
`verify_divergence.py` fetch a PDB entry and run two sequence alignment , `stats_structural.py`
read a table of SMILES and compute a maximum common substructure against a decoy set . None of
them contain any generative logic and none of them need it , because the molecules are already
recorded inside `data/campaigns/` .

If you want to verify a claim of the paper everything required for that is in this repository .

---

## Reproducing every number

Install the dependencies :

```
pip install -r requirements.txt
```

Then from the repository root :

| claim in the paper | command | expected result |
|---|---|---|
| scaffold retained in every final lead | `python analysis/stats_structural.py` | 17 of 17 , Wilson 82% to 100% |
| recognition chemistry preserved | same command | 19 of 44 , decoy rate 12.6% , p < 0.0001 |
| shared core against a decoy null | same command | 11 of 44 at p < 0.05 |
| covalent warhead at divergent position | same command | 4 of 6 directed rounds keep a warhead , 3 divergent |
| the drug to drug scale | `python analysis/stats_chemotype_scale.py` | 66 pairs , median 33% |
| leads above that scale | same command | 4 of 12 , ABL1 78% and EGFR 59% at p = 0.016 |
| ABL1 Cys388 diverge , EGFR Cys797 does not | `python analysis/verify_divergence.py` | ABL1 align to SRC Ile , EGFR conserved |
| which molecules each figure use | `python analysis/datasets.py` | the five populations , counted |
| the receptor preparation defect | `python figures_code/fig_prep_defect.py` | 1066 records in 6SBL , 9 of the 12 receptors |

`verify_divergence.py` and `fig_prep_defect.py` download the structures from the RCSB PDB so
these two need a network connection . Everything else run offline from the files of `data/` .

---

## Layout

```
analysis/          the statistics , the divergence check and the population manifest
figures_code/      the code of every figure of the paper
data/campaigns/    18 campaign records , one CSV per receptor pair
data/tables/       the 44 reference compounds with their InChIKey , and the comparison tables
results/           the outputs , so that yours can be diffed against ours
```

`results/` contain the statistics exactly as they were produced for the manuscript . Re running
the analysis should reproduce them . If it does not then that is a finding and the author would
like to be informed .

---

## The campaign records

`data/campaigns/panel_plan_<TARGET>_<ANTITARGET>.csv` is one row per round per campaign . The
columns that matter for checking the paper are the following .

| column | what it is |
|---|---|
| `round` | the design round , from 1 upward . Round 0 is the mined scaffold and it is in `scaffold_smiles` |
| `smiles` | the carried lead after that round |
| `scaffold_smiles` | the mined starting scaffold that every later round have to contain |
| `goal_residue` | the pocket residue that the round was directed at |
| `goal_feature` | the interaction type that it aimed to form |
| `goal_axis` | the engine own classification of that position as divergent or conserved |
| `heavy_atoms` | the molecular size |
| `affinity_target` , `affinity_counter` , `margin` | docking values , reported as context only |
| `sel_gate_reason` | the verdict of the per round counter screen |

Two remarks about this last column , both of them are stated in the paper and they are repeated
here so that nobody is misled by the data alone . The counter screen run **after** the lead have
already been carried and it record a verdict rather than blocking anything , so a round that is
marked REJECTED is still carried and the next round still build on it . The selectivity enter the
search earlier than that , through the divergent contacts axis of the Pareto front and through a
margin based eligibility preference .

---

## About the docking values

They are inside the records and they are used nowhere in the structural claims . The paper report
that the rank correlation of the docking score against the measured pChEMBL is positive on three
of the seven receptors of this panel , that it rise with the heavy atom count , and that it carry
no reactive term so it can not represent a covalent warhead at all . No statistic of `analysis/`
is computed on a docking score . They are supplied anyway because withholding them would be
worse .

---

## Corrections that the code carry

The analysis code recompute its own claims rather than trusting a stored label , because several
of them were found to be wrong during the study .

EGFR Cys797 is **conserved** and not divergent . The anti target of T790M is the wild type EGFR
and the two sequences differ at the position 790 and nowhere else . An earlier draft claimed the
opposite , `verify_divergence.py` now recompute every statement of that kind from the structure
and from the two alignments .

The residue numbering is not the same between the conventions . 6SBL omit the initiator
methionine so the deposited residue 206 is the canonical residue 205 , and 6NPV use the ABL1
isoform 1b numbering which run 19 residues ahead of the canonical one so the deposited residue
388 is the canonical residue 369 . Both numbers are reported .

Three strain values were withdrawn . An earlier version of the reference list recorded GenBench3D
medians of 56 , 133 and 2006 kcal/mol . None of those three number appear in that paper . The
verified values are 102.5 , 194.9 , 592.2 , 1241.7 and 1243.1 and they are in
`results/octara_references_sbdd.md` .

A receptor preparation defect corrupted 1066 ATOM records in 6SBL and affect 9 of the 12
receptors of the panel , inflating every reference docking score in the same direction . It is
reproducible with `figures_code/fig_prep_defect.py` against the structures that you fetch
yourself .

---

## Citation

If you use the campaign records or the analysis code please cite the paper . This repository is
archived with a DOI , use the DOI rather than the repository URL because a URL can move .

---

## Licence

The code and the data of this repository are released under the MIT Licence , see `LICENSE` . The
design engine is not a part of this release and it is not covered by it .
