# ALKYL

```
 ░▒▓██████▓▒░░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░
░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░
░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░
░▒▓████████▓▒░▒▓█▓▒░      ░▒▓███████▓▒░ ░▒▓██████▓▒░░▒▓█▓▒░
░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░  ░▒▓█▓▒░   ░▒▓█▓▒░
░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░  ░▒▓█▓▒░   ░▒▓█▓▒░
░▒▓█▓▒░░▒▓█▓▒░▒▓████████▓▒░▒▓█▓▒░░▒▓█▓▒░  ░▒▓█▓▒░   ░▒▓████████▓▒░
```

**Computational chemistry plugin for [Claude Code](https://claude.ai/code).**

Transforms Claude Code into a specialized computational chemistry assistant — molecular modeling, drug discovery, quantum chemistry, MD simulations, and machine learning for molecules. No separate command, no wrapper, just `claude`.

## Install

```bash
git clone https://github.com/YOUR_USERNAME/alkyl
cd alkyl
bash install.sh
```

Open a new `claude` session. ALKYL is active.

## Uninstall

```bash
bash uninstall.sh
```

---

## What ALKYL provides

### Identity & behavior
- Specialized identity as **ALKYL**, a computational chemistry assistant
- Defaults to IUPAC nomenclature, SMILES notation, and computational cost awareness
- Adapts language to user context (English/French)

### 22 chemistry skills

Skills are loaded on demand in Claude Code sessions. Each skill covers a specific domain with practical code patterns and theoretical context.

| Skill | Description |
|-------|-------------|
| `rdkit` | Complete RDKit usage: molecule I/O, descriptors (MW, cLogP, TPSA, QED), Morgan and MACCS fingerprints, 2D/3D conformer generation, substructure search, reactions, SVG/PNG visualization |
| `ase` | Atomic Simulation Environment: Atoms objects, geometry optimization (BFGS/LBFGS/FIRE), NVE/NVT/NPT molecular dynamics (Langevin, Berendsen), NEB/AutoNEB transition states, vibrational analysis, thermochemistry, interface to ORCA/xTB/GPAW/LAMMPS |
| `mdanalysis` | MD trajectory analysis: Universe/AtomGroup selection language, RMSD/RMSF/alignment, contact analysis (Q-value), hydrogen bond analysis, Ramachandran/DSSP, PCA free energy landscapes, RDF, MSD diffusion, protein-ligand interaction workflow |
| `openbabel` | Open Babel format conversion (146 formats), 3D structure generation (MMFF94/UFF/GAFF), conformer search, protonation state at pH, SMARTS-based filtering, FP2/FP3/FP4/MACCS fingerprints, molecular descriptors, RDKit interoperability |
| `deepchem` | DeepChem molecular ML: MoleculeNet datasets, featurization (ECFP/GraphConv/Weave), graph neural networks, QSAR/QSPR model training and evaluation |
| `docking` | Protein-ligand docking and virtual screening: receptor preparation with pdbfixer/propka3, AutoDock Vina Python API and CLI, Gnina CNN rescoring, meeko PDBQT preparation, batch parallel docking, ProLIF interaction fingerprints, RMSD pose clustering, ensemble docking on MD snapshots |
| `force-fields` | Molecular mechanics force fields: AMBER/CHARMM/OPLS-AA/SMIRNOFF families, OpenMM simulation setup (LangevinMiddleIntegrator, NPT barostat, DCD/XTC reporters), OpenFF Sage 2.2 with SMIRNOFF, GAFF2 parameterization (antechamber/acpype), AM1-BCC/RESP charge methods, water models (TIP3P/OPC), HMR |
| `qm-dft` | Quantum chemistry: DFT functional/basis set selection (Jacob's ladder, D3BJ dispersion), ORCA 6.0 input/output/parsing (Opt, Freq, TS, TD-DFT, NMR, DLPNO, solvation), xTB/GFN2 (CLI, tblite Python API, CREST conformers/tautomers, pKa), PySCF (HF/DFT/MP2/CCSD, GIAO NMR, ESP/CHELPG), standard workflows (opt→freq→SP, barriers, UV-Vis) |
| `homology-modeling` | Protein structure modeling: template search with HHblits/BLAST, pairwise alignment (BLOSUM62/PIR format), MODELLER 10 automodel/loopmodel/multi-template/DOPE ranking, AlphaFold2 via ColabFold CLI, ESMFold Python API, structure quality (pLDDT, Ramachandran, MolProbity), structure preparation (pdbfixer, propka3, HIS tautomers, disulfides, terminal capping) |
| `free-energy` | Alchemical free energy calculations: thermodynamic cycles (RBFE/ABFE), FEP/TI/BAR/MBAR estimators, OpenMMTools AlchemicalFactory/MultiStateSampler/HREX, RBFE network design with LOMAP, perses and openfe campaign management, ABFE double-decoupling with Boresch restraints, pymbar (overlap matrix diagnostics, convergence, autocorrelation) |
| `pharmacophore` | Pharmacophore modeling: feature types (HBD/HBA/AR/HYD/POS/NEG), FDEF format, RDKit Pharm2D Gobbi fingerprints, Pharm3D 3D matching with EmbedLib, structure-based pharmacophore from ProLIF/PLIP interactions, ligand-based alignment with O3A and DBSCAN clustering, virtual screening pipeline (conformers → scoring → exclusion volumes → enrichment EF/ROC) |
| `generative-design` | De novo molecular generation: SELFIES always-valid grammar, SMILES language models (LSTM/GPT2/ChemGPT), REINVENT 4 RL with scoring components (QED, SA, docking oracle, custom Python), JT-VAE latent space Bayesian optimization, structure-based generation (DiffSBDD, TargetDiff, Pocket2Mol, DiffLinker), MOSES/GuacaMol evaluation metrics |
| `mmpa` | Matched Molecular Pair Analysis: MMP definition and fragmentation (Hussain-Rea algorithm), SMIRKS transform notation, mmpdb 4 CLI workflow (fragment→index→loadprops→transform→analyze), RDKit programmatic MMP generation, activity cliff detection (SALI score), bioisostere table, focused library generation, REINVENT/docking integration |
| `uncertainty-qsar` | Reliable QSAR predictions: epistemic vs. aleatoric decomposition, conformal prediction with MAPIE (split/CV+, coverage guarantee, Mondrian stratification), Gaussian processes with GPyTorch TanimotoKernel, MC Dropout (T=50 inference passes), deep ensembles (M=5 seeds), heteroscedastic head (NLL loss), Laplace approximation, applicability domain (kNN Tanimoto, Williams plot leverage, Mahalanobis), OECD Principle 3 compliance |
| `active-learning` | Active molecular learning: query strategies (UCB/EI/BALD/QBC/Core-Set), batch mode (DPP/greedy submodular/cluster-then-rank), molecular representations for acquisition, RF/GP/conformal uncertainty, docking oracle integration (Vina/Gnina, ~50× screening speedup), BEDROC/EF evaluation, DMTA cycle management (Design-Make-Test-Analyze, batch composition, stopping criteria, round reports) |
| `py3Dmol` | Interactive 3D molecular visualization (3Dmol.js): loading PDB/SDF/SMILES, cartoon/stick/sphere/surface styles (SES/SAS/VDW), selection language (chain/resi/resn/within), color schemes (spectrum/b-factor/rainbow/carbonColor), docking pose batch viewer with ipywidgets slider, pharmacophore overlay with spheres and labels, conformer animation, PNG/HTML export, NGLview for MD trajectories |
| `coarse-grained` | Coarse-grained MD simulations (MARTINI 3): CG resolution levels and mapping schemes, protein CG with martinize2 (ElNeDyn elastic network, Go-MARTINI), lipid membrane assembly with insane.py (POPC/POPE/POPS/CHOL/PIP2, asymmetric bilayers), protein-membrane embedding, GROMACS CG workflows (semiisotropic pressure coupling), backward.py backmapping CG→AA with OpenMM refinement, membrane analysis (thickness, APL, Scd order parameter, lateral diffusion, density profiles) |
| `binding-kinetics` | Drug-target binding kinetics: kon/koff/KD/residence time theory, kinetic selectivity (Copeland framework), two-state/induced-fit/conformational-selection models, SPR data fitting (Langmuir 1:1, two-state, Biacore CSV parsing), ITC data analysis (Wiseman isotherm, ΔG/ΔH/ΔS/ΔCp fitting, van't Hoff, enthalpy-entropy compensation), residence time by MD (τRAMD with HTMD, funnel metadynamics PLUMED, unbinding pathway analysis), kinetic QSAR (RF/GP koff models, koff cliffs, MMPA kinetic SAR) |
| `fbdd` | Fragment-Based Drug Design: Hann complexity model, Rule of 3 filters (MW/cLogP/HBD/HBA), Ligand Efficiency metrics (LE/LLE/LLEAT/BEI/SEI/GE/LELP), fragment library design (Ro3+PAINS+reactive filters, Fsp3/3D diversity, ESOL solubility, MaxMin selection), fragment docking (high exhaustiveness, Gnina rescoring, RMSD clustering, hotspot validation), growing/linking/merging strategies (R-group enumeration, SMARTS reactions, MCS merging, REINVENT scaffold constraint), efficiency plot (Abad-Zapatero) and round reports |
| `chem-brainstorm` | Flexible computational chemistry brainstorming guide: classify problem type → audit available data → map tools → generate directions → sanity checks → literature. Includes 4 rigid protocols (molecule evaluation, SAR hypothesis, reaction design, pipeline architecture). Integrates ALKYL scripts and all MCP tools (ChEMBL, OpenTargets, bioRxiv, ClinicalTrials) |
| `daylight-theory` | Molecular informatics theory (based on Daylight Theory Manual): complete SMILES specification (atoms, bonds, branches, rings, isomeric, chirality), SMARTS query language (all primitives D/H/R/r/v/X/x/+/-/#n/a/A, operators !&;,, recursive SMARTS, reaction queries), SMIRKS transform grammar, path fingerprints and similarity metrics (Tanimoto/Dice/Tversky/Cosine + 15 variants), molecular graph model and aromaticity (Hückel 4N+2) |
| `lit-rescue` | Literature search of last resort when hallucination risk is >20%: Perplexity→bioRxiv→PubMed waterfall, 7 query types (METHOD/PARAM/BUG/THEORY/PROTOCOL/BENCHMARK/DOMAIN), structured prompts with anti-hallucination framing, confidence reporting (★★★ to ☆☆☆), mandatory negative result block when no source found |

### 23 standalone Python scripts

All scripts in `scripts/`. Each requires only RDKit (and stdlib). Run with any Python ≥ 3.9 environment that has RDKit installed.

| Script | Description |
|--------|-------------|
| `chem_convert.py` | Convert molecules between SMILES, SDF, InChI, InChIKey, and SVG. Supports single molecules and batch files. |
| `chem_props.py` | Compute molecular properties: MW, cLogP, TPSA, HBD, HBA, RotBonds, QED. Checks Lipinski Ro5 and PAINS alerts. Computes Morgan (ECFP4) and MACCS fingerprints. |
| `chem_fetch.py` | Fetch molecular data from PubChem (by name, CID, or InChI) and ChEMBL (by name or ChEMBL ID). Returns SMILES, synonyms, and properties. Uses stdlib urllib only. |
| `chem_3d.py` | Generate 3D conformers using ETKDGv3, then minimize with MMFF94 or UFF. Outputs SDF. Configurable number of conformers and energy window. |
| `chem_qm.py` | Generate ORCA or Gaussian input files from SMILES (with automatic 3D embedding). Parse ORCA output files for energy, frequencies, thermochemical data, and IR spectrum (--parse-ir flag). |
| `chem_batch.py` | Batch-process SDF/SMI/CSV files: compute descriptors, check Lipinski Ro5, flag PAINS. Supports --skip-invalid for robust handling of malformed inputs. Outputs CSV. |
| `chem_search.py` | Search a molecular library by substructure (SMILES or SMARTS), Tanimoto similarity (configurable threshold), or exact match. Inputs: SDF or SMILES file. Outputs matches as SMILES or SDF. |
| `chem_standardize.py` | Standardize molecules: desalt (largest fragment), neutralize charges, canonicalize SMILES. Uses RDKit MolStandardize LargestFragmentChooser and Uncharger. |
| `chem_analyze.py` | Comprehensive single-molecule analysis: molecular formula, 16 functional group SMARTS, ring systems, stereocenters, QED, Synthetic Accessibility (SA) score, Bertz complexity. |
| `chem_scaffold.py` | Extract Murcko scaffold, generic scaffold (all non-hydrogen atoms replaced by carbon), and BRICS fragments from a molecule. |
| `chem_compare.py` | Compare two molecules: Maximum Common Substructure (rdFMCS), Tanimoto fingerprint similarity, and delta properties (B minus A for MW, cLogP, TPSA, HBD, HBA). |
| `chem_filter.py` | Filter a molecular library by drug-likeness rules: Lipinski Ro5, Veber (oral bioavailability), Egan (absorption), Ghose (requires ≥20 heavy atoms), and PAINS alerts. |
| `chem_react.py` | Apply SMARTS reaction transforms to reactants using RunReactants. Deduplicates and sanitizes products. Supports multi-step and multi-reactant reactions. |
| `chem_tautomers.py` | Enumerate tautomers of a molecule using RDKit TautomerEnumerator. Returns the canonical tautomer and the full enumerated list with counts. |
| `chem_enum.py` | Enumerate stereoisomers of a molecule (unique=True). Caps output at configurable max_isomers to prevent combinatorial explosion. |
| `chem_pka.py` | Estimate pKa values using SMARTS-based ionization rules (acids and bases). Computes Henderson-Hasselbalch pH-speciation and dominant protonation state at target pH via RDKit MolStandardize. |
| `chem_metabolism.py` | Predict CYP450 metabolic soft spots using 12 SMARTS rules across five isoforms (CYP3A4, CYP2D6, CYP2C9, CYP1A2, UGT/SULT). Highlights vulnerable positions on the molecule. |
| `chem_diversity.py` | Select a maximally diverse subset of molecules using the MaxMin algorithm (O(k·n) complexity). Supports Morgan (ECFP4) and MACCS fingerprints. Handles cases where k ≥ library size. |
| `chem_cluster.py` | Butina (Taylor-Butina) clustering of a molecular library by Tanimoto distance. Returns cluster IDs, centroids, and members. Configurable distance cutoff and fingerprint type. |
| `chem_rgroup.py` | R-group decomposition around a SMARTS core scaffold. Returns a table of R1/R2/... substituents per molecule, matched and unmatched counts. Uses RDKit RGroupDecompose. |
| `chem_admet.py` | Heuristic ADMET profiling: ESOL aqueous solubility (Delaney 2004), BBB penetration score, hERG cardiac risk (SMARTS alerts), P-gp substrate likelihood, plasma protein binding estimate. |
| `chem_highlight.py` | Draw a molecule as SVG or PNG with a SMARTS-matched substructure highlighted. Outputs SVG to stdout or saves SVG/PNG to file. |
| `chem_lead.py` | Compute ligand efficiency metrics (LE, LLE, BEI, LELP) from a CSV of structures and activities. Tracks metric evolution across optimization rounds with mean/delta per round. |

### Optional: Perplexity MCP

For grounded web search in chemistry literature workflows:

```bash
bash setup-perplexity.sh YOUR_API_KEY
```

This adds the `@perplexity-ai/mcp-server` to your Claude Code MCP settings. The `lit-rescue` skill uses it as its primary search source.

---

## Requirements

- [Claude Code](https://claude.ai/code)
- Python ≥ 3.9 with RDKit (for scripts and certain skills)
- Optional per workflow: ORCA, xTB, GROMACS, OpenMM, MODELLER

### Setup

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install rdkit pytest
```

## Tests

```bash
# Unit tests only (no network)
python -m pytest tests/ -m "not network" -v

# All tests including network
python -m pytest tests/ -v
```

---

## Project structure

```
ALKYL/
├── install.sh              # inject chemistry context into ~/.claude/CLAUDE.md
├── uninstall.sh            # remove the injected block
├── setup-perplexity.sh     # configure Perplexity MCP (optional)
├── config/
│   └── CLAUDE.md           # ALKYL identity, behavior, and full skill index
├── scripts/
│   ├── chem_convert.py     # format conversion (SMILES/SDF/InChI/SVG)
│   ├── chem_props.py       # molecular properties and fingerprints
│   ├── chem_fetch.py       # PubChem and ChEMBL data retrieval
│   ├── chem_3d.py          # ETKDGv3 conformer generation and minimization
│   ├── chem_qm.py          # ORCA/Gaussian input generation and output parsing
│   ├── chem_batch.py       # batch processing of molecular libraries
│   ├── chem_search.py      # substructure, similarity, and exact search
│   ├── chem_standardize.py # desalting, neutralization, canonicalization
│   ├── chem_analyze.py     # comprehensive single-molecule analysis
│   ├── chem_scaffold.py    # Murcko scaffold and BRICS fragmentation
│   ├── chem_compare.py     # MCS, Tanimoto, and property delta
│   ├── chem_filter.py      # Lipinski/Veber/Egan/Ghose/PAINS filters
│   ├── chem_react.py       # SMARTS reaction application
│   ├── chem_tautomers.py   # tautomer enumeration and canonicalization
│   ├── chem_enum.py        # stereoisomer enumeration
│   ├── chem_pka.py         # pKa estimation and protonation state
│   ├── chem_metabolism.py  # CYP450 metabolic soft spot prediction
│   ├── chem_diversity.py   # MaxMin diversity selection
│   ├── chem_cluster.py     # Butina clustering by Tanimoto distance
│   ├── chem_rgroup.py      # R-group decomposition around a SMARTS core
│   ├── chem_admet.py       # ADMET heuristics (ESOL, BBB, hERG, P-gp, PPB)
│   ├── chem_highlight.py   # SVG/PNG with SMARTS-matched atoms highlighted
│   └── chem_lead.py        # ligand efficiency metrics (LE/LLE/BEI/LELP) per round
└── skills/
    ├── rdkit/              # RDKit cheminformatics
    ├── ase/                # Atomic Simulation Environment
    ├── mdanalysis/         # MD trajectory analysis
    ├── openbabel/          # format conversion and filtering
    ├── deepchem/           # molecular machine learning
    ├── docking/            # virtual screening and docking
    ├── force-fields/       # AMBER/OpenMM/OpenFF/GAFF2
    ├── qm-dft/             # ORCA/xTB/PySCF quantum chemistry
    ├── homology-modeling/  # MODELLER/ColabFold/ESMFold
    ├── free-energy/        # FEP/MBAR/RBFE/ABFE
    ├── pharmacophore/      # pharmacophore modeling and VS
    ├── generative-design/  # de novo molecular generation
    ├── mmpa/               # matched molecular pair analysis
    ├── uncertainty-qsar/   # conformal prediction and GP uncertainty
    ├── active-learning/    # DMTA loop and active screening
    ├── py3Dmol/            # interactive 3D visualization
    ├── coarse-grained/     # MARTINI 3 and membrane simulations
    ├── binding-kinetics/   # SPR/ITC/τRAMD/kinetic QSAR
    ├── fbdd/               # fragment-based drug design
    ├── chem-brainstorm/    # workflow brainstorming guide
    ├── daylight-theory/    # SMILES/SMARTS/SMIRKS/fingerprints theory
    └── lit-rescue/         # literature search of last resort
```

---

## License

MIT — see [LICENSE](LICENSE).

## Acknowledgments

The skills in this repository draw on and are informed by the following works and their authors:

- **Daylight Theory Manual** — Daylight Chemical Information Systems (SMILES, SMARTS, SMIRKS, fingerprints)
- **RDKit documentation** — Greg Landrum and RDKit contributors
- **ASE documentation** — Ask Hjorth Larsen, Jens Jørgen Mortensen, and ASE contributors
- **MDAnalysis documentation** — Oliver Beckstein, Richard Gowers, and MDAnalysis contributors
- **MARTINI force field** — Siewert-Jan Marrink, Xavier Periole, D. Peter Tieleman, and CGMD community
- **OpenFF Sage / SMIRNOFF** — Open Force Field Initiative contributors
- **REINVENT** — AstraZeneca Molecular AI team
- **mmpdb** — Andrew Dalke and contributors
- **AlphaFold / ColabFold** — DeepMind, Sergey Ovchinnikov, Martin Steinegger
- **ORCA** — Frank Neese and the ORCA development team
- **Copeland binding kinetics framework** — Robert A. Copeland (*Evaluation of Enzyme Inhibitors in Drug Discovery*, Wiley)
- **Hussain-Rea fragmentation** — Jameed Hussain, Ceara Rea (J. Chem. Inf. Model., 2010)
- **Haussler Tanimoto kernel** — David Haussler (1999)
- **Conformal prediction** — Vladimir Vovk, Alexander Gammerman, Glenn Shafer (*Algorithmic Learning in a Random World*, Springer)
- All open-source tool authors and scientific communities whose work these skills build upon
