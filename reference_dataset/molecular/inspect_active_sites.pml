# M7 Independent Layer 5 Molecular Verification Script for PyMOL
# Usage: pymol inspect_active_sites.pml
#
# Epistemic & Structural Provenance Notes:
# 1. dsrB_2V4J_catalytic_cluster.pdb is a COMPOSITE ACTIVE-SITE EXCERPT of selected
#    chain fragments from three separate source crystal structures:
#      - Chain A: DsrB dissimilatory sulfite reductase (PDB 2V4J, UniProt P45575)
#      - Chain B: OahA oxaloacetate acetylhydrolase (PDB 3M0J, UniProt A2QSW8)
#      - Chain C: OxDC bicupin oxalate decarboxylase (PDB 1J58, UniProt O34714)
#    It is NOT a single crystallographic co-assembly of all three proteins.
# 2. For independent verification of the complete original crystallographic assemblies
#    from the RCSB PDB archive, execute in PyMOL:
#      fetch 2v4j 3m0j 1j58, async=0
# 3. csgA_curli_repeat_model.pdb is badged as [THEORETICAL_MODEL_FRAGMENT]: an 18-residue
#    (Ser43..Gln60) theoretical beta-solenoid repeat fragment of E. coli CsgA (UniProt P28307),
#    contrasted with full-length experimental cryo-EM curli structures PDB 8ENQ / 8ENR:
#      fetch 8enq 8enr, async=0
# 4. Visual inspection verifies structural coordinates and residue numbering,
#    not reaction kinetics or metabolic flux.

reinitialize
bg_color black

# 1. Load composite active-site excerpt (selected fragments from PDB 2V4J, 3M0J, 1J58)
load dsrB_2V4J_catalytic_cluster.pdb, m7_experimental_sites
# 2. Load 18-residue CsgA repeat [THEORETICAL_MODEL_FRAGMENT] (quarantined from experimental PDBs; contrast cryo-EM 8ENQ/8ENR)
load csgA_curli_repeat_model.pdb, csgA_theoretical_model

# Optional online full-assembly retrieval (uncomment when connected to RCSB PDB):
# fetch 2v4j 3m0j 1j58, async=0
# fetch 8enq 8enr, async=0

hide everything
show cartoon, m7_experimental_sites
show cartoon, csgA_theoretical_model
color marine, m7_experimental_sites and chain A
color teal, m7_experimental_sites and chain B
color slate, m7_experimental_sites and chain C
color gold, csgA_theoretical_model

# DsrB (PDB 2V4J excerpt, Chain A): C-X(2)-C-X(3)-C octapeptide CGNCYTMC (residues 261-268)
select dsrb_motif, m7_experimental_sites and chain A and resi 261-268
select dsrb_active_cys, m7_experimental_sites and chain A and resi 261+264+268
select dsrb_cofactors, m7_experimental_sites and chain A and resn FS4+SRM
show sticks, dsrb_active_cys or dsrb_cofactors

# OahA (PDB 3M0J excerpt, Chain B): M-G-F-R-I-M-I-F-S-F-A (residues 18-28, catalytic Arg21 + Ser26)
select oaha_motif, m7_experimental_sites and chain B and resi 18-28
select oaha_active_site, m7_experimental_sites and chain B and (resi 21+26 or resn MG+OXA)
show sticks, oaha_active_site

# OxdC (PDB 1J58 excerpt, Chain C): H-X-H-X(3)-E bicupin Mn-binding triad (His95, His97, Glu101)
select oxdc_active_site, m7_experimental_sites and chain C and (resi 95+97+101 or resn MN+OXL)
show sticks, oxdc_active_site

# CsgA (UniProt P28307, Chain A, [THEORETICAL_MODEL_FRAGMENT]; contrast cryo-EM 8ENQ/8ENR): 18-aa R1 repeat (43-60) & 12-aa core (49-60)
select csga_r1_repeat, csgA_theoretical_model and chain A and resi 43-60
select csga_12aa_core, csgA_theoretical_model and chain A and resi 49-60
select csga_conserved_gln_asn, csgA_theoretical_model and chain A and resi 49+52+54+60
show sticks, csga_conserved_gln_asn

# Measure key metal-ligand coordination distances
distance d_dsrb_cys261_fe, (m7_experimental_sites and chain A and resi 261 and name SG), (m7_experimental_sites and chain A and resn FS4 and name FE1)
distance d_oaha_ser26_mg, (m7_experimental_sites and chain B and resi 26 and name OG), (m7_experimental_sites and chain B and resn MG and name MG)
distance d_oxdc_his95_mn, (m7_experimental_sites and chain C and resi 95 and name NE2), (m7_experimental_sites and chain C and resn MN and name MN)

zoom dsrb_active_cys
