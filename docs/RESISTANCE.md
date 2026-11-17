# Resistance mechanisms of the handbook Sanger assays

Hand-written on 2026-11-17 (not generated, unlike [`PRIMERS.md`](PRIMERS.md)). Source: the "Prinzip", "Referenz" and
"Auswertung" sections of the protocols in `data/Arbeitshandbücher` that are listed in `PRIMERS.md`, read through
`cosmic_seq.handbook.handbook_text`. Only the assays that look for resistance are covered; the identification assays
(16S 5' part and 3' end, rpoB-NTM, hsp65, 16S-23S ITS) are not.

Two kinds of mutations are kept apart:

- **Protocol**: written in the protocol text (mutation, codon or nucleotide position). Tables or figures that exist in
  the protocol only as pictures (eis, folC, NZM 21 probes) cannot be read from the text and are marked as such.
- **Literature**: what I add from the WHO catalogue of mutations in *M. tuberculosis* complex, 2nd edition 2023 (which
  32.43616 cites) and from the references the protocols name. Not in the handbook; check before use.

Numbering: *M. tuberculosis* H37Rv (gene coordinates: nucleotide 1 = first base of the start codon, negative = upstream)
unless stated. rRNA positions are given as in the protocols, often in *E. coli* numbering. The protocol positions
were checked against the codon they name (e.g. embB nt 916 = codon 306, thyA nt 604 = codon 202).

## Overview

| Drug | Gene / region | Mechanism | Mutations named in the protocol | Main mutations (literature) | Scope | Amplicon(s) | Protocol |
|---|---|---|---|---|---|---|---|
| Rifampicin | rpoB 81-bp core (RRDR) | target alteration (RNA polymerase β) | core region codons 507–533 (*E. coli* numbering), H37Rv rpoB nt 1276–1356 | S450L, H445Y/D/R, D435V/Y, L430P, L452P (= S531L, H526Y, D516V, L511P, L533P *E. coli*) | Mtb complex | rpoB_TR1_TR2b | 32.43612 |
| Rifampicin | rpoB 5' region | target alteration | V146F (*E. coli*) = "V176" (old Mtb numbering) | V170F (current H37Rv numbering), I491F lies outside both amplicons | Mtb complex | rpoB_Tb176f_Tb176r | 32.43619 |
| Rifampicin | rpoB, Xpert MTB/RIF Ultra region | target alteration | probes rpo1–4 over the RRDR (figure only) | as for the 81-bp core | Mtb complex, from the Xpert cartridge | rpoB_Xpert_F1_R, _F2_R2, _F2_R | NZM 21 |
| Isoniazid | katG codons 291–523 | loss of prodrug activation (catalase-peroxidase) | aa 315; aa 279 lies outside the amplicon | S315T (also S315N/I/R/G) | Mtb complex | katG_core | 32.43610 |
| Isoniazid, ethionamide | fabG1–inhA promoter | target overexpression (InhA) | promoter positions −8 and −15 ("most frequent") | fabG1 c.−15C>T, c.−8T>C/A, c.−17G>T; fabG1 L203L (c.609G>A) | Mtb complex | inhAprom_for_rev | 32.43609 |
| Isoniazid, ethionamide | inhA gene | target alteration (enoyl-ACP reductase) | "10 different mutations", via TBDReamDB (not listed) | S94A, I21T/V, I194T | Mtb complex | inhA_for_rev | 32.43609 |
| Ethambutol | embB | target alteration (arabinosyltransferase) | ~70 % of EMB resistance in embB, codons 306, 406, 497: M306L/V/I (nt 916, 918), G406S/D (nt 1216, 1217), Q497R (nt 1490) | M306V/I/L, G406A/C/D/S, Q497R/K/P, D354A, Y319C | Mtb complex | embB_1f_1r, embB_2f_2r | 32.43603 |
| Ethambutol | embB codon 306 (alternative) | target alteration | M306 ATG → ATA/ATC/ATT/GTG/CTG (probe set) | as above | Mtb complex (marked "not finished") | embB_151_131 | NZM 28 |
| Pyrazinamide | pncA + promoter | loss of prodrug activation (pyrazinamidase) | amino-acid changes in the gene or promoter; evaluation with the WHO catalogue 2023 | no hotspot: missense, nonsense and frameshifts spread over the 561-bp gene and promoter (e.g. H57D, D49N, Q10P, c.−11A>G) | Mtb complex | pncA_Ber | 32.43616 |
| Pyrazinamide | panD | target alteration (aspartate decarboxylase) | M117I (G351A), L136R (T407G), E126* (G376T), del 393G, V138E/A (T413A/C); in-vitro mutants of H37Ra (Shi 2014) | not graded as resistance markers in the WHO catalogue | Mtb complex | panD_F_R | NZM 5 |
| Pyrazinamide | rpsA | target (ribosomal protein S1), disputed | amino-acid changes "may be associated" | role disputed (Alexander 2012, cited in the protocol); not a WHO marker | Mtb complex | rpsA_1 … rpsA_4, rpsA_M5 | NZM 29 |
| Fluoroquinolones | gyrA QRDR | target alteration (DNA gyrase A) | nt 220, 238, 263, 265, 269, 271, 280, 281, 284, 305, 325; table: A74S, A90V, S91P, D94N/H/Y/A/G, P102H | A90V, S91P, D94G/A/N/Y/H, G88C/A; **S95T and E21Q are lineage polymorphisms, not resistance** | Mtb complex | gyrA_174f1_72r2 | 32.43611 |
| Fluoroquinolones | gyrA QRDR | target alteration | none ("may explain variable resistance"); refers to the Mtb positions above | intrinsic resistance is mostly not gyrA-based (like Mtb, *M. abscessus* has Ala at the Mtb A90 position, *E. coli* S83); acquired changes at the Mtb A90/D94 equivalents raise the MIC further | *M. abscessus* | gyrA_Mabs | NZM 16 |
| Fluoroquinolones | gyrA | target alteration | none (primers only) | QRDR changes at the Mtb A90/S91/D94 equivalents | *M. chimaera* | gyrA_Mchim | NZM 15 |
| Fluoroquinolones | gyrA QRDR | target alteration | none (primers only) | as above | genus *Mycobacterium* (degenerate primers) | gyrA_Pri8_Pri9, gyrA_Pri8_Pri10 | NZM 3 |
| Amikacin, kanamycin, capreomycin | rrs (16S rRNA) 3' part | target alteration (ribosomal A site) | Mtb rrs 1401, 1402, 1403, 1483, 1484 (= *E. coli* 1408, 1409, 1410, 1490, 1491); table: A1401G (KAN/AMK/CAP), C1402T/A (KAN, CAP), G1484T (KAN) | rrs a1401g (high-level AMK/KAN), c1402t, g1484t | genus *Mycobacterium* (primers), evaluation for Mtb | rrs_289_261 | 32.43613 |
| Kanamycin (low level) | eis promoter and gene | enzymatic modification: overexpression of the aminoglycoside acetyltransferase Eis | promoter mutations of Zaunbrecher 2009 and one more of Engström 2011 (in a figure, not readable from the text) | eis c.−10G>A, c.−12C>T, c.−14C>T, c.−37G>T (KAN; not AMK) | Mtb complex | eisprom_F_R, eis_F_R | NZM 2 |
| Capreomycin | tlyA | loss of function (2'-O-methyltransferase of 16S/23S rRNA) | compare with the wild type; no mutation named | no hotspot: frameshifts, nonsense and missense spread over the gene | Mtb complex | tlyA_F3_R3, tlyA_1f_1r, tlyA_1fmod_1r, tlyA_2f_2r, tlyA_2fmod_2rmod | 32.43621 |
| Streptomycin (high level) | rpsL | target alteration (ribosomal protein S12) | aa 43 and 88 ("high confidence") | K43R, K88R/Q/T | Mtb complex | rpsL_S13_S16 | 32.43615 |
| Streptomycin (low level) | gidB (gid) | loss of function (16S rRNA methyltransferase) | "several mutations"; via TBDReamDB | spread over the gene; several changes are lineage polymorphisms (e.g. E92D, A205A) | Mtb complex | gidB_5_3 | NZM 14 |
| Ethionamide | ethA | loss of prodrug activation (monooxygenase) | compare with the wild type; "many mutations" (TBDReamDB) | no hotspot: frameshifts, nonsense and missense over the 1470-bp gene | Mtb complex | ethA_f1_r1, ethA_f2_r2 | 32.43607 |
| Ethionamide | ethR | repression of ethA (regulator) | compare with the wild type | rare; minor contribution | Mtb complex | ethR_F02_R07 | 32.43608 |
| Para-aminosalicylic acid | thyA | loss of function (thymidylate synthase) | "high confidence": nt 604 ACC/GCC T202A, nt 664 CGC/GGC R222G | loss-of-function changes over the gene; **T202A is a lineage marker (LAM) and not resistance** (Mathys 2009, cited in the protocol) | Mtb complex | thyA_F1_R1 | 32.43614 |
| Para-aminosalicylic acid | folC | target alteration (dihydrofolate synthase) | table/figure of Zhang 2015 (not readable from the text) | I43A/T, R49W, S150G, F152L/S, E153A/G/K | Mtb complex | folC_F1_R1, folC_F2_R2 | NZM 30 |
| Linezolid | rplC | target alteration (ribosomal protein L3) | T460C (Beckert 2012) | C154R (= T460C) | Mtb complex | rplC_5_3 | 32.43617 |
| Linezolid | rrl (23S rRNA), peptidyl transferase centre | target alteration | G2032 (*E. coli*) = Mtb 2270; double mutations G2032A+C2499A, G2032A+U2504G, C2055A+U2504G, C2055A+A2572U (Long 2010) | rrl g2270t/c, g2814t (Mtb numbering) | genus *Mycobacterium* (primers), evaluation for Mtb | 23S_SP1F_SP2Rmod | 32.43604 |
| Bedaquiline | atpE | target alteration (ATP synthase subunit c) | drug-binding residues D28, L59, E61, A63, I66 (Segala 2012) | A63P, I66M, D28V/G/A, E61D | Mtb complex | atpE_F_R | 32.43618 |
| Bedaquiline, clofazimine | Rv0678 (mmpR5) | derepression of the MmpS5–MmpL5 efflux pump | S63R (nt 189 C>A), R134* (nt 400 C>T) (Hartkoorn 2014) | no hotspot: frameshifts (e.g. c.138_139insG, c.141_142insC), nonsense and missense over the gene and promoter | Mtb complex | Rv0678_F4_R1, Rv0678_F1_R2, Rv0678_F2_R2 | 32.43620 |
| Delamanid | ddn (Rv3547) | loss of prodrug activation (deazaflavin-dependent nitroreductase) | none (primers and run only) | loss-of-function changes spread over ddn; also in fgd1, fbiA/B/C/D (not covered) | Mtb complex | ddn_Rv3547_F_R | NZM 1 |
| Clarithromycin (macrolides) | rrl (23S rRNA) domain V | target alteration (peptidyl transferase loop) | A2058G/C (*E. coli* numbering), positions 2058/2059 marked | A2058G/C/T, A2059G/C/T (*E. coli*) | NTM (*M. intracellulare*, other NTM after checking the primers) | 23S_MCLR19F_MCLR21R | 32.43605 |
| Macrolides (inducible) | erm(41) | target modification: 23S rRNA methylation at A2058 (*E. coli*) | T/C polymorphism at position 28; also subspecies differentiation (*abscessus*, *bolletii*, *massiliense*) | T28 (Trp10) = functional erm(41), inducible resistance; C28 (Arg10) = non-functional; *massiliense*: truncated erm(41), no inducible resistance | *M. abscessus* complex | erm41_f2_r2, erm41_f_r | 32.43622 |
| Macrolides (inducible) | erm(41) promoter region | as above | position 28; hairpin in the promoter needs a denaturing step before sequencing | as above | rapidly growing mycobacteria | erm41prom_for_rev | NZM erm41 promoter |
| Macrolides (inducible) | erm(38) / erm(39) / erm(40) | target modification: 23S rRNA methylation at A2058 (*E. coli*) | presence of the gene (band on the gel) = inducible resistance, reported even when the phenotype is susceptible; sequencing optional | presence/absence, not a point mutation | *M. fortuitum* group | erm39_CME1Y_CME2 | NZM 18 |
