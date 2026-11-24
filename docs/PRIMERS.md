# Primers of the mycobacterial Sanger assays

Written by `scripts/02_handbook_primers.py` on 2026-11-24 from `primers.tsv` and `amplicons.tsv`. Do not edit by hand.

Source: 38 protocols in `data/Arbeitshandbücher` (25 chapters of the molecular handbook, 13 research protocols kept outside it, protocol `NZM <n>`). Every sequence below was found in the text of its protocol by the script.

**115 primers, 59 primer pairs.** Only assays that end in Sanger sequencing are listed; probes and real-time or fragment-length assays are in the last section.

## Primers

| Primer | Target | Scope | Direction | Sequence 5'→3' | Role | Protocol |
|---|---|---|---|---|---|---|
| Mbak283f | 16S rRNA gene (rrs), 5' part | genus Mycobacterium | forward | `GAGTTTGATCCTGGCTCAGGA` | amplification | 32.43304 |
| Mbak264r | 16S rRNA gene (rrs), 5' part | genus Mycobacterium | reverse | `TGCACACAGGCCACAAGGGA` | amplification | 32.43304 |
| Mbak259r | 16S rRNA gene (rrs), 5' part | genus Mycobacterium | reverse | `TTTCACGAACAACGCGACAA` | amplification and sequencing | 32.43304 |
| Mbak14 | 16S rRNA gene (rrs), 5' part | genus Mycobacterium | forward | `GRGRTACTCGAGTGGCGAAC` | sequencing | 32.43304 |
| IMM F | 16S rRNA gene (rrs), 5' part | genus Mycobacterium, direct detection | forward | `TTAACACATGCAAGTCGAACGGAA` | amplification | 32.43325 |
| IMM R2 | 16S rRNA gene (rrs), 5' part | genus Mycobacterium, direct detection | reverse | `CAGTTAAGCCGTGAGATTTCACGAACA` | amplification | 32.43325 |
| Bact-1098s | 16S rRNA gene (rrs), 3' end | bacteria (broad range) | forward | `CGCAACGAGCGCAACCC` | amplification and sequencing | 32.43506 |
| Bact-1525a | 16S rRNA gene (rrs), 3' end | bacteria (broad range) | reverse | `AAGGAGGTGATCCARCC` | amplification and sequencing | 32.43506 |
| rpoB-2573f | rpoB | NTM | forward | `GGCAAGGTCACCCCGAAGGG` | amplification and sequencing | 32.43504 |
| rpoB-3337r | rpoB | NTM | reverse | `AGCGGCTGCTGGGTGATCATC` | amplification | 32.43504 |
| TB-11 | hsp65 (groEL2) | NTM | forward | `ACCAACGATGGTGTGTCCAT` | amplification and sequencing | 32.43505 |
| TB-12w | hsp65 (groEL2) | NTM | reverse | `CTTGTCGAASCGCATACCCT` | amplification and sequencing | 32.43505 |
| TB-12 | hsp65 (groEL2) | NTM | reverse | `CTTGTCGAACCGCATACCCT` | amplification and sequencing | 32.43505 |
| 16-23S-ITS-1f | 16S-23S ITS | NTM (M. fortuitum group) | forward | `ACCTCCTTTCTAAGGAGCACC` | amplification and sequencing | 32.43606 |
| 16-23S-ITS-2r | 16S-23S ITS | NTM (M. fortuitum group) | reverse | `GATGCTCGCAACCACTATCCA` | amplification and sequencing | 32.43606 |
| MCLR 19 F | 23S rRNA gene (rrl), clarithromycin | NTM | forward | `GTAGCGAAATTCCTTGTCGG` | amplification and sequencing | 32.43605 |
| MCLR 21R | 23S rRNA gene (rrl), clarithromycin | NTM | reverse | `TTCCCGCTTAGATGCTTTCAG` | amplification and sequencing | 32.43605 |
| erm41f2 | erm(41) | NTM (M. abscessus / M. chelonae) | forward | `GACCGGGGCCTTCTTCGTGAT` | amplification and sequencing | 32.43622 |
| erm41r2 | erm(41) | NTM (M. abscessus / M. chelonae) | reverse | `GACTTCCCCGCACCGATTCC` | amplification and sequencing | 32.43622 |
| erm41f | erm(41) | NTM (M. abscessus / M. chelonae) | forward | `TGGTATCCGCTCACTGATGA` | amplification and sequencing | 32.43622 |
| erm41r | erm(41) | NTM (M. abscessus / M. chelonae) | reverse | `GCGGTGGATGATGGAAAG` | amplification and sequencing | 32.43622 |
| embB-1f_mod2 | embB, ethambutol | Mtb complex | forward | `CAACTTCGTCGGGCTCAAG` | amplification and sequencing | 32.43603 |
| embB-1r | embB, ethambutol | Mtb complex | reverse | `CAGTAGGCGGGTTTGCTG` | amplification | 32.43603 |
| embB-2f | embB, ethambutol | Mtb complex | forward | `GCTAGTGTGCTGGCTGCTG` | amplification and sequencing | 32.43603 |
| embB-2r | embB, ethambutol | Mtb complex | reverse | `GGTGGGCAGGATGAGGTAG` | amplification | 32.43603 |
| SP1 F | 23S rRNA gene (rrl), clarithromycin / linezolid | genus Mycobacterium | forward | `CCTGCACGAATGGCGTAACG` | amplification and sequencing | 32.43604 |
| SP2 R mod | 23S rRNA gene (rrl), clarithromycin / linezolid | genus Mycobacterium | reverse | `CACCAGAGGTTCGTCCGTC` | amplification and sequencing | 32.43604 |
| ethAf1_mod | ethA, ethionamide | Mtb complex | forward | `TGGCAGCTTACTACGTGTCG` | amplification and sequencing | 32.43607 |
| ethAr1_mod | ethA, ethionamide | Mtb complex | reverse | `ACTTCTGGCAGGCGCTGTA` | amplification and sequencing | 32.43607 |
| ethAf2_mod | ethA, ethionamide | Mtb complex | forward | `GGAGACCATGGCCTACACC` | amplification and sequencing | 32.43607 |
| ethAr2_mod | ethA, ethionamide | Mtb complex | reverse | `GGCATCATCGTCGTCTGACT` | amplification and sequencing | 32.43607 |
| F02 | ethR, ethionamide | Mtb complex | forward | `CCACCGGACCTCCGTGAG` | amplification and sequencing | 32.43608 |
| R07 | ethR, ethionamide | Mtb complex | reverse | `GGCCGCGGTGCTCTTGAC` | amplification and sequencing | 32.43608 |
| inhAprom-for | inhA promoter (fabG1 upstream) and inhA gene, isoniazid | Mtb complex | forward | `ACAAACGTCACGAGCGTAACC` | amplification and sequencing | 32.43609 |
| inhAprom-rev | inhA promoter (fabG1 upstream) and inhA gene, isoniazid | Mtb complex | reverse | `GTTGGCGTTGATGACCTTCTC` | amplification and sequencing | 32.43609 |
| inhA-for | inhA promoter (fabG1 upstream) and inhA gene, isoniazid | Mtb complex | forward | `GAGGATGCGAGCTATATCTCC` | amplification and sequencing | 32.43609 |
| inhA-rev | inhA promoter (fabG1 upstream) and inhA gene, isoniazid | Mtb complex | reverse | `CGAAATGCAGGTAGTGCTCG` | amplification and sequencing | 32.43609 |
| katG_core_for | katG, isoniazid | Mtb complex | forward | `CTCCGCTGGAGCAGATGGGC` | amplification and sequencing | 32.43610 |
| katG_core_rev | katG, isoniazid | Mtb complex | reverse | `TCTTCCAGGGTGCGAATGACCT` | amplification and sequencing | 32.43610 |
| gyrA 174-f1 | gyrA, fluoroquinolones | Mtb complex | forward | `CGAGCGCGACGGTCTGCT` | amplification and sequencing | 32.43611 |
| gyrA 72-r2 | gyrA, fluoroquinolones | Mtb complex | reverse | `AATCCTCTTCTACCTCAACAACT` | amplification and sequencing | 32.43611 |
| TR1 | rpoB 81-bp core region, rifampicin | Mtb complex | forward | `TACGGTCGGCGAGCTGATCC` | amplification and sequencing | 32.43612 |
| TR2b | rpoB 81-bp core region, rifampicin | Mtb complex | reverse | `TACGGCGTTTCGATGAACC` | amplification and sequencing | 32.43612 |
| Tb176-f | rpoB 5' region, rifampicin | Mtb complex | forward | `CTTCTCCGGGTCGATGTCGTTG` | amplification and sequencing | 32.43619 |
| Tb176-r | rpoB 5' region, rifampicin | Mtb complex | reverse | `CGCGCTTGTCGACGTCAAACTC` | amplification and sequencing | 32.43619 |
| #289 | 16S rRNA gene (rrs) 3' end, amikacin / capreomycin | genus Mycobacterium | forward | `AAGTCGGAGTCGCTAGTAATCG` | amplification and sequencing | 32.43613 |
| #261 | 16S rRNA gene (rrs) 3' end, amikacin / capreomycin | genus Mycobacterium | reverse | `AAAGGAGGTGATCCAGCCGC` | amplification and sequencing | 32.43613 |
| thyAF1 | thyA, PAS | Mtb complex | forward | `CTGGGTGTATTCCCGTCGACT` | amplification and sequencing | 32.43614 |
| thyAR1 | thyA, PAS | Mtb complex | reverse | `ATCGTGTGCCCCATGGTGATCT` | amplification and sequencing | 32.43614 |
| thyAI1 | thyA, PAS | Mtb complex | forward | `GTTCAACATCGCCAGCTATGC` | sequencing | 32.43614 |
| S13 | rpsL, streptomycin | Mtb complex | forward | `GGCCGACAAACAGAACGT` | amplification and sequencing | 32.43615 |
| S16 | rpsL, streptomycin | Mtb complex | reverse | `GTTCACCAACTGGGTGAC` | amplification and sequencing | 32.43615 |
| pncA-Ber-fw | pncA, pyrazinamide | Mtb complex | forward | `CGCTGACCTGGGTGGAG` | amplification and sequencing | 32.43616 |
| pncA-Ber-rv | pncA, pyrazinamide | Mtb complex | reverse | `GCCGCCAACAGTTCATCC` | amplification and sequencing | 32.43616 |
| rplC_5` | rplC, linezolid | Mtb complex | forward | `GCTGCGGCTGGACGACTC` | amplification and sequencing | 32.43617 |
| rplC_3` | rplC, linezolid | Mtb complex | reverse | `CTCTTGCGCAGCCATCACTTC` | amplification and sequencing | 32.43617 |
| atpE Forward | atpE, bedaquiline | Mtb complex | forward | `GGACCCCACTATCGCTGCCGGCG` | amplification and sequencing | 32.43618 |
| atpE Reverse | atpE, bedaquiline | Mtb complex | reverse | `TTACTTGACGGGTGTAGCGAA` | amplification and sequencing | 32.43618 |
| Rv0678-F4 | Rv0678 (mmpR5), bedaquiline / clofazimine | Mtb complex | forward | `TGAAGTTCACGCCGGTCTGG` | amplification and sequencing | 32.43620 |
| Rv0678-R1 | Rv0678 (mmpR5), bedaquiline / clofazimine | Mtb complex | reverse | `AGTCGTCCTCTCCGGTTCGCT` | amplification and sequencing | 32.43620 |
| Rv0678-F1 | Rv0678 (mmpR5), bedaquiline / clofazimine | Mtb complex | forward | `AGCGTCAACGACGGGGTCG` | amplification and sequencing | 32.43620 |
| Rv0678-F2 | Rv0678 (mmpR5), bedaquiline / clofazimine | Mtb complex | forward | `GCGCCGAGCCCGACATC` | amplification and sequencing | 32.43620 |
| Rv0678-R2 | Rv0678 (mmpR5), bedaquiline / clofazimine | Mtb complex | reverse | `TCTCCGGTTCGCTGGCTGTATC` | amplification and sequencing | 32.43620 |
| tlyAF3 | tlyA, capreomycin | Mtb complex | forward | `AAGGCATCGCACGTCGTCTTTCC` | amplification and sequencing | 32.43621 |
| tlyAR3 | tlyA, capreomycin | Mtb complex | reverse | `TGTCGCCCAATACTTTTTCTACGC` | amplification and sequencing | 32.43621 |
| Mtb_tlyA-1f | tlyA, capreomycin | Mtb complex | forward | `TCTCTGGCCGAACTCGAA` | amplification and sequencing | 32.43621 |
| Mtb_tlyA-1r | tlyA, capreomycin | Mtb complex | reverse | `ACGGTAGCCAACGAGATGAA` | amplification and sequencing | 32.43621 |
| Mtb_tlyA-1f-mod | tlyA, capreomycin | Mtb complex | forward | `CAACGCTGGTCGTTGATG` | amplification and sequencing | 32.43621 |
| Mtb_tlyA-2r | tlyA, capreomycin | Mtb complex | reverse | `GCGAAGCGCAATTTTATTGT` | amplification and sequencing | 32.43621 |
| Mtb_tlyA-2f | tlyA, capreomycin | Mtb complex | forward | `GTACTGCTGGACCGTGGTG` | amplification and sequencing | 32.43621 |
| Mtb_tlyA-2r-mod | tlyA, capreomycin | Mtb complex | reverse | `CGTCAACCACCTCGATCTC` | amplification and sequencing | 32.43621 |
| Mtb_tlyA-2f-mod | tlyA, capreomycin | Mtb complex | forward | `GGAGCGCACAAACTAGTCG` | amplification and sequencing | 32.43621 |
| Mtb_tlyA-3f | tlyA, capreomycin | Mtb complex | forward | `GTTCCACTGGTGAAGCCGCAGTT` | amplification and sequencing | 32.43621 |
| Rv3547_F | ddn (Rv3547), delamanid | Mtb complex | forward | `CGAGCGCACCGACCAGAGC` | amplification and sequencing | NZM 1 |
| Rv3547_R | ddn (Rv3547), delamanid | Mtb complex | reverse | `GCATGGCCCGCAGGTGGACAA` | amplification and sequencing | NZM 1 |
| eispromtbF | eis promoter, kanamycin | Mtb complex | forward | `GCGTAACGTCACGGCGAAATTC` | amplification and sequencing | NZM 2 |
| eispromtbR | eis promoter, kanamycin | Mtb complex | reverse | `GTCAGCTCATGCAAGGTG` | amplification and sequencing | NZM 2 |
| eistbF | eis, kanamycin | Mtb complex | forward | `GTGACTGTGACCCTGTGTAG` | amplification and sequencing | NZM 2 |
| eistbR | eis, kanamycin | Mtb complex | reverse | `CACCGCGATCGAGAGTGTAC` | amplification and sequencing | NZM 2 |
| eistbseq | eis, kanamycin | Mtb complex | forward | `GCGATCTACGAGCGATGG` | sequencing | NZM 2 |
| Pri8 | gyrA (QRDR) | genus Mycobacterium | forward | `YGGTGGRTCRTTRCCYGGCGA` | amplification and sequencing | NZM 3 |
| Pri9 | gyrA (QRDR) | genus Mycobacterium | reverse | `CGCCGCGTGCTSTATGCRATG` | amplification and sequencing | NZM 3 |
| Pri10 | gyrA (QRDR) | genus Mycobacterium | reverse | `GACTCGGGCTTCCGTCCGGA` | amplification and sequencing | NZM 3 |
| panD_F | panD, pyrazinamide | Mtb complex | forward | `TCAACGGTTCCGGTCGGCTGCT` | amplification and sequencing | NZM 5 |
| panD_R | panD, pyrazinamide | Mtb complex | reverse | `TATCCGCCACTGCTGCACGACCTT` | amplification and sequencing | NZM 5 |
| gidB 5 | gidB (gid), streptomycin | Mtb complex | forward | `CGCCGAGTCGTTGTGCT` | amplification and sequencing | NZM 14 |
| gidB 3 | gidB (gid), streptomycin | Mtb complex | reverse | `AGCCTGGCCCGACCTTA` | amplification and sequencing | NZM 14 |
| gyrA-Mchim fw | gyrA | M. chimaera | forward | `GCGATTCAGACCTGCCACCG` | amplification and sequencing | NZM 15 |
| gyrA-Mchim rev | gyrA | M. chimaera | reverse | `CGGCCGTCGTAGTTAGGGATG` | amplification and sequencing | NZM 15 |
| mabs-gyrA f1 | gyrA | M. abscessus | forward | `TCGGACGTCAGAGTGGTCATGAT` | amplification and sequencing | NZM 16 |
| mabs-gyrA r1 | gyrA | M. abscessus | reverse | `CTTCGGCCAACTCACGCAGGT` | amplification and sequencing | NZM 16 |
| CME-1Y | erm(38) / erm(39) / erm(40) | M. fortuitum group | forward | `ACGTGGTGGTGGGCAAYCTG` | amplification and sequencing | NZM 18 |
| CME-2 | erm(38) / erm(39) / erm(40) | M. fortuitum group | reverse | `AATTCGAACCACGGCCACCACT` | amplification and sequencing | NZM 18 |
| rpoB-F1 | rpoB (Xpert MTB/RIF region), rifampicin | Mtb complex | forward | `CGTGGAGGCGATCACACCGCAGAC` | amplification | NZM 21 |
| rpoB-R | rpoB (Xpert MTB/RIF region), rifampicin | Mtb complex | reverse | `AGCTCCAGCCCGGCACGCTCACGT` | amplification and sequencing | NZM 21 |
| rpoB-F2 | rpoB (Xpert MTB/RIF region), rifampicin | Mtb complex | forward | `ATCAACATCCGGCCGGTGGTCGCC` | amplification and sequencing | NZM 21 |
| rpoB-R2 | rpoB (Xpert MTB/RIF region), rifampicin | Mtb complex | reverse | `TCACGTGACAGACCGCCGGGC` | amplification | NZM 21 |
| emb 151 | embB, ethambutol (alternative) | genus Mycobacterium | forward | `CGGCATGCGCCGGCTGATTC` | amplification and sequencing | NZM 28 |
| emb 131 | embB, ethambutol (alternative) | genus Mycobacterium | reverse | `TCCACAGACTGGCGTCGCTG` | amplification and sequencing | NZM 28 |
| rspA-1F | rpsA, pyrazinamide | Mtb complex | forward | `ATGCCGAGTCCCACCGTC` | amplification and sequencing | NZM 29 |
| rspA-1R | rpsA, pyrazinamide | Mtb complex | reverse | `ACCCTTGACGACCTCGATGA` | amplification and sequencing | NZM 29 |
| rspA-2F | rpsA, pyrazinamide | Mtb complex | forward | `AAACGCGCGCAGTACG` | amplification and sequencing | NZM 29 |
| rspA-2R | rpsA, pyrazinamide | Mtb complex | reverse | `GTGACCTCGTCACCAACCTG` | amplification and sequencing | NZM 29 |
| rspA-3F | rpsA, pyrazinamide | Mtb complex | forward | `GACGGTCTGGTGCATGTCT` | amplification and sequencing | NZM 29 |
| rspA-3R | rpsA, pyrazinamide | Mtb complex | reverse | `CTGCTTTTCGAATCCCTCAA` | amplification and sequencing | NZM 29 |
| rspA-4F | rpsA, pyrazinamide | Mtb complex | forward | `ATGGCTTGAGGGATTCGAAAA` | amplification and sequencing | NZM 29 |
| rspA-4R | rpsA, pyrazinamide | Mtb complex | reverse | `AGCGCTGCCGGCGAGT` | amplification and sequencing | NZM 29 |
| rspA-M5F | rpsA, pyrazinamide | Mtb complex | forward | `AAACCAACGAATGGCTTGAG` | amplification and sequencing | NZM 29 |
| rspA-M5R | rpsA, pyrazinamide | Mtb complex | reverse | `CTACTGGCCGACGACTGAT` | amplification and sequencing | NZM 29 |
| folCF1 | folC, PAS | Mtb complex | forward | `CGACGTGGCCGCCGAGCGTCG` | amplification and sequencing | NZM 30 |
| folCR1 | folC, PAS | Mtb complex | reverse | `CGGCGAACCGTCGGGGGCCCG` | amplification and sequencing | NZM 30 |
| folCF2 | folC, PAS | Mtb complex | forward | `GCGGACGCGCCCGTCGACGTG` | amplification and sequencing | NZM 30 |
| folCR2 | folC, PAS | Mtb complex | reverse | `GCGTCGCCGCCATCACCGCGC` | amplification and sequencing | NZM 30 |
| ermpromfor | erm(41) promoter region | M. abscessus | forward | `ACCGTGCAGATGGAGAAGTC` | amplification and sequencing | NZM erm41 promoter |
| ermpromrev | erm(41) promoter region | M. abscessus | reverse | `CCAGATCCACAACGAACTCC` | amplification and sequencing | NZM erm41 promoter |

## Amplicons

"Protocol" is the product length written in the protocol; "H37Rv" is the in-silico product on *M. tuberculosis* H37Rv (GCF_000195955.2), primers included, at most 3 mismatches per primer. Pairs for other species give no product there.

| Amplicon | Target | Forward + reverse | Protocol | H37Rv | Use | Sequenced with | Protocol no. |
|---|---|---|---|---|---|---|---|
| 16S_Mbak283f_Mbak264r | 16S rRNA gene (rrs), 5' part | Mbak283f + Mbak264r | – | 1026 bp | routine (protocol 3; also 32.43522) | Mbak14, possibly Mbak259r | 32.43304 |
| 16S_Mbak283f_Mbak259r | 16S rRNA gene (rrs), 5' part | Mbak283f + Mbak259r | – | 589 bp | reamplification (protocol 4) | Mbak14, possibly Mbak259r | 32.43304 |
| 16S_IMMF_IMMR2 | 16S rRNA gene (rrs), 5' part | IMM F + IMM R2 | – | 566 bp | routine (LC480 direct detection) | genus-positive products are sequenced according to 32.43103 | 32.43325 |
| 16S_Bact1098s_Bact1525a | 16S rRNA gene (rrs), 3' end | Bact-1098s + Bact-1525a | – | 446 bp | routine | amplification primers | 32.43506 |
| rpoB_2573f_3337r | rpoB | rpoB-2573f + rpoB-3337r | 764 bp | 761 bp | routine | forward primer | 32.43504 |
| hsp65_TB11_TB12w | hsp65 (groEL2) | TB-11 + TB-12w | 439 bp | 441 bp | routine | amplification primers | 32.43505 |
| hsp65_TB11_TB12 | hsp65 (groEL2) | TB-11 + TB-12 | 439 bp | 441 bp | fallback if TB-12w gives no product | amplification primers | 32.43505 |
| ITS_1f_2r | 16S-23S ITS | 16-23S-ITS-1f + 16-23S-ITS-2r | 210 bp | 221 bp | routine | amplification primers | 32.43606 |
| 23S_MCLR19F_MCLR21R | 23S rRNA gene (rrl), clarithromycin | MCLR 19 F + MCLR 21R | 836 bp | 836 bp | routine (816 bp in E. coli, 836 bp in M. abscessus) | amplification primers | 32.43605 |
| erm41_f2_r2 | erm(41) | erm41f2 + erm41r2 | 673 bp | no product | routine | forward and reverse primer | 32.43622 |
| erm41_f_r | erm(41) | erm41f + erm41r | 451 bp | no product | alternative, after consultation | forward and reverse primer | 32.43622 |
| embB_1f_1r | embB, ethambutol | embB-1f_mod2 + embB-1r | 678 bp | 678 bp | routine | forward primers | 32.43603 |
| embB_2f_2r | embB, ethambutol | embB-2f + embB-2r | 523 bp | 523 bp | routine | forward primers | 32.43603 |
| 23S_SP1F_SP2Rmod | 23S rRNA gene (rrl), clarithromycin / linezolid | SP1 F + SP2 R mod | 728 bp | 728 bp | routine | amplification primers | 32.43604 |
| ethA_f1_r1 | ethA, ethionamide | ethAf1_mod + ethAr1_mod | 839 bp | 839 bp | routine | forward and reverse primers | 32.43607 |
| ethA_f2_r2 | ethA, ethionamide | ethAf2_mod + ethAr2_mod | 834 bp | 834 bp | routine | forward and reverse primers | 32.43607 |
| ethR_F02_R07 | ethR, ethionamide | F02 + R07 | 1067 bp | 1070 bp | routine | amplification primers | 32.43608 |
| inhAprom_for_rev | inhA promoter (fabG1 upstream) and inhA gene, isoniazid | inhAprom-for + inhAprom-rev | 451 bp | 451 bp | routine (promoter) | amplification primers | 32.43609 |
| inhA_for_rev | inhA promoter (fabG1 upstream) and inhA gene, isoniazid | inhA-for + inhA-rev | 1044 bp | 1044 bp | routine (whole gene) | amplification primers | 32.43609 |
| katG_core | katG, isoniazid | katG_core_for + katG_core_rev | 697 bp | 697 bp | routine | amplification primers | 32.43610 |
| gyrA_174f1_72r2 | gyrA, fluoroquinolones | gyrA 174-f1 + gyrA 72-r2 | 1095 bp | 1095 bp | routine | amplification primers | 32.43611 |
| rpoB_TR1_TR2b | rpoB 81-bp core region, rifampicin | TR1 + TR2b | 411 bp | 411 bp | routine | amplification primers | 32.43612 |
| rpoB_Tb176f_Tb176r | rpoB 5' region, rifampicin | Tb176-f + Tb176-r | 365 bp | 365 bp | routine | amplification primers | 32.43619 |
| rrs_289_261 | 16S rRNA gene (rrs) 3' end, amikacin / capreomycin | #289 + #261 | 212 bp | 212 bp | routine | amplification primers | 32.43613 |
| thyA_F1_R1 | thyA, PAS | thyAF1 + thyAR1 | 1117 bp | 1116 bp | routine | forward and reverse primer, thyAI1 in addition | 32.43614 |
| rpsL_S13_S16 | rpsL, streptomycin | S13 + S16 | 504 bp | 504 bp | routine | amplification primers | 32.43615 |
| pncA_Ber | pncA, pyrazinamide | pncA-Ber-fw + pncA-Ber-rv | 823 bp | 823 bp | routine | forward primer, reverse after consultation | 32.43616 |
| rplC_5_3 | rplC, linezolid | rplC_5` + rplC_3` | 419 bp | 420 bp | routine | forward and reverse primer | 32.43617 |
| atpE_F_R | atpE, bedaquiline | atpE Forward + atpE Reverse | 244 bp | 244 bp | routine | forward and reverse primer | 32.43618 |
| Rv0678_F4_R1 | Rv0678 (mmpR5), bedaquiline / clofazimine | Rv0678-F4 + Rv0678-R1 | 560 bp | 560 bp | routine | forward and reverse primer | 32.43620 |
| Rv0678_F1_R2 | Rv0678 (mmpR5), bedaquiline / clofazimine | Rv0678-F1 + Rv0678-R2 | – | 485 bp | alternative primers (pairing not written in the protocol) | forward and reverse primer | 32.43620 |
| Rv0678_F2_R2 | Rv0678 (mmpR5), bedaquiline / clofazimine | Rv0678-F2 + Rv0678-R2 | – | 457 bp | alternative primers (pairing not written in the protocol) | forward and reverse primer | 32.43620 |
| tlyA_F3_R3 | tlyA, capreomycin | tlyAF3 + tlyAR3 | 981 bp | 981 bp | routine | amplification primers | 32.43621 |
| tlyA_1f_1r | tlyA, capreomycin | Mtb_tlyA-1f + Mtb_tlyA-1r | – | 585 bp | alternative primers | amplification primers | 32.43621 |
| tlyA_1fmod_1r | tlyA, capreomycin | Mtb_tlyA-1f-mod + Mtb_tlyA-1r | – | 725 bp | alternative primers | amplification primers | 32.43621 |
| tlyA_2f_2r | tlyA, capreomycin | Mtb_tlyA-2f + Mtb_tlyA-2r | – | 623 bp | alternative primers | amplification primers | 32.43621 |
| tlyA_2fmod_2rmod | tlyA, capreomycin | Mtb_tlyA-2f-mod + Mtb_tlyA-2r-mod | – | 813 bp | alternative primers | amplification primers | 32.43621 |
| ddn_Rv3547_F_R | ddn (Rv3547), delamanid | Rv3547_F + Rv3547_R | – | 846 bp | research protocol | forward primer, possibly reverse | NZM 1 |
| eisprom_F_R | eis promoter, kanamycin | eispromtbF + eispromtbR | 567 bp | 567 bp | research protocol | amplification primers | NZM 2 |
| eis_F_R | eis, kanamycin | eistbF + eistbR | 1368 bp | 1368 bp | research protocol | amplification primers, eistbseq in addition | NZM 2 |
| gyrA_Pri8_Pri9 | gyrA (QRDR) | Pri8 + Pri9 | – | 216 bp (0 + 1 mismatches) | research protocol (file "3-Mycobacterium gyrA" names Pri 8 as forward and Pri 9 as reverse primer) | forward primer, possibly reverse | NZM 3 |
| gyrA_Pri8_Pri10 | gyrA (QRDR) | Pri8 + Pri10 | – | 192 bp (0 + 3 mismatches) | research protocol (pairing not written: Pri10 binds on the same side as Pri9 on H37Rv) | forward primer, possibly reverse | NZM 3 |
| panD_F_R | panD, pyrazinamide | panD_F + panD_R | 663 bp | 663 bp | research protocol | amplification primers | NZM 5 |
| gidB_5_3 | gidB (gid), streptomycin | gidB 5 + gidB 3 | 886 bp | 892 bp | research protocol | amplification primers | NZM 14 |
| gyrA_Mchim | gyrA | gyrA-Mchim fw + gyrA-Mchim rev | 522 bp | no product | research protocol | amplification primers | NZM 15 |
| gyrA_Mabs | gyrA | mabs-gyrA f1 + mabs-gyrA r1 | 822 bp | no product | research protocol | forward or reverse primer | NZM 16 |
| erm39_CME1Y_CME2 | erm(38) / erm(39) / erm(40) | CME-1Y + CME-2 | 175 bp | no product | research protocol | forward or reverse primer | NZM 18 |
| rpoB_Xpert_F1_R | rpoB (Xpert MTB/RIF region), rifampicin | rpoB-F1 + rpoB-R | 192 bp | 192 bp (0 + 1 mismatches) | research protocol (Xpert phase 1) | rpoB-F2 and rpoB-R | NZM 21 |
| rpoB_Xpert_F2_R2 | rpoB (Xpert MTB/RIF region), rifampicin | rpoB-F2 + rpoB-R2 | 146 bp | 146 bp | research protocol (Xpert phase 2) | rpoB-F2 and rpoB-R | NZM 21 |
| rpoB_Xpert_F2_R | rpoB (Xpert MTB/RIF region), rifampicin | rpoB-F2 + rpoB-R | 164 bp | 164 bp (0 + 1 mismatches) | research protocol (main Xpert product, sequenced) | rpoB-F2 and rpoB-R | NZM 21 |
| embB_151_131 | embB, ethambutol (alternative) | emb 151 + emb 131 | – | 260 bp | research protocol (marked not finished) | emb 151, possibly emb 131 | NZM 28 |
| rpsA_1 | rpsA, pyrazinamide | rspA-1F + rspA-1R | 408 bp | 408 bp | research protocol | amplification primers | NZM 29 |
| rpsA_2 | rpsA, pyrazinamide | rspA-2F + rspA-2R | 464 bp | 464 bp | research protocol | amplification primers | NZM 29 |
| rpsA_3 | rpsA, pyrazinamide | rspA-3F + rspA-3R | 536 bp | 537 bp | research protocol | amplification primers | NZM 29 |
| rpsA_4 | rpsA, pyrazinamide | rspA-4F + rspA-4R | 241 bp | 238 bp | research protocol | amplification primers | NZM 29 |
| rpsA_M5 | rpsA, pyrazinamide | rspA-M5F + rspA-M5R | 163 bp | 163 bp | research protocol | amplification primers | NZM 29 |
| folC_F1_R1 | folC, PAS | folCF1 + folCR1 | 887 bp | 882 bp | research protocol | amplification primers | NZM 30 |
| folC_F2_R2 | folC, PAS | folCF2 + folCR2 | 1039 bp | 1038 bp | research protocol | amplification primers | NZM 30 |
| erm41prom_for_rev | erm(41) promoter region | ermpromfor + ermpromrev | 400 bp | no product | research protocol | amplification primers | NZM erm41 promoter |

Product length differs from the protocol for: rpoB_2573f_3337r (764 bp in the protocol, 761 bp on H37Rv); hsp65_TB11_TB12w (439 bp in the protocol, 441 bp on H37Rv); hsp65_TB11_TB12 (439 bp in the protocol, 441 bp on H37Rv); ITS_1f_2r (210 bp in the protocol, 221 bp on H37Rv); ethR_F02_R07 (1067 bp in the protocol, 1070 bp on H37Rv); thyA_F1_R1 (1117 bp in the protocol, 1116 bp on H37Rv); rplC_5_3 (419 bp in the protocol, 420 bp on H37Rv); gidB_5_3 (886 bp in the protocol, 892 bp on H37Rv); rpsA_3 (536 bp in the protocol, 537 bp on H37Rv); rpsA_4 (241 bp in the protocol, 238 bp on H37Rv); folC_F1_R1 (887 bp in the protocol, 882 bp on H37Rv); folC_F2_R2 (1039 bp in the protocol, 1038 bp on H37Rv). Some protocols give the length for another species (16S-23S ITS: *M. fortuitum*; 23S MCLR: *M. abscessus*).

## Which protocol says what

| Protocol | File | Primers |
|---|---|---|
| 32.43304 | 32.43304-5-Genus-Amplifikation seltener Mykobakterien.docx | Mbak283f, Mbak264r, Mbak259r, Mbak14 |
| 32.43325 | 32.43325-6-MTB-Direktnachweis (LC480-II).docx | IMM F, IMM R2 |
| 32.43506 | 32.43506-2-16S rRNA 3-End.docx | Bact-1098s, Bact-1525a |
| 32.43504 | 32.43504-5-Mykobakterien NTM rpob-Seq.docx | rpoB-2573f, rpoB-3337r |
| 32.43505 | 32.43505-4-Mykobakterien NTM 65kDa.docx | TB-11, TB-12w, TB-12 |
| 32.43606 | 32.43606-1-NTM 16S-23S-ITS-Seq. (M. fortuitum-Gruppe).doc | 16-23S-ITS-1f, 16-23S-ITS-2r |
| 32.43605 | 32.43605-2-NTM 23S CLR-Seq. Clarithromycin.doc | MCLR 19 F, MCLR 21R |
| 32.43622 | 32.43622-1-NTM erm(41)-Endpunkt-PCR.doc | erm41f2, erm41r2, erm41f, erm41r |
| 32.43603 | 32.43603-4-Mtb-Komplex embB-Seq. Ethambutol.doc | embB-1f_mod2, embB-1r, embB-2f, embB-2r |
| 32.43604 | 32.43604-2-Mtb-Komplex 23S CLR-LIN-Seq. Linezolid.doc | SP1 F, SP2 R mod |
| 32.43607 | 32.43607-2-Mtb-Komplex ethA-Seq. Ethionamid.doc | ethAf1_mod, ethAr1_mod, ethAf2_mod, ethAr2_mod |
| 32.43608 | 32.43608-2-Mtb-Komplex ethR-Seq. Ethionamid.doc | F02, R07 |
| 32.43609 | 32.43609-3-Mtb-Komplex inhA-Seq. Isoniazid.doc | inhAprom-for, inhAprom-rev, inhA-for, inhA-rev |
| 32.43610 | 32.43610-2-Mtb-Komplex katG-Seq. Isoniazid.doc | katG_core_for, katG_core_rev |
| 32.43611 | 32.43611-3-Mtb-Komplex gyrA-Seq.doc | gyrA 174-f1, gyrA 72-r2 |
| 32.43612 | 32.43612-2-Mtb-Komplex rpoB-81bp-Coreregion Rifampicin.doc | TR1, TR2b |
| 32.43619 | 32.43619-2-Mtb-Komplex rpoB-5`End-Region.doc | Tb176-f, Tb176-r |
| 32.43613 | 32.43613-2-Mtb-Komplex rrs-Seq. Amikacin Capreomycin.doc | #289, #261 |
| 32.43614 | 32.43614-3-Mtb-Komplex thyA-Seq. PAS.doc | thyAF1, thyAR1, thyAI1 |
| 32.43615 | 32.43615-3-Mtb-Komplex rpsL-Seq Streptomycin.doc | S13, S16 |
| 32.43616 | 32.43616-4-Mtb-Komplex pncA-Seq. PZA.doc | pncA-Ber-fw, pncA-Ber-rv |
| 32.43617 | 32.43617-2-Mtb-Komplex rplC-Seq. Lineozolid.doc | rplC_5`, rplC_3` |
| 32.43618 | 32.43618-2-Mtb-Komplex atpE-Seq  Diarylquinoline TMC207.doc | atpE Forward, atpE Reverse |
| 32.43620 | 32.43620-3-Mtb-Komplex Rv0678-Seq.doc | Rv0678-F4, Rv0678-R1, Rv0678-F1, Rv0678-F2, Rv0678-R2 |
| 32.43621 | 32.43621-2-Mtb-Komplex tlyA-Seq Capreomycin.doc | tlyAF3, tlyAR3, Mtb_tlyA-1f, Mtb_tlyA-1r, Mtb_tlyA-1f-mod, Mtb_tlyA-2r, Mtb_tlyA-2f, Mtb_tlyA-2r-mod, Mtb_tlyA-2f-mod, Mtb_tlyA-3f |
| NZM 1 | 1-Rv3547 (Delamanid-Resistenz).doc | Rv3547_F, Rv3547_R |
| NZM 2 | 2-MtbKomplex-eis-seq.doc | eispromtbF, eispromtbR, eistbF, eistbR, eistbseq |
| NZM 3 | 3-Primer gyrA Mycobacterium.doc | Pri8, Pri9, Pri10 |
| NZM 5 | 5-MtbKomplexl-panD-seq.doc | panD_F, panD_R |
| NZM 14 | 14-MtbKomplex-gidB-seq.doc | gidB 5, gidB 3 |
| NZM 15 | 15-gyrA-M.chimaera.doc | gyrA-Mchim fw, gyrA-Mchim rev |
| NZM 16 | 16-PCR-Mabs_gyrA.doc | mabs-gyrA f1, mabs-gyrA r1 |
| NZM 18 | 18-PCR-M_fortuitum-Komplex_Makrolid-Resistenz_erm39.doc | CME-1Y, CME-2 |
| NZM 21 | 21-32.435XX-MTB-Komplex rpob-seq aus MTB-RIF GeneXpert.doc | rpoB-F1, rpoB-R, rpoB-F2, rpoB-R2 |
| NZM 28 | 28-not finished embB-Gen, welche eine Ethambutol-alternativ.doc | emb 151, emb 131 |
| NZM 29 | 29-MtbKomplexl-rpsA-seq.doc | rspA-1F, rspA-1R, rspA-2F, rspA-2R, rspA-3F, rspA-3R, rspA-4F, rspA-4R, rspA-M5F, rspA-M5R |
| NZM 30 | 30-MtbKomplexl-folC-seq (1).doc | folCF1, folCR1, folCF2, folCR2 |
| NZM erm41 promoter | erm41_Promoterregion_EndpunktPCR-im Ordner2-MutationsSeq..doc | ermpromfor, ermpromrev |

## Other oligos in the mycobacterial protocols

64 protocols were searched (routine chapters with a mycobacterial title and the research folder). Oligos written there with 5'/3' ends that are not in `primers.tsv`: probes, real-time PCR and fragment-length assays, variants of a listed primer.

| Protocol file | Line | Sequence |
|---|---|---|
| 32.43325-6-MTB-Direktnachweis (LC480-II).docx | 6 | `AAAGCGCTTTAGCGGTGTGG` |
| 32.43325-6-MTB-Direktnachweis (LC480-II).docx | 7 | `CGACAAACCACCTACGAGCTCTTTACGC` |
| 32.43522-1-Mtb-Komplex NTM Kulturen Differenzierung (LC480-II).docx | 6 | `AAAGCGCTTTAGCGGTGTGG` |
| 32.43522-1-Mtb-Komplex NTM Kulturen Differenzierung (LC480-II).docx | 7 | `CGACAAACCACCTACGAGCTCTTTACGC` |
| 15-M.chimaera-Genom-typisierung.doc | 9 | `ACATTCACCCTGTCCATTC` |
| 15-M.chimaera-Genom-typisierung.doc | 9 | `CCTCCTTACGGAGCAGGAA` |
| 15-M.chimaera-Genom-typisierung.doc | 10 | `GCCGAACCATTTGGCGAAC` |
| 15-M.chimaera-Genom-typisierung.doc | 10 | `GGATTCGGCCGCGCAATTC` |
| 15-M.chimaera-Genom-typisierung.doc | 11 | `CATGGTTCGCCCTCTACAC` |
| 15-M.chimaera-Genom-typisierung.doc | 11 | `TAGGGGCAGGTCATCGAAG` |
| 15-M.chimaera-Genom-typisierung.doc | 12 | `GCTGAGCTACAGCCTCGAC` |
| 15-M.chimaera-Genom-typisierung.doc | 12 | `CGACGCCGATGACGTAAAC` |
| 15-M.chimaera-Genom-typisierung.doc | 13 | `TCAGGAATGGGTCCGGTTC` |
| 15-M.chimaera-Genom-typisierung.doc | 13 | `AGCTCGTGACGACGGAAAC` |
| 15-M.chimaera-Genom-typisierung.doc | 14 | `CGACCGCATCCAGAAACAG` |
| 15-M.chimaera-Genom-typisierung.doc | 14 | `GCTCTATGACGACCTCAAG` |
| 15-M.chimaera-Genom-typisierung.doc | 15 | `GTGCAGTTCAACCACGAAC` |
| 15-M.chimaera-Genom-typisierung.doc | 15 | `GGCGTTGAACACGTTGGTG` |
| 22-2-M.leprae RLEP RT PCR (LC480-II).docx | 18 | `GCAGCAGTATCGTGTTAGTGAA` |
| 22-2-M.leprae RLEP RT PCR (LC480-II).docx | 21 | `CGCTAGAAGGTTGCCGTAT` |

## Research protocols without a primer sequence

Protocols of the research folder in which no primer sequence was found (none written, or written without 5'/3' ends), so no reference is cut for them. "Primers named" are the names in the reagent table of the protocol.

| Protocol file | Primers named |
|---|---|
| 10-rpoB für M. marinum.doc | – |
| 11-RD_Denmark.doc | RD Denmark-flank F, RD Denmark-flank R |
| 12-hsp65 mac.doc | – |
| 15-PCR Amplifikation-RAPD M chimaera.doc | IS986-FP oder OPA 18 |
| 17-dnaJ1.doc | J10F, J335R |
| 23-M.chimaera HCU RT PCR (LC480-II).doc | – |
| 25-gyrA gyrB Mtb.doc | MTBU f, MTBU r |
| 6-katG-5`Seq. Mtb.doc | – |
| 7-NTM_ITS-Seq.doc | EC16S, Mb23S.44n |
| 8-Rd1mic- Mtb.doc | – |
