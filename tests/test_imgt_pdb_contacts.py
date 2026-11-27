"""Crystallographic MHC-peptide contact sanity checks across solved PDB structures.

Verifies that ``align_groove_to_imgt()`` aligns canonical peptide-contacting
residues onto their exact 3D-homologous coordinates in both class-specific
``fixed_seq`` and the unified 202-slot ``pan_class_fixed_seq`` across 31
well-characterized vertebrate Class I and Class II MHC molecules (human,
mouse, rat, macaque, pig, dog, bat, marsupial, and chicken).
"""

from __future__ import annotations

import pytest

from mhcseqs.imgt import (
    G_ALPHA1_INSERTION_HOTSPOTS,
    G_ALPHA2_INSERTION_HOTSPOTS,
    G_ALPHA_INSERTION_HOTSPOTS,
    G_BETA_INSERTION_HOTSPOTS,
    IMGT_PAN_CLASS_GROOVE_SLOTS,
    PAN_CLASS_HALF1_INSERTION_HOTSPOTS,
    PAN_CLASS_HALF2_INSERTION_HOTSPOTS,
    InsertionHotspot,
    align_groove_to_imgt,
    validate_imgt_frame_schemas,
)

# Verified full-groove sequences (from mhcseqs v2.24.0 / IPD-IMGT/HLA / IPD-MHC / UniProt)
# paired with their solved PDB reference and key experimentally established peptide-contact
# residues in pan-class IMGT coordinates ("HALF1"/"HALF2", slot).
PDB_STRUCTURES: tuple[tuple[str, str, str, str, str, dict[tuple[str, str], str]], ...] = (
    # ---- Human Class I ----
    (
        "HLA-A*02:01",
        "PDB 1HHK (Tax / GILGFVFTL; A/B/F pockets)",
        "I",
        "GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEA",
        "GSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQRTD",
        {
            ("HALF1", "7"): "Y",
            ("HALF1", "9"): "F",
            ("HALF1", "24"): "A",
            ("HALF1", "45"): "M",
            ("HALF1", "59"): "Y",
            ("HALF1", "63"): "E",
            ("HALF1", "66"): "K",
            ("HALF1", "67"): "V",
            ("HALF1", "70"): "H",
            ("HALF1", "77"): "D",
            ("HALF1", "80"): "T",
            ("HALF1", "81"): "L",
            ("HALF1", "84"): "Y",
            ("HALF2", "5"): "V",
            ("HALF2", "9"): "Y",
            ("HALF2", "11"): "C",
            ("HALF2", "26"): "Y",
            ("HALF2", "55"): "T",
            ("HALF2", "58"): "K",
            ("HALF2", "59"): "W",
            ("HALF2", "70"): "Y",
            ("HALF2", "74"): "C",
            ("HALF2", "81"): "Y",
        },
    ),
    (
        "HLA-B*27:05",
        "PDB 1HSA (ankylosing spondylitis; acidic P2-Arg B-pocket E45/C67)",
        "I",
        "GSHSMRYFHTSVSRPGRGEPRFITVGYVDDTLFVRFDSDAASPREEPRAPWIEQEGPEYWDRETQICKAKAQTDREDLRTLLRYYNQSEA",
        "GSHTLQNMYGCDVGPDGRLLRGYHQDAYDGKDYIALNEDLSSWTAADTAAQITQRKWEAARVAEQLRAYLEGECVEWLRRYLENGKETLQRAD",
        {
            ("HALF1", "24"): "T",
            ("HALF1", "45"): "E",
            ("HALF1", "63"): "E",
            ("HALF1", "67"): "C",
            ("HALF1", "70"): "K",
            ("HALF1", "84"): "Y",
            ("HALF2", "7"): "N",
            ("HALF2", "26"): "D",
        },
    ),
    (
        "HLA-B*57:01",
        "PDB 3UPR (HIV controller / abacavir F-pocket D114/S116)",
        "I",
        "GSHSMRYFYTAMSRPGRGEPRFIAVGYVDDTQFVRFDSDAASPRMAPRAPWIEQEGPEYWDGETRNMKASAQTYRENLRIALRYYNQSEA",
        "GSHIIQVMYGCDVGPDGRLLRGHDQSAYDGKDYIALNEDLSSWTAADTAAQITQRKWEAARVAEQLRAYLEGLCVEWLRRYLENGKETLQRAD",
        {
            ("HALF1", "45"): "M",
            ("HALF1", "67"): "M",
            ("HALF1", "74"): "Y",
            ("HALF2", "24"): "D",
            ("HALF2", "26"): "S",
            ("HALF2", "59"): "W",
        },
    ),
    (
        "HLA-C*07:01",
        "KIR2DL2/3 C1 epitope (S77/N80)",
        "I",
        "CSHSMRYFDTAVSRPGRGEPRFISVGYVDDTQFVRFDSDAASPRGEPRAPWVEQEGPEYWDRETQNYKRQAQADRVSLRNLRGYYNQSED",
        "GSHTLQRMYGCDLGPDGRLLRGYDQSAYDGKDYIALNEDLRSWTAADTAAQITQRKLEAARAAEQLRAYLEGTCVEWLRRYLENGKETLQRAE",
        {
            ("HALF1", "77"): "S",
            ("HALF1", "80"): "N",
            ("HALF1", "84"): "Y",
            ("HALF2", "59"): "L",
        },
    ),
    (
        "HLA-C*04:01",
        "PDB 1QQD (KIR2DL1 C2 epitope N77/K80)",
        "I",
        "GSHSMRYFSTSVSWPGRGEPRFIAVGYVDDTQFVRFDSDAASPRGEPREPWVEQEGPEYWDRETQKYKRQAQADRVNLRKLRGYYNQSED",
        "GSHTLQRMFGCDLGPDGRLLRGYNQFAYDGKDYIALNEDLRSWTAADTAAQITQRKWEAAREAEQRRAYLEGTCVEWLRRYLENGKETLQRAE",
        {
            ("HALF1", "77"): "N",
            ("HALF1", "80"): "K",
            ("HALF1", "84"): "Y",
            ("HALF2", "59"): "W",
        },
    ),
    (
        "HLA-E*01:01",
        "PDB 1MHE (VMAPRTLVL leader peptide; S66/T70/S143/S147)",
        "I",
        "GSHSLKYFHTSVSRPGRGEPRFISVGYVDDTQFVRFDNDAASPRMVPRAPWMEQEGSEYWDRETRSARDTAQIFRVNLRTLRGYYNQSEA",
        "GSHTLQWMHGCELGPDRRFLRGYEQFAYDGKDYLTLNEDLRSWTAVDTAAQISEQKSNDASEAEHQRAYLEDTCVEWLHKYLEKGKETLLHLE",
        {
            ("HALF1", "9"): "H",
            ("HALF1", "66"): "S",
            ("HALF1", "70"): "T",
            ("HALF1", "73"): "I",
            ("HALF2", "7"): "W",
            ("HALF2", "55"): "S",
            ("HALF2", "59"): "S",
        },
    ),
    (
        "HLA-G*01:01",
        "PDB 1YDP (non-classical; S143/C147)",
        "I",
        "GSHSMRYFSAAVSRPGRGEPRFIAMGYVDDTQFVRFDSDSACPRMEPRAPWVEQEGPEYWEEETRNTKAHAQTDRMNLQTLRGYYNQSEA",
        "SSHTLQWMIGCDLGSDGRLLRGYEQYAYDGKDYLALNEDLRSWTAADTAAQISKRKCEAANVAEQRRAYLEGTCVEWLHRYLENGKEMLQRAD",
        {
            ("HALF1", "9"): "S",
            ("HALF1", "66"): "N",
            ("HALF1", "77"): "N",
            ("HALF1", "84"): "Y",
            ("HALF2", "55"): "S",
            ("HALF2", "59"): "C",
        },
    ),
    (
        "HLA-F*01:01",
        "PDB 5IUE (open-groove Class I with R84 replacing Y84)",
        "I",
        "GSHSLRYFSTAVSRPGRGEPRYIAVEYVDDTQFLRFDSDAAIPRMEPREPWVEQEGPQYWEWTTGYAKANAQTDRVALRNLLRRYNQSEA",
        "GSHTLQGMNGCDMGPDGRLLRGYHQHAYDGKDYISLNEDLRSWTAADTVAQITQRFYEAEEYAEEFRTYLEGECLELLRRYLENGKETLQRAD",
        {
            ("HALF1", "66"): "Y",
            ("HALF1", "77"): "A",
            ("HALF1", "80"): "N",
            ("HALF1", "84"): "R",
            ("HALF2", "26"): "H",
            ("HALF2", "59"): "Y",
        },
    ),
    # ---- Non-human Class I ----
    (
        "H2-Kb",
        "PDB 1VAC (SIINFEKL; deep C-pocket for P5-Phe and F-pocket for P8-Leu)",
        "I",
        "GPHSLRYFVTAVSRPGLGEPRYMEVGYVDDTEFVRFDSDAENPRYEPRARWMEQEGPEYWERETQKAKGNEQSFRVDLRTLLGYYNQSKG",
        "GSHTIQVISGCEVGSDGRLLRGYQQYAYDGCDYIALNEDLKTWTAADMAALITKHKWEQAGEAERLRAYLEGTCVEWLRRYLKNGNATLLRTD",
        {
            ("HALF1", "70"): "N",
            ("HALF1", "73"): "S",
            ("HALF1", "74"): "F",
            ("HALF1", "77"): "D",
            ("HALF1", "80"): "T",
            ("HALF1", "81"): "L",
            ("HALF1", "84"): "Y",
            ("HALF2", "7"): "V",
            ("HALF2", "9"): "S",
            ("HALF2", "26"): "Y",
            ("HALF2", "55"): "T",
            ("HALF2", "59"): "W",
        },
    ),
    (
        "H2-Db",
        "PDB 1HOC (ASNENMETM; P5-Asn ridge W73/Q97/F116)",
        "I",
        "GPHSMRYFETAVSRPGLEEPRYISVGYVDNKEFVRFDSDAENPRYEPRAPWMEQEGPEYWERETQKAKGQEQWFRVSLRNLLGYYNQSAG",
        "GSHTLQQMSGCDLGSDWRLLRGYLQFAYEGRDYIALNEDLKTWTAADMAAQITRRKWEQSGAAEHYKAYLEGECVEWLHRYLKNGNATLLRTD",
        {
            ("HALF1", "70"): "Q",
            ("HALF1", "73"): "W",
            ("HALF1", "84"): "Y",
            ("HALF2", "7"): "Q",
            ("HALF2", "26"): "F",
        },
    ),
    (
        "Mamu-A1*001:01",
        "PDB 1ZVS (SIV Gag CM9; P2-Thr/Ser B-pocket)",
        "I",
        "GSHSMKYFYTSMSRPGRGQPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWVEQEGPEYWDRETRNMKTETQNAPVNLRTLLRYYNQSEA",
        "GSHTLQRMVGCDLGPDGRLLRGYEQYAYDGKDYIALNEDLRSWTAADVAAQNTQRKWEAADVAESMRAYLEGQCVEWLPRYLEKGKETLQRTD",
        {
            ("HALF1", "9"): "Y",
            ("HALF1", "66"): "N",
            ("HALF1", "77"): "N",
            ("HALF1", "84"): "Y",
            ("HALF2", "26"): "Y",
            ("HALF2", "59"): "W",
        },
    ),
    (
        "SLA-1*04:01",
        "PDB 3QQ3 (swine MHC-I)",
        "I",
        "GPHSLSYFYTAVSRPDRGDSRFIAVGYVDDTQFVRFDNYAPNPRMEPRVPWIQQEGQEYWDRETRNVKETAQTYGVGLNTLRGYYNQSEA",
        "GSHTLQSMYGCYLGPDGLLLHGYRQDAYDGADYIALNEDLRSWTAADMAAQITKRKWEAADEAERRRSYLQGLCVESLRRYLEMGKDTLQRAE",
        {
            ("HALF1", "9"): "Y",
            ("HALF1", "66"): "N",
            ("HALF1", "77"): "G",
            ("HALF1", "84"): "Y",
            ("HALF2", "26"): "D",
            ("HALF2", "59"): "W",
        },
    ),
    (
        "DLA-88*508:01",
        "PDB 5F1I (dog MHC-I with alpha2-helix L155 insertion at 65.1)",
        "I",
        "SHSLRYFYTSVSRPGRGDPRFIAVGYVDDTQFVRFDSDAATGRTEPRAPWVEQEGPEYWDGETRKVKETAQVYRVDLDTLRGYYNQSEA",
        "GSHTIQTMYGCDLGPGGRLLRGYRQDAYDGADYIALNEDLRSWTAADTAAQITRRKWEAAGVAELQWRNYLETTCVEWLQRYLEMGKETLLRA",
        {
            ("HALF1", "7"): "Y",
            ("HALF1", "45"): "T",
            ("HALF1", "66"): "K",
            ("HALF1", "77"): "D",
            ("HALF1", "84"): "Y",
            ("HALF2", "26"): "D",
            ("HALF2", "59"): "W",
            ("HALF2", "65"): "E",
            ("HALF2", "66"): "Q",
        },
    ),
    (
        "Gaga-BF2*21:01",
        "PDB 3BEV (promiscuous avian Class I; charge-switch R9/D24 and avian R84)",
        "I",
        "ELHTLRYIRTAMTDPGPGLPWFVDVGYVDGELFMHYNSTARRAVPRTEWIAANTDQQYWDRETQIVQGSEQINRENLDILRRRYNQTG",
        "GSHTVQWMSGCDILEDGTIRGYHQAAYDGRDFVAFDKGTMTLTAAVPEAVPTKRKWEEGGYAEGLKQYLEETCVEWLRRYVEYGKAELGRRE",
        {
            ("HALF1", "9"): "R",
            ("HALF1", "24"): "D",
            ("HALF1", "77"): "N",
            ("HALF1", "80"): "I",
            ("HALF1", "81"): "L",
            ("HALF1", "84"): "R",
            ("HALF1", "92A"): "-",
            ("HALF2", "59"): "W",
        },
    ),
    # ---- Human Class II ----
    (
        "HLA-DRA*01:01/DRB1*01:01",
        "PDB 1DLH (HA 306-318; P1 F24/I31/G86, P4 QRRAA, conserved backbone H-bonds)",
        "II",
        "IKEEHVIIQAEFYLNPDQSGEFMFDFDGDEIFHVDMAKKETVWRLEEFGRFASFEAQGALANIAVDKANLEIMTKRSNYTPIT",
        "GDTRPRFLWQLKFECHFFNGTERVRLLERCIYNQEESVRFDSDVGEYRAVTELGRPDAEYWNSQKDLLEQRRAAVDTYCRHNYGVGESFTVQR",
        {
            ("HALF1", "7"): "Q",
            ("HALF1", "26"): "F",
            ("HALF1", "33"): "I",
            ("HALF1", "56"): "A",
            ("HALF1", "61"): "S",
            ("HALF1", "70"): "N",
            ("HALF1", "77"): "N",
            ("HALF1", "84"): "R",
            ("HALF2", "7"): "L",
            ("HALF2", "9"): "F",
            ("HALF2", "22"): "L",
            ("HALF2", "24"): "E",
            ("HALF2", "55"): "D",
            ("HALF2", "59"): "W",
            ("HALF2", "66"): "Q",
            ("HALF2", "67"): "R",
            ("HALF2", "70"): "A",
            ("HALF2", "73"): "Y",
            ("HALF2", "76"): "H",
            ("HALF2", "77"): "N",
            ("HALF2", "80"): "V",
            ("HALF2", "81"): "G",
            ("HALF2", "84"): "F",
            ("HALF2", "85"): "T",
        },
    ),
    (
        "HLA-DRA*01:01/DRB1*04:01",
        "PDB 1J8H (RA shared epitope QKRAA; K71/G86)",
        "II",
        "IKEEHVIIQAEFYLNPDQSGEFMFDFDGDEIFHVDMAKKETVWRLEEFGRFASFEAQGALANIAVDKANLEIMTKRSNYTPIT",
        "GDTRPRFLEQVKHECHFFNGTERVRFLDRYFYHQEEYVRFDSDVGEYRAVTELGRPDAEYWNSQKDLLEQKRAAVDTYCRHNYGVGESFTVQR",
        {
            ("HALF1", "84"): "R",
            ("HALF2", "55"): "D",
            ("HALF2", "66"): "Q",
            ("HALF2", "67"): "K",
            ("HALF2", "70"): "A",
            ("HALF2", "81"): "G",
        },
    ),
    (
        "HLA-DRA*01:01/DRB1*03:01",
        "PDB 1A6A (CLIP; K71/R74/V86)",
        "II",
        "IKEEHVIIQAEFYLNPDQSGEFMFDFDGDEIFHVDMAKKETVWRLEEFGRFASFEAQGALANIAVDKANLEIMTKRSNYTPIT",
        "GDTRPRFLEYSTSECHFFNGTERVRYLDRYFHNQEENVRFDSDVGEFRAVTELGRPDAEYWNSQKDLLEQKRGRVDNYCRHNYGVVESFTVQR",
        {
            ("HALF1", "84"): "R",
            ("HALF2", "55"): "D",
            ("HALF2", "66"): "Q",
            ("HALF2", "67"): "K",
            ("HALF2", "70"): "R",
            ("HALF2", "81"): "V",
        },
    ),
    (
        "HLA-DRA*01:01/DRB1*15:01",
        "PDB 1BX2 (MS MBP; A71/V86)",
        "II",
        "IKEEHVIIQAEFYLNPDQSGEFMFDFDGDEIFHVDMAKKETVWRLEEFGRFASFEAQGALANIAVDKANLEIMTKRSNYTPIT",
        "GDTRPRFLWQPKRECHFFNGTERVRFLDRYFYNQEESVRFDSDVGEFRAVTELGRPDAEYWNSQKDILEQARAAVDTYCRHNYGVVESFTVQR",
        {
            ("HALF1", "84"): "R",
            ("HALF2", "55"): "D",
            ("HALF2", "66"): "Q",
            ("HALF2", "67"): "A",
            ("HALF2", "70"): "A",
            ("HALF2", "81"): "V",
        },
    ),
    (
        "HLA-DRA*01:01/DRB5*01:01",
        "PDB 1FV1 (MBP; D70/R71/G86)",
        "II",
        "IKEEHVIIQAEFYLNPDQSGEFMFDFDGDEIFHVDMAKKETVWRLEEFGRFASFEAQGALANIAVDKANLEIMTKRSNYTPIT",
        "GDTRPRFLQQDKYECHFFNGTERVRFLHRDIYNQEEDLRFDSDVGEYRAVTELGRPDAEYWNSQKDFLEDRRAAVDTYCRHNYGVGESFTVQR",
        {
            ("HALF1", "84"): "R",
            ("HALF2", "55"): "D",
            ("HALF2", "66"): "D",
            ("HALF2", "67"): "R",
            ("HALF2", "81"): "G",
        },
    ),
    (
        "HLA-DQA1*05:01/DQB1*02:01",
        "PDB 1S9V (DQ2.5 celiac gliadin; alpha-56 deletion, non-Asp beta-A57, K71)",
        "II",
        "EDIVADHVASYGVNLYQSYGPSGQYTHEFDGDEQFYVDLGRKETVWCLPVLRQFRFDPQFALTNIAVLKHNLNSLIKRSNSTAAT",
        "RDSPEDFVYQFKGMCYFTNGTERVRLVSRSIYNREEIVRFDSDVGEFRAVTLLGLPAAEYWNSQKDILERKRAAVDRVCRHNYQLELRTTLQR",
        {
            ("HALF1", "7A"): "G",
            ("HALF1", "56"): "-",
            ("HALF1", "70"): "N",
            ("HALF1", "77"): "N",
            ("HALF1", "84"): "R",
            ("HALF2", "55"): "A",
            ("HALF2", "67"): "K",
            ("HALF2", "70"): "A",
            ("HALF2", "76"): "H",
            ("HALF2", "77"): "N",
        },
    ),
    (
        "HLA-DQA1*03:01/DQB1*03:02",
        "PDB 1JK8 (DQ8 T1D insulin B9-23; non-Asp beta-A57 exposing alpha-R76)",
        "II",
        "EDIVADHVASYGVNLYQSYGPSGQYSHEFDGDEEFYVDLERKETVWQLPLFRRFRRFDPQFALTNIAVLKHNLNIVIKRSNSTAAT",
        "RDSPEDFVYQFKGMCYFTNGTERVRLVTRYIYNREEYARFDSDVGVYRAVTPLGPPAAEYWNSQKEVLERTRAELDTVCRHNYQLELRTTLQR",
        {
            ("HALF1", "7A"): "G",
            ("HALF1", "84"): "R",
            ("HALF2", "55"): "A",
            ("HALF2", "76"): "H",
            ("HALF2", "77"): "N",
        },
    ),
    (
        "HLA-DQA1*01:02/DQB1*06:02",
        "PDB 1UVQ (DQ6 T1D-protective Asp beta-D57)",
        "II",
        "EDIVADHVASCGVNLYQFYGPSGQYTHEFDGDEQFYVDLERKETAWRWPEFSKFGGFDPQGALRNMAVAKHNLNIMIKRYNSTAAT",
        "RDSPEDFVFQFKGMCYFTNGTERVRLVTRYIYNREEYARFDSDVGVYRAVTPQGRPDAEYWNSQKEVLEGTRAELDTVCRHNYEVAFRGILQR",
        {
            ("HALF1", "7A"): "G",
            ("HALF1", "84"): "R",
            ("HALF2", "55"): "D",
            ("HALF2", "76"): "H",
            ("HALF2", "77"): "N",
        },
    ),
    (
        "HLA-DPA1*01:03/DPB1*02:01",
        "PDB 3LQZ (DP2 chronic beryllium disease E69 + GGPM P1 motif)",
        "II",
        "IKADHVSTYAAFVQTHRPTGEFMFEFDEDEMFYVDLDKKETVWHLEEFGQAFSFEAQGGLANIAILNNNLNTLIQRSNHTQAT",
        "RATPENYLFQGRQECYAFNGTQRFLERYIYNREEFVRFDSDVGEFRAVTELGRPDEEYWNSQKDILEEERAVPDRMCRHNYELGGPMTLQR",
        {
            ("HALF1", "70"): "N",
            ("HALF1", "77"): "N",
            ("HALF1", "84"): "R",
            ("HALF2", "55"): "D",
            ("HALF2", "67"): "E",
            ("HALF2", "76"): "H",
            ("HALF2", "77"): "N",
            ("HALF2", "81"): "G",
            ("HALF2", "82"): "G",
            ("HALF2", "83"): "P",
            ("HALF2", "84"): "M",
        },
    ),
    (
        "HLA-DPA1*01:03/DPB1*04:01",
        "DP4 (K69 + GGPM P1 motif)",
        "II",
        "IKADHVSTYAAFVQTHRPTGEFMFEFDEDEMFYVDLDKKETVWHLEEFGQAFSFEAQGGLANIAILNNNLNTLIQRSNHTQAT",
        "RATPENYLFQGRQECYAFNGTQRFLERYIYNREEFARFDSDVGEFRAVTELGRPAAEYWNSQKDILEEKRAVPDRMCRHNYELGGPMTLQR",
        {
            ("HALF1", "84"): "R",
            ("HALF2", "55"): "A",
            ("HALF2", "67"): "K",
            ("HALF2", "81"): "G",
            ("HALF2", "82"): "G",
            ("HALF2", "83"): "P",
            ("HALF2", "84"): "M",
        },
    ),
    (
        "HLA-DPA1*01:03/DPB1*01:01",
        "DP1 (K69 + DEAV P1 motif)",
        "II",
        "IKADHVSTYAAFVQTHRPTGEFMFEFDEDEMFYVDLDKKETVWHLEEFGQAFSFEAQGGLANIAILNNNLNTLIQRSNHTQAT",
        "RATPENYVYQGRQECYAFNGTQRFLERYIYNREEYARFDSDVGEFRAVTELGRPAAEYWNSQKDILEEKRAVPDRVCRHNYELDEAVTLQR",
        {
            ("HALF1", "84"): "R",
            ("HALF2", "55"): "A",
            ("HALF2", "67"): "K",
            ("HALF2", "81"): "D",
            ("HALF2", "82"): "E",
            ("HALF2", "83"): "A",
            ("HALF2", "84"): "V",
        },
    ),
    # ---- Mouse & Rat Class II ----
    (
        "H2-IAg7",
        "PDB 1ESO (NOD T1D; H56/S57 non-Asp + 2-aa mouse I-A beta-helix deletion)",
        "II",
        "EDDIEADHVGFYGTTVYQSPGDIGQYTHEFDGDELFYVDLDKKKTVWRLPEFGQLILFEPQGGLQNIAAEKHNLGILTKRSNFTPAT",
        "GDSERHFVHQFKGECYFTNGTQRIRLVTRYIYNREEYLRFDSDVGEYRAVTELGRHSAEYYNKQYLERTRAELDTACRHNYEETEVPTSLRR",
        {
            ("HALF1", "7A"): "G",
            ("HALF1", "84"): "R",
            ("HALF2", "54"): "H",
            ("HALF2", "55"): "S",
            ("HALF2", "61B"): "-",
            ("HALF2", "62"): "-",
            ("HALF2", "76"): "H",
            ("HALF2", "77"): "N",
        },
    ),
    (
        "H2-IAb",
        "PDB 1MUJ (C57BL/6 protective P56/D57)",
        "II",
        "EDDIEADHVGTYGISVYQSPGDIGQYTFEFDGDELFYVDLDKKETVWMLPEFGQLASFDPQGGLQNIAVVKHNLGVLTKRSNSTPAT",
        "GDSERHFVYQFMGECYFTNGTQRIRYVTRYIYNREEYVRYDSDVGEHRAVTELGRPDAEYWNSQPEILERTRAELDTVCRHNYEGPETHTSLRR",
        {
            ("HALF1", "7A"): "G",
            ("HALF1", "84"): "R",
            ("HALF2", "54"): "P",
            ("HALF2", "55"): "D",
            ("HALF2", "76"): "H",
            ("HALF2", "77"): "N",
        },
    ),
    (
        "H2-IEk",
        "PDB 1IEB (mouse I-Ek; E70/K71/R72 and F86)",
        "II",
        "IKEEHTIIQAEFYLLPDKRGEFMFDFDGDEIFHVDIEKSETIWRLEEFAKFASFEAQGALANIAVDKANLDVMKERSNNTPDA",
        "VRDSRPWFLEYCKSECHFYNGTQRVRFLKRYFYNLEENLRFDSDVGEFRAVTELGRPDAENWNSQPEILDEKRAAVDTYCRHNYEIFDNFLVPR",
        {
            ("HALF1", "84"): "R",
            ("HALF2", "55"): "D",
            ("HALF2", "66"): "E",
            ("HALF2", "67"): "K",
            ("HALF2", "68"): "R",
            ("HALF2", "81"): "F",
        },
    ),
    (
        "H2-EA*k/EB*g7",
        "NOD I-Eg7 (Q70/K71/R72 and F86)",
        "II",
        "IKEEHTIIQAEFYLLPDKRGEFMFDFDGDEIFHVDIEKSETIWRLEEFAKFASFEAQGALANIAVDKANLDVMKERSNNTPDA",
        "VRDSRPWFLEYCKSECHFYNGTQRVRFLKRYFYNLEENLRFDSDVGEFRAVTELGRPDAENWNSQPEILEQKRAAVDTYCRHNYEIFDNFLVRR",
        {
            ("HALF1", "84"): "R",
            ("HALF2", "55"): "D",
            ("HALF2", "66"): "Q",
            ("HALF2", "67"): "K",
            ("HALF2", "68"): "R",
            ("HALF2", "81"): "F",
        },
    ),
    (
        "RT1-B*l",
        "Rat Lewis RT1-B (class II A-ortholog)",
        "II",
        "QDDIEADHVAAYGINMYQYYESRGQFTHEFDGDEEFYVDLDKKETIWRIPEFGQLTSFDPQGGLQNIAIIKHNLEILMKRSNSTPAV",
        "RDSPRDFVYQFKGLCYYTNGTQRIRDVIRYIYNQEEYLRYDSDVGEYRALTELGRPSAEYFNKQYLEQTRAELDTVCRHNYEGSEVRTSLRR",
        {
            ("HALF1", "7A"): "G",
            ("HALF1", "70"): "N",
            ("HALF1", "77"): "N",
            ("HALF1", "84"): "R",
            ("HALF2", "55"): "S",
            ("HALF2", "76"): "H",
            ("HALF2", "77"): "N",
        },
    ),
    (
        "RT1-D*l",
        "Rat Lewis RT1-D (class II E-ortholog)",
        "II",
        "IREEHTIIQAEFYLSPDQNGEFMFDFDGDEIFHVDIKKSETIWRLEEFAQFASFEAQGALANIAVDKANLDIMIKRSNNTPDA",
        "VRDPTPRFLEQVKGECHFYNGTQRVRFLARYIYNREEYTRFDSDVGEFRAVTELGRPSAEYYNKQKEYMEQLRATVDTACKHDYEISESFLVPR",
        {
            ("HALF1", "70"): "N",
            ("HALF1", "77"): "N",
            ("HALF1", "84"): "R",
            ("HALF2", "55"): "S",
            ("HALF2", "76"): "H",
        },
    ),
)


def _is_subseq(sub: str, full: str) -> bool:
    it = iter(full)
    return all(ch in it for ch in sub)


def test_imgt_frame_and_hotspot_schemas_validate_cleanly():
    """Issue #180: every native and pan-class hotspot schema passes strict integrity & lineage checks."""
    import mhcseqs

    validate_imgt_frame_schemas()

    all_tables = (
        G_ALPHA1_INSERTION_HOTSPOTS,
        G_ALPHA2_INSERTION_HOTSPOTS,
        G_ALPHA_INSERTION_HOTSPOTS,
        G_BETA_INSERTION_HOTSPOTS,
        *PAN_CLASS_HALF1_INSERTION_HOTSPOTS.values(),
        *PAN_CLASS_HALF2_INSERTION_HOTSPOTS.values(),
    )
    for table in all_tables:
        assert len(table) >= 2
        for hs in table:
            assert isinstance(hs, InsertionHotspot)
            assert hs.name and hs.lineage

    # Verify clean spaced names and backwards-compatible compact aliases refer to the exact same tuples
    assert mhcseqs.G_ALPHA1_POSITIONS is mhcseqs.GALPHA1_POSITIONS
    assert mhcseqs.G_ALPHA2_POSITIONS is mhcseqs.GALPHA2_POSITIONS
    assert mhcseqs.G_ALPHA_POSITIONS is mhcseqs.GALPHA_POSITIONS
    assert mhcseqs.G_BETA_POSITIONS is mhcseqs.GBETA_POSITIONS
    assert mhcseqs.IMGT_G_ALPHA1_SLOTS is mhcseqs.IMGT_GALPHA1_SLOTS
    assert mhcseqs.IMGT_G_ALPHA2_SLOTS is mhcseqs.IMGT_GALPHA2_SLOTS
    assert mhcseqs.IMGT_G_ALPHA_SLOTS is mhcseqs.IMGT_GALPHA_SLOTS
    assert mhcseqs.IMGT_G_BETA_SLOTS is mhcseqs.IMGT_GBETA_SLOTS
    assert mhcseqs.IMGT_PAN_CLASS_HALF1_SLOTS is mhcseqs.IMGT_PAN_HALF1_SLOTS
    assert mhcseqs.IMGT_PAN_CLASS_HALF2_SLOTS is mhcseqs.IMGT_PAN_HALF2_SLOTS


@pytest.mark.parametrize(
    ("allele", "pdb_note", "mhc_class", "groove1", "groove2", "expected_contacts"),
    PDB_STRUCTURES,
    ids=[row[0] for row in PDB_STRUCTURES],
)
def test_crystallographic_peptide_contacts_and_structural_invariants(
    allele: str,
    pdb_note: str,
    mhc_class: str,
    groove1: str,
    groove2: str,
    expected_contacts: dict[tuple[str, str], str],
):
    """Verify exact PDB peptide contacts, disulfide lock, class geometry, and monotonic residue order."""
    aln = align_groove_to_imgt(groove1, groove2, mhc_class=mhc_class)
    assert aln.resolved_mhc_class == mhc_class
    assert aln.alpha1 is not None and aln.alpha2 is not None
    assert len(aln.pan_class_fixed_seq) == len(IMGT_PAN_CLASS_GROOVE_SLOTS) == 202

    # 1. Strict monotonic subsequence invariant on both class-specific and pan-class frames
    assert _is_subseq(aln.alpha1.fixed_seq.replace("-", "").replace("X", ""), groove1), f"{allele} ({pdb_note}) H1 fixed scrambled"
    assert _is_subseq(aln.alpha1.pan_class_fixed_seq.replace("-", "").replace("X", ""), groove1), f"{allele} ({pdb_note}) H1 pan scrambled"
    assert _is_subseq(aln.alpha2.fixed_seq.replace("-", "").replace("X", ""), groove2), f"{allele} ({pdb_note}) H2 fixed scrambled"
    assert _is_subseq(aln.alpha2.pan_class_fixed_seq.replace("-", "").replace("X", ""), groove2), f"{allele} ({pdb_note}) H2 pan scrambled"

    # 2. Universal conserved G-domain disulfide bridge at HALF2:11 and HALF2:74
    assert aln.residue_at_imgt("HALF2", "11") == "C", f"{allele} ({pdb_note}) lost Cys11"
    assert aln.residue_at_imgt("HALF2", "74") == "C", f"{allele} ({pdb_note}) lost Cys74"

    # 3. Class-defining A-pocket / P1-end geometry at HALF1:57-60
    cap = aln.extract_pan_class_slots([("HALF1", "57"), ("HALF1", "58"), ("HALF1", "59"), ("HALF1", "60")])
    if mhc_class == "I":
        assert "-" not in cap and "X" not in cap, f"{allele} ({pdb_note}) expected closed Class I cap at HALF1:57-60, got {cap!r}"
    else:
        assert cap == "----", f"{allele} ({pdb_note}) expected open Class II P1 end at HALF1:57-60, got {cap!r}"

    # 4. Exact experimentally verified PDB peptide-contacting residues
    for (half, pos), expected_aa in expected_contacts.items():
        got = aln.residue_at_imgt(half, pos)
        assert got == expected_aa, f"{allele} ({pdb_note}) mismatch at {half}:{pos}: expected {expected_aa!r}, got {got!r}"
