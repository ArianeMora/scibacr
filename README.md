# Detecting RNA modifications in bacteria by leveraging the variation of sequence and error content
Hello all, this is possibly (definitely) my worst repo. So if this is the first time you meet my code I am sorry.

I have made this public for the ISME conference but it is all still a work in progress. As such, things aren't packaged
 nicely, and I have just put the code into functions in the `scibacr` folder. 

The workflows that I used to generate figures for the poster are in the `notebooks` folder. I generally try and make 
things reproducible while I am writing my code, but I haven't checked yet so there are probably some sneaky hard coded 
paths in there. I will be updating this later this year, but if you're really keen to use it let me know and I can 
help you out!

### Notebooks used to generate the data
The names should be self-explanatory:

notebooks/N1_SRA_data_download.html  
notebooks/N2_mapping.html  
notebooks/N3_data_QC.html  
notebooks/N4_othertools-differr.html  
notebooks/N4_quality.html  
notebooks/N5_ecoli_cVAE.html  
notebooks/N5_pseudo_cVAE.html  
notebooks/N5_xyella_cVAE.html  
notebooks/N6_tool_comparison.html  
notebooks/N7_16S.html  
notebooks/N7_16S.ipynb  

## Not reproducible yet soz
I'll be making everything into a proper package as I'm writing my thesis (which is over the next few months).  
If you want to use my code/data please let me know and I am happy to give you direction or prioritise packaging this
(I have a number of projects that I'm finishing up.)

## References
[158] I. Barbieri, T. Kouzarides, Role of rna modifications in cancer, Nature Reviews. Cancer 20 (6)
(2020) 303–322. doi:10.1038/s41568-020-0253-2.  

[159] S. Blanco, M. Frye, Role of rna methyltransferases in tissue renewal and pathology, Current
Opinion in Cell Biology 31 (2014) 1–7. doi:10.1016/j.ceb.2014.06.006.  

[160] H. Liu, O. Begik, M. C. Lucas, J. M. Ramirez, C. E. Mason, D. Wiener, S. Schwartz, J. S.
Mattick, M. A. Smith, E. M. Novoa, Accurate detection of m6a rna modifications in native rna sequences,
Nature Communications 10 (1) (2019) 4079. doi:10.1038/s41467-019-11713-9.
BIBLIOGRAPHY 108  

[161] H. Zhang, X. Shi, T. Huang, X. Zhao, W. Chen, N. Gu, R. Zhang, Dynamic landscape and
evolution of m6a methylation in human, Nucleic Acids Research 48 (11) (2020) 6251–6264.
doi:10.1093/nar/gkaa347.  

[162] P. N. Pratanwanich, F. Yao, Y. Chen, C. W. Q. Koh, Y. K. Wan, C. Hendra, P. Poon, Y. T.
Goh, P. M. L. Yap, J. Y. Chooi, W. J. Chng, S. B. Ng, A. Thiery, W. S. S. Goh, J. G¨oke,
Identification of differential rna modifications from nanopore direct rna sequencing with xpore,
Nature Biotechnology 39 (11) (2021) 1394–1402. doi:10.1038/s41587-021-00949-w.  

[163] Y. Fu, D. Dominissini, G. Rechavi, C. He, Gene expression regulation mediated through
reversible m6a rna methylation, Nature Reviews. Genetics 15 (5) (2014) 293–306. doi:
10.1038/nrg3724.  

[164] R. Garc´ıa-V´ılchez, A. Sevilla, S. Blanco, Post-transcriptional regulation by cytosine-5 methylation
of rna, Biochimica Et Biophysica Acta. Gene Regulatory Mechanisms 1862 (3) (2019)
240–252. doi:10.1016/j.bbagrm.2018.12.003.  

[165] T. Huang, W. Chen, J. Liu, N. Gu, R. Zhang, Genome-wide identification of mrna 5-
methylcytosine in mammals, Nature Structural and Molecular Biology 26 (5) (2019) 380–388.
doi:10.1038/s41594-019-0218-x.  

[166] X. Chen, A. Li, B.-F. Sun, Y. Yang, Y.-N. Han, X. Yuan, R.-X. Chen, W.-S. Wei, Y. Liu, C.-C.
Gao, Y.-S. Chen, M. Zhang, X.-D. Ma, Z.-W. Liu, J.-H. Luo, C. Lyu, H.-L. Wang, J. Ma,
Y.-L. Zhao, F.-J. Zhou, Y. Huang, D. Xie, Y.-G. Yang, 5-methylcytosine promotes pathogenesis
of bladder cancer through stabilizing mrnas, Nature Cell Biology 21 (8) (2019) 978–990.
doi:10.1038/s41556-019-0361-y.  

[167] A. F. Lovejoy, D. P. Riordan, P. O. Brown, Transcriptome-wide mapping of pseudouridines:
pseudouridine synthases modify specific mrnas in s. cerevisiae, PloS One 9 (10) (2014) e110799.
doi:10.1371/journal.pone.0110799.  

[168] T. M. Carlile, M. F. Rojas-Duran, B. Zinshteyn, H. Shin, K. Bartoli, W. Gilbert, Pseudouridine
profiling reveals regulated mrna pseudouridylation in yeast and human cells, Nature (2014).
doi:10.1038/nature13802.  

[169] D. R. Garalde, E. A. Snell, D. Jachimowicz, B. Sipos, J. H. Lloyd, M. Bruce, N. Pantic,
T. Admassu, P. James, A. Warland, M. Jordan, J. Ciccone, S. Serra, J. Keenan, S. Martin,
L. McNeill, E. J. Wallace, L. Jayasinghe, C. Wright, J. Blasco, S. Young, D. Brocklebank,
S. Juul, J. Clarke, A. J. Heron, D. J. Turner, Highly parallel direct rna sequencing on an array of
nanopores, Nature Methods 15 (3) (2018) 201–206. doi:10.1038/nmeth.4577.  

[170] L. Xu, M. Seki, Recent advances in the detection of base modifications using the
nanopore sequencer, Journal of Human Genetics 65 (1) (2020) 25–33. doi:10.1038/
s10038-019-0679-0.
BIBLIOGRAPHY 109  

[171] M. Stoiber, J. Quick, R. Egan, J. E. Lee, S. Celniker, R. K. Neely, N. Loman, L. A. Pennacchio,
J. Brown, De novo identification of dna modifications enabled by genome-guided nanopore
signal processing (2017) 094672doi:10.1101/094672.
URL https://www.biorxiv.org/content/10.1101/094672v2  

[172] S. Schwartz, M. R. Mumbach, M. Jovanovic, T. Wang, K. Maciag, G. G. Bushkin, P. Mertins,
D. Ter-Ovanesyan, N. Habib, D. Cacchiarelli, N. E. Sanjana, E. Freinkman, M. E. Pacold,
R. Satija, T. S. Mikkelsen, N. Hacohen, F. Zhang, S. A. Carr, E. S. Lander, A. Regev, Perturbation
of m6a writers reveals two distinct classes of mrna methylation at internal and 5’ sites, Cell
Reports 8 (1) (2014) 284–296. doi:10.1016/j.celrep.2014.05.048.  

[173] M. Furlan, A. Delgado-Tejedor, L. Mulroney, M. Pelizzola, E. M. Novoa, T. Leonardi, Computational
methods for rna modification detection from nanopore direct rna sequencing data, RNA
biology 18 (sup1) (2021) 31–40. doi:10.1080/15476286.2021.1978215.  

[174] A. Leger, P. P. Amaral, L. Pandolfini, C. Capitanchik, F. Capraro, V. Miano, V. Migliori,
P. Toolan-Kerr, T. Sideri, A. J. Enright, K. Tzelepis, F. J. van Werven, N. M. Luscombe,
I. Barbieri, J. Ule, T. Fitzgerald, E. Birney, T. Leonardi, T. Kouzarides, Rna modifications
detection by comparative nanopore direct rna sequencing, Nature Communications 12 (11)
(2021) 7198. doi:10.1038/s41467-021-27393-3.  

[175] L. Antoine, R. Bahena-Ceron, H. Devi Bunwaree, M. Gobry, V. Loegler, P. Romby, S. Marzi,
Rna modifications in pathogenic bacteria: Impact on host adaptation and virulence, Genes
12 (88) (2021) 1125. doi:10.3390/genes12081125.  

[176] Y. Fang, A. Changavi, M. Yang, L. Sun, A. Zhang, D. Sun, Z. Sun, B. Zhang, M.-Q. Xu,
Nanopore whole transcriptome analysis and pathogen surveillance by a novel solid-phase
catalysis approach, Advanced Science 9 (3) (2022) 2103373. doi:10.1002/advs.202103373.  

[177] H. Zhang, H. Zhong, X. Wang, S. Zhang, X. Shao, H. Hu, Z. Yu, Z. Cai, X. Chen, Y. Xia,
Use of nad tagseq ii to identify growth phase-dependent alterations in e. coli rna nad+ capping,
Proceedings of the National Academy of Sciences of the United States of America 118 (14)
(2021) e2026183118. doi:10.1073/pnas.2026183118.