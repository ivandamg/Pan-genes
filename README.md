# Pan-genome
Identification of the core-genome of a group of strains.



# 1. Calculate the core-genes of a group of strains

- Re anotate the genome assemblies with Prokka-1.12.
  - In local 
  
        for i in $(ls *NCBI.fna); do echo $i ; ~/software/prokka-1.12/prokka/bin/prokka --outdir Annotation_$(echo $i | cut -d'.' -f1) --genus Vibrio --species cholerae --strain $(echo $i | cut -d'_' -f1) --locustag VC-$(echo $i | cut -d'_' -f1) --prefix $(echo $i | cut -d'.' -f1)_Prokka --rfam --usegenus $i ; done

  - In cluster

        for i in $(ls Vibrio*.fna); do echo $i; bsub -q normal -L /bin/bash -J $(echo $i | cut -d'_' -f3) -u ivan.mateusgonzalez@epfl.ch -n 8 -R "rusage[mem=12000]" -M 12000000 -N "module add UHTS/Analysis/prokka/1.12; module add UHTS/Analysis/rnammer/1.2; module add UHTS/Analysis/LMAT/1.2.6; module add SequenceAnalysis/HMM-Profile/hmmer/3.1b2; prokka --outdir Annotation_$(echo $i | cut -d'_' -f3) --genus Vibrio --species cholerae --strain $(echo $i | cut -d'_' -f3) --cpus 8 --locustag VC_$(echo $i | cut -d'_' -f3) --prefix $(echo $i | cut -d'_' -f1,2,3)_Prokka --rfam --usegenus $i "; done

# 2. Iterative blast between all the strains
Filter blastp 80% identity.  difference in length of the sequence max. 20%

    # Rename files names Vibrio_Cholera_STAIN_BLA.fna

    # create db
    
    a=0;for i in $(ls *.faa); do echo $(echo $i | cut -d'_' -f3) ;makeblastdb -dbtype prot -in $i -parse_seqids -out db_prot_genomes/$(echo $i | cut -d'_' -f3)_db ; done

    # Sequencial blast

    # CHECK order
      ls -1 *.faa
      
    # LOOP OVER FILES TO ADD NEW HOMOLOGS
    
    # modify name of first strain to ->
    cat Vibrio_cholerae_YB1A01_Prokka.faa > Core_YB1A01.faa

    for OTHER in $(ls db_prot_genomes/*_db.pin) ; do echo $(echo $OTHER | cut -d'/' -f2 | cut -d'_' -f1) ;

    for i in $(ls db_prot_genomes/*_db.pin) ; do echo $(echo $i | cut -d'.' -f1) ;

    blastp -db $(echo $i | cut -d'.' -f1) -outfmt 6 -evalue 1e-8 -show_gis -num_alignments 1 -max_hsps 20 -num_threads 30 -out db_prot_genomes/blastProt_$(echo $i | cut -d'/' -f2 | cut -d'_' -f1)_in_$(echo $OTHER | cut -d'/' -f2 | cut -d'_' -f1)core.xml -query Core_$(echo $OTHER | cut -d'/' -f2 | cut -d'_' -f1).faa

    # Filter blastp 80% identity aminoacid sequence <20%
    cat db_prot_genomes/blastProt_$(echo $i | cut -d'/' -f2 | cut -d'_' -f1)_in_$(echo $OTHER | cut -d'/' -f2 | cut -d'_' -f1)core.xml | awk '$3 > 79.9' | awk '!seen[$1]++'  | awk '!seen[$2]++'  > $(echo $i | cut -d'/' -f2 | cut -d'_' -f1)_core.tmp

    # Filter core proteins of all strains
    awk '{a=$8-$7;print $0,a;}' $(echo $i | cut -d'/' -f2 | cut -d'_' -f1)_core.tmp > $(echo $i | cut -d'/' -f2 | cut -d'_' -f1)_core.tmp2
    awk '{a=$10-$9;print $0,a;}' $(echo $i | cut -d'/' -f2 | cut -d'_' -f1)_core.tmp2 > $(echo $i | cut -d'/' -f2 | cut -d'_' -f1)_core.tmp3
    awk '{a=$13/$14;print $0,a;}' $(echo $i | cut -d'/' -f2 | cut -d'_' -f1)_core.tmp3 > $(echo $i | cut -d'/' -f2 | cut -d'_' -f1)_core.tmp4
    cat $(echo $i | cut -d'/' -f2 | cut -d'_' -f1)_core.tmp4 | awk '$15 > 0.8' | awk '$15 < 1.2' > $(echo $i | cut -d'/' -f2 | cut -d'_' -f1)_core.tmp5

    # extract names of sequences.

    cat $(echo $i | cut -d'/' -f2 | cut -d'_' -f1)_core.tmp5 | cut -f2 > $(echo $i | cut -d'/' -f2 | cut -d'_' -f1)_core.txt

    cat $(echo $i | cut -d'/' -f2 | cut -d'_' -f1)_core.tmp5 | cut -f1 > $(echo $OTHER | cut -d'/' -f2 | cut -d'_' -f1)_core.txt

    # Core proteins of N16961

    xargs faidx -d ' ' Vibrio_cholerae_$(echo $i | cut -d'/' -f2 | cut -d'_' -f1)_Prokka.faa < $(echo $i | cut -d'/' -f2 | cut -d'_' -f1)_core.txt > Core_$(echo $i | cut -d'/' -f2 | cut -d'_' -f1).faa

    xargs faidx -d ' ' Vibrio_cholerae_$(echo $OTHER | cut -d'/' -f2 | cut -d'_' -f1)_Prokka.faa < $(echo $OTHER | cut -d'/' -f2 | cut -d'_' -f1)_core.txt > Core_$(echo $OTHER | cut -d'/' -f2 | cut -d'_' -f1).faa

    # delete temp data
    rm *tmp
    rm *tmp*
    rm *core.txt; done;done

Check that all the strains have same number of proteins
  
     for i in $(ls Core*.faa); do echo $i; cat $i | grep ">" | wc -l; done


#3. Phylogeny construction based on core-genes.

- Concatenate the sequences (no sequence ID.) and put all the info in a single file

      for i in $(ls Core*.faa); do echo $i; cat $i | grep -v ">" ; done | sed 's/Core/>Core/g' > AllStrains_CoreGenes.faa

- Align with MAFFT online 
https://mafft.cbrc.jp/alignment/server/cgi-bin/mafft5-lsf.cgi

- Discard poorly aligned regions with BMGE

      java -jar ~/software/BMGE-1.12/BMGE.jar -i ALLStrains_CoreGenes_MAFFT.pir -t AA -o ALLStrains_CoreGenes_MAFFT.phy

- Phyml Phylogeny

      bsub -q long -L /bin/bash -J Core-Genomes_PHY -u ivan.mateusgonzalez@epfl.ch -N -R "rusage[mem=60000]" -M 60000000 "module add Phylogeny/prottest/3.4.1; PhyML_3.0_linux64 -i ALLStrains_CoreGenes_MAFFT.phy -d aa -b -1 -m WAG -c 4 -a e"

# 4. Analysis of Core-genes, Variable-genes and Pan-genes for populations of n individuals.



