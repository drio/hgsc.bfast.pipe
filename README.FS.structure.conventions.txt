
This document specifies what is the FS layout for all the next-gen sequencing
pipeline data @Ardmore.

  ANY SINGLE BIT OF DATA WE DROP ON THE FILESYSTEM HAS TO FOLLOW THE STANDARDS
  DEFINED HERE. 

  LET ME SAY IT ONE MORE TIME: ANY SINGLE BIT OF DATA WE DROP ON THE FILESYSTEM 
  HAS TO FOLLOW THE STANDARDS DEFINED HERE. 

  NOTE: I have never written a sentence in uppercase. 

Overview:
--------

Currently (work in progress to clean up)

/stornext/snfs1/next-gen/solid
  analysis/             : (2)
  capture_stats/        : (M)
  cmaps/                : references for corona_lite
  egenotype_scratch/    : (M) 
  hgsc.solid.pipeline/  : repository for SOLiD pipeline code/projects
  reference/            : (M)
  results/              : (1)
  solid_offline/        : (?) links to analysis based on SE name
  solid_runs/           : (?) links to raw data based on SE name
  solid_snps/           : (4)

/stornext/snfs4/next-gen/solid
  analysis/             : (2)
  bf.references/        : references for Bfast. One single place (this one)
  capture_stats/        : (M)
  merged_bams/          : (3)
  results/              : (1)
  solid.bf.references/  : (D)
  solid_snps/           : (4)

M : delete, clean up, not needed (we can move to another location if we 
    don't want to remove)
\d: A detailed explanation follows.

MOVING FORWARD STANDARDS:

/stornext/snfs4/next-gen/solid
  analysis/             : (2)
  merged_bams/          : (3)
  results/              : (1)
  solid_snps/           : (4)

-----------
1. RAW DATA
-----------

  /stornext/snfs#/next-gen/solid/results/INST/(*)

  INST = solid0044

  (*) = mimics the FS structure of the instrument (__for the moment__)
  Moving forward there should not be any intensity files. If any that
  are non-1kg project, it can be removed.

  The final dirs, should have only 1 or 2(if MP) primary dirs and
  only csfastas + qual + stats (because it is small and cheap).

---
SEA (sequence event analysis)
---

When SEA done, we'll see:

  -rw-------  1 p-solid next-gen 2.7K Mar 29 14:10 bf.config.yaml
  drwxrwsr-x  2 p-solid next-gen 2.1K Mar 30 04:23 cap_stats
  -rwxr-xr-x  1 p-solid next-gen  48K Mar 29 14:36 cluster_JOBS.sh
  -rw-rw-r--  1 p-solid next-gen 2.0K Mar 29 14:12 go.sh
  drwxrwsr-x  2 p-solid next-gen 2.1K Mar 29 14:12 input
  drwxrwsr-x 11 p-solid next-gen 6.1K Mar 30 04:24 lsf_logs
  -rw-rw-r--  1 p-solid next-gen  919 Mar 30 03:57 marked.stats.txt
  -rw-rw-r--  1 p-solid next-gen  931 Mar 30 03:06 metric_file.picard
  drwxrwsr-x 11 p-solid next-gen 2.1K Mar 30 03:57 output
  drwxrwsr-x  2 p-solid next-gen 2.1K Mar 29 14:35 reads
  -rw-rw-r--  1 p-solid next-gen  157 Mar 29 14:36 rg.txt
  drwxrwsr-x  2 p-solid next-gen 8.1K Mar 30 04:23 track_jobs


After cleaning:

  -rw-------  1 p-solid next-gen 2.7K Mar 29 14:10 bf.config.yaml
  drwxrwsr-x  2 p-solid next-gen 2.1K Mar 30 04:23 cap_stats
  -rwxr-xr-x  1 p-solid next-gen  48K Mar 29 14:36 cluster_JOBS.sh
    -> REMOVE
  -rw-rw-r--  1 p-solid next-gen 2.0K Mar 29 14:12 go.sh
    -> REMOVE
  drwxrwsr-x  2 p-solid next-gen 2.1K Mar 29 14:12 input
    -> Symlinks!!
  drwxrwsr-x 11 p-solid next-gen 6.1K Mar 30 04:24 lsf_logs
  -rw-rw-r--  1 p-solid next-gen  919 Mar 30 03:57 marked.stats.txt
  -rw-rw-r--  1 p-solid next-gen  931 Mar 30 03:06 metric_file.picard
    -> REMOVE
  drwxrwsr-x 11 p-solid next-gen 2.1K Mar 30 03:57 output
    -> only BAM + txts
  drwxrwsr-x  2 p-solid next-gen 2.1K Mar 29 14:35 reads
    -> REMOVE
  -rw-rw-r--  1 p-solid next-gen  157 Mar 29 14:36 rg.txt
    -> REMOVE
  drwxrwsr-x  2 p-solid next-gen 8.1K Mar 30 04:23 track_jobs
    -> REMOVE

  we'll look like:

  -rw-------  1 p-solid next-gen 2.7K Mar 29 14:10 bf.config.yaml
  drwxrwsr-x  2 p-solid next-gen 2.1K Mar 30 04:23 cap_stats
  drwxrwsr-x  2 p-solid next-gen 2.1K Mar 29 14:12 input
  drwxrwsr-x 11 p-solid next-gen 6.1K Mar 30 04:24 lsf_logs
  -rw-rw-r--  1 p-solid next-gen  919 Mar 30 03:57 marked.stats.txt
  drwxrwsr-x 11 p-solid next-gen 2.1K Mar 30 03:57 output


  For old bfast stuff at least we should have: We don't have exactly the same filesystem structure 
  or files like in the new SEA, but, for programming purposes we should try to make this old SEA
  as similar as the new ones as possible. Meaning:

  output/*.bam
         *.txt (stats)
  input/symlinks to raw data
  config 
  ...

---
3. MERGED BAMs
---

  NOTE: projects dir is redudant and should go away.
  /stornext/snfs#/next-gen/solid/merged_bams/projects/project_name/sample_id/

  + sample_id is optional for small projects where project = sample 
  + if multiple samples, sample_id is necessary

  in the directories we should have:

  1 bam     : (sort + marked)
  cap_stats : optional if not capture
  logs      : logging info from the commands necessary to the the merge

---
4. SOLiD SNPs
---
 
  /stornext/snfs#/next-gen/solid/solid_snps/project_name/(*) 

  1 txt file with the samtools output. The name should indicate clearly
  what is the sample.

