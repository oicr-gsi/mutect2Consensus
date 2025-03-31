# mutect2Consensus

The Mutect2Consensus workflow will process umiConsensus outputs for the tumour data through mutect2 in tumour only mode to call variants then use information from the matched normal to identify likely germline variants.

## Overview

## Dependencies

* [gatk 3.6-0](https://gatk.broadinstitute.org)
* [python 3.9](https://www.python.org/downloads/)
* [vep 105.0](https://useast.ensembl.org/info/docs/tools/vep/)
* [gatk 4.1.6.0](https://gatk.broadinstitute.org/)
* [tabix 0.2.6](https://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2/download)
* [vcf2maf 1.6](https://github.com/mskcc/vcf2maf)
* [pandas 1.4.2](https://pandas.pydata.org/)


## Usage

### Cromwell
```
java -jar cromwell.jar run mutect2Consensus.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`tumorInputGroup`|InputGroup|partitioned bam files from umiConsensus outputs for tumor sample
`intervalFile`|String|interval file to subset variant calls
`inputIntervalsToParalellizeBy`|String|intervals for parallelization
`tumorName`|String|Name of the tumor sample
`reference`|String|reference version
`gatk`|String|gatk version to be used
`filterMafFile`|Boolean|whether filter the maf file
`combineVariants.workflows`|Array[String]|array of ids of producer workflows
`normalCombineVariants.workflows`|Array[String]|array of ids of producer workflows
`somaticCombineVariants.workflows`|Array[String]|array of ids of producer workflows


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`normalInputGroup`|InputGroup?|None|partitioned bam files from umiConsensus outputs for normal sample
`normalName`|String?|None|name of the normal sample


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`mutect2.filter_timeout`|Int|12|Hours before task timeout
`mutect2.filter_memory`|Int|16|Memory allocated for job
`mutect2.filter_filterExtraArgs`|String?|None|placehoulder for extra arguments
`mutect2.mergeStats_timeout`|Int|5|Hours before task timeout
`mutect2.mergeStats_memory`|Int|4|Memory allocated for job
`mutect2.mergeVCFs_timeout`|Int|12|Hours before task timeout
`mutect2.mergeVCFs_memory`|Int|4|Memory allocated for job
`mutect2.runMutect2_timeout`|Int|24|Maximum amount of time (in hours) the task can run for.
`mutect2.runMutect2_memory`|Int|32|Memory allocated to job (in GB).
`mutect2.runMutect2_threads`|Int|4|Requested CPU threads
`mutect2.runMutect2_mutect2ExtraArgs`|String?|None|placehoulder for extra arguments
`mutect2.runMutect2_mutectTag`|String|"mutect2"|version tag for mutect
`mutect2.splitStringToArray_timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`mutect2.splitStringToArray_memory`|Int|1|Memory allocated to job (in GB)
`mutect2.splitStringToArray_lineSeparator`|String|","|Interval group separator - these are the intervals to split by.
`mutect2.normalBam`|File?|None|Input normal file (bam or sam).
`mutect2.normalBai`|File?|None|Index for noramlBam
`mutect2.pon`|File?|None|panel of normal
`mutect2.ponIdx`|File?|None|index of pon
`mutect2.gnomad`|File?|None|Genome Aggregation Database
`mutect2.gnomadIdx`|File?|None|Index of gnomad
`combineVariants.jobMemory`|Int|24|memory allocated to preprocessing, in GB
`combineVariants.timeout`|Int|20|timeout in hours
`combineVariants.threads`|Int|8|number of cpu threads to be used
`annotation.modules`|String|"samtools/1.9 bcftools/1.9 htslib/1.9 tabix/1.9"|module for running preprocessing
`annotation.jobMemory`|Int|24|memory allocated to preprocessing, in GB
`annotation.timeout`|Int|20|timeout in hours
`annotation.threads`|Int|8|number of cpu threads to be used
`tumorVep.mergeVcfs_timeout`|Int|24|Maximum amount of time (in hours) the task can run for.
`tumorVep.mergeVcfs_threads`|Int|4|Requested CPU threads.
`tumorVep.mergeVcfs_overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`tumorVep.mergeVcfs_jobMemory`|Int|24|Memory allocated to job (in GB).
`tumorVep.mergeVcfs_extraArgs`|String?|None|Additional arguments to be passed directly to the command.
`tumorVep.mergeVcfs_modules`|String|"gatk/4.1.7.0"|Required environment modules.
`tumorVep.mergeMafs_timeout`|Int|24|Maximum amount of time (in hours) the task can run for.
`tumorVep.mergeMafs_threads`|Int|4|Requested CPU threads.
`tumorVep.mergeMafs_jobMemory`|Int|24|Memory allocated to job (in GB).
`tumorVep.mergeMafs_modules`|String|"tabix/0.2.6"|Required environment modules
`tumorVep.vcf2maf_timeout`|Int|48|Hours before task timeout
`tumorVep.vcf2maf_threads`|Int|4|Requested CPU threads
`tumorVep.vcf2maf_jobMemory`|Int|32|Memory allocated for this job (GB)
`tumorVep.vcf2maf_bufferSize`|Int|200|The buffer size
`tumorVep.vcf2maf_minHomVaf`|Float|0.7|The minimum vaf for homozygous calls
`tumorVep.vcf2maf_vepStats`|Boolean|true|If vepStats is true, remove flag '--no_stats' from vep. If vepStats is false, running vep with flag '--no_stats'
`tumorVep.vcf2maf_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`tumorVep.tumorOnlyAlign_timeout`|Int|6|Hours before task timeout
`tumorVep.tumorOnlyAlign_threads`|Int|4|Requested CPU threads
`tumorVep.tumorOnlyAlign_jobMemory`|Int|32|Memory allocated for this job (GB)
`tumorVep.tumorOnlyAlign_modules`|String|"bcftools/1.9 tabix/0.2.6"|Required environment modules
`tumorVep.tumorOnlyAlign_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`tumorVep.vep_timeout`|Int|16|Hours before task timeout
`tumorVep.vep_threads`|Int|4|Requested CPU threads
`tumorVep.vep_jobMemory`|Int|32|Memory allocated for this job (GB)
`tumorVep.vep_vepStats`|Boolean|true|If vepStats is true, remove flag '--no_stats' from vep. If vepStats is false, running vep with flag '--no_stats'
`tumorVep.vep_addParam`|String?|None|Additional vep parameters
`tumorVep.vep_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`tumorVep.subsetVcf_timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`tumorVep.subsetVcf_threads`|Int|4|Requested CPU threads.
`tumorVep.subsetVcf_jobMemory`|Int|32|Memory allocated to job (in GB).
`tumorVep.subsetVcf_modules`|String|"bcftools/1.9"|Required environment modules
`tumorVep.subsetVcf_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`tumorVep.chromosomeArray_timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`tumorVep.chromosomeArray_threads`|Int|4|Requested CPU threads.
`tumorVep.chromosomeArray_jobMemory`|Int|1|Memory allocated to job (in GB).
`tumorVep.getSampleNames_timeout`|Int|1|Hours before task timeout
`tumorVep.getSampleNames_threads`|Int|4|Requested CPU threads
`tumorVep.getSampleNames_jobMemory`|Int|1|Memory allocated for this job (GB)
`tumorVep.targetBedTask_timeout`|Int|6|Hours before task timeout
`tumorVep.targetBedTask_threads`|Int|4|Requested CPU threads
`tumorVep.targetBedTask_jobMemory`|Int|32|Memory allocated for this job (GB)
`tumorVep.targetBedTask_modules`|String|"bedtools/2.27 tabix/0.2.6"|Required environment modules
`tumorVep.targetBedTask_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`tumorVep.targetBed`|String?|None|Target bed file
`tumorVep.normalName`|String?|None|Name of the normal sample
`normalMutect2.filter_timeout`|Int|12|Hours before task timeout
`normalMutect2.filter_memory`|Int|16|Memory allocated for job
`normalMutect2.filter_filterExtraArgs`|String?|None|placehoulder for extra arguments
`normalMutect2.mergeStats_timeout`|Int|5|Hours before task timeout
`normalMutect2.mergeStats_memory`|Int|4|Memory allocated for job
`normalMutect2.mergeVCFs_timeout`|Int|12|Hours before task timeout
`normalMutect2.mergeVCFs_memory`|Int|4|Memory allocated for job
`normalMutect2.runMutect2_timeout`|Int|24|Maximum amount of time (in hours) the task can run for.
`normalMutect2.runMutect2_memory`|Int|32|Memory allocated to job (in GB).
`normalMutect2.runMutect2_threads`|Int|4|Requested CPU threads
`normalMutect2.runMutect2_mutect2ExtraArgs`|String?|None|placehoulder for extra arguments
`normalMutect2.runMutect2_mutectTag`|String|"mutect2"|version tag for mutect
`normalMutect2.splitStringToArray_timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`normalMutect2.splitStringToArray_memory`|Int|1|Memory allocated to job (in GB)
`normalMutect2.splitStringToArray_lineSeparator`|String|","|Interval group separator - these are the intervals to split by.
`normalMutect2.normalBam`|File?|None|Input normal file (bam or sam).
`normalMutect2.normalBai`|File?|None|Index for noramlBam
`normalMutect2.pon`|File?|None|panel of normal
`normalMutect2.ponIdx`|File?|None|index of pon
`normalMutect2.gnomad`|File?|None|Genome Aggregation Database
`normalMutect2.gnomadIdx`|File?|None|Index of gnomad
`normalCombineVariants.jobMemory`|Int|24|memory allocated to preprocessing, in GB
`normalCombineVariants.timeout`|Int|20|timeout in hours
`normalCombineVariants.threads`|Int|8|number of cpu threads to be used
`normalAnnotation.modules`|String|"samtools/1.9 bcftools/1.9 htslib/1.9 tabix/1.9"|module for running preprocessing
`normalAnnotation.jobMemory`|Int|24|memory allocated to preprocessing, in GB
`normalAnnotation.timeout`|Int|20|timeout in hours
`normalAnnotation.threads`|Int|8|number of cpu threads to be used
`normalVep.mergeVcfs_timeout`|Int|24|Maximum amount of time (in hours) the task can run for.
`normalVep.mergeVcfs_threads`|Int|4|Requested CPU threads.
`normalVep.mergeVcfs_overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`normalVep.mergeVcfs_jobMemory`|Int|24|Memory allocated to job (in GB).
`normalVep.mergeVcfs_extraArgs`|String?|None|Additional arguments to be passed directly to the command.
`normalVep.mergeVcfs_modules`|String|"gatk/4.1.7.0"|Required environment modules.
`normalVep.mergeMafs_timeout`|Int|24|Maximum amount of time (in hours) the task can run for.
`normalVep.mergeMafs_threads`|Int|4|Requested CPU threads.
`normalVep.mergeMafs_jobMemory`|Int|24|Memory allocated to job (in GB).
`normalVep.mergeMafs_modules`|String|"tabix/0.2.6"|Required environment modules
`normalVep.vcf2maf_timeout`|Int|48|Hours before task timeout
`normalVep.vcf2maf_threads`|Int|4|Requested CPU threads
`normalVep.vcf2maf_jobMemory`|Int|32|Memory allocated for this job (GB)
`normalVep.vcf2maf_bufferSize`|Int|200|The buffer size
`normalVep.vcf2maf_minHomVaf`|Float|0.7|The minimum vaf for homozygous calls
`normalVep.vcf2maf_vepStats`|Boolean|true|If vepStats is true, remove flag '--no_stats' from vep. If vepStats is false, running vep with flag '--no_stats'
`normalVep.vcf2maf_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`normalVep.tumorOnlyAlign_timeout`|Int|6|Hours before task timeout
`normalVep.tumorOnlyAlign_threads`|Int|4|Requested CPU threads
`normalVep.tumorOnlyAlign_jobMemory`|Int|32|Memory allocated for this job (GB)
`normalVep.tumorOnlyAlign_modules`|String|"bcftools/1.9 tabix/0.2.6"|Required environment modules
`normalVep.tumorOnlyAlign_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`normalVep.vep_timeout`|Int|16|Hours before task timeout
`normalVep.vep_threads`|Int|4|Requested CPU threads
`normalVep.vep_jobMemory`|Int|32|Memory allocated for this job (GB)
`normalVep.vep_vepStats`|Boolean|true|If vepStats is true, remove flag '--no_stats' from vep. If vepStats is false, running vep with flag '--no_stats'
`normalVep.vep_addParam`|String?|None|Additional vep parameters
`normalVep.vep_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`normalVep.subsetVcf_timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`normalVep.subsetVcf_threads`|Int|4|Requested CPU threads.
`normalVep.subsetVcf_jobMemory`|Int|32|Memory allocated to job (in GB).
`normalVep.subsetVcf_modules`|String|"bcftools/1.9"|Required environment modules
`normalVep.subsetVcf_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`normalVep.chromosomeArray_timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`normalVep.chromosomeArray_threads`|Int|4|Requested CPU threads.
`normalVep.chromosomeArray_jobMemory`|Int|1|Memory allocated to job (in GB).
`normalVep.getSampleNames_timeout`|Int|1|Hours before task timeout
`normalVep.getSampleNames_threads`|Int|4|Requested CPU threads
`normalVep.getSampleNames_jobMemory`|Int|1|Memory allocated for this job (GB)
`normalVep.targetBedTask_timeout`|Int|6|Hours before task timeout
`normalVep.targetBedTask_threads`|Int|4|Requested CPU threads
`normalVep.targetBedTask_jobMemory`|Int|32|Memory allocated for this job (GB)
`normalVep.targetBedTask_modules`|String|"bedtools/2.27 tabix/0.2.6"|Required environment modules
`normalVep.targetBedTask_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`normalVep.targetBed`|String?|None|Target bed file
`normalVep.normalName`|String?|None|Name of the normal sample
`filterMaf.freqList`|String|"$MAF_FILTERING_ROOT/TGL.frequency.20210609.annot.txt"|frequency list used in maf annotation
`filterMaf.modules`|String|"python/3.9 pandas/1.4.2 maf-filtering/2024-07-10"|module for running preprocessing
`filterMaf.jobMemory`|Int|8|memory allocated to preprocessing, in GB
`filterMaf.timeout`|Int|1|timeout in hours
`filterMaf.threads`|Int|1|number of cpu threads to be used
`somaticMutect2.filter_timeout`|Int|12|Hours before task timeout
`somaticMutect2.filter_memory`|Int|16|Memory allocated for job
`somaticMutect2.filter_filterExtraArgs`|String?|None|placehoulder for extra arguments
`somaticMutect2.mergeStats_timeout`|Int|5|Hours before task timeout
`somaticMutect2.mergeStats_memory`|Int|4|Memory allocated for job
`somaticMutect2.mergeVCFs_timeout`|Int|12|Hours before task timeout
`somaticMutect2.mergeVCFs_memory`|Int|4|Memory allocated for job
`somaticMutect2.runMutect2_timeout`|Int|24|Maximum amount of time (in hours) the task can run for.
`somaticMutect2.runMutect2_memory`|Int|32|Memory allocated to job (in GB).
`somaticMutect2.runMutect2_threads`|Int|4|Requested CPU threads
`somaticMutect2.runMutect2_mutect2ExtraArgs`|String?|None|placehoulder for extra arguments
`somaticMutect2.runMutect2_mutectTag`|String|"mutect2"|version tag for mutect
`somaticMutect2.splitStringToArray_timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`somaticMutect2.splitStringToArray_memory`|Int|1|Memory allocated to job (in GB)
`somaticMutect2.splitStringToArray_lineSeparator`|String|","|Interval group separator - these are the intervals to split by.
`somaticMutect2.pon`|File?|None|panel of normal
`somaticMutect2.ponIdx`|File?|None|index of pon
`somaticMutect2.gnomad`|File?|None|Genome Aggregation Database
`somaticMutect2.gnomadIdx`|File?|None|Index of gnomad
`somaticCombineVariants.jobMemory`|Int|24|memory allocated to preprocessing, in GB
`somaticCombineVariants.timeout`|Int|20|timeout in hours
`somaticCombineVariants.threads`|Int|8|number of cpu threads to be used
`somaticAnnotation.modules`|String|"samtools/1.9 bcftools/1.9 htslib/1.9 tabix/1.9"|module for running preprocessing
`somaticAnnotation.jobMemory`|Int|24|memory allocated to preprocessing, in GB
`somaticAnnotation.timeout`|Int|20|timeout in hours
`somaticAnnotation.threads`|Int|8|number of cpu threads to be used
`somaticVep.mergeVcfs_timeout`|Int|24|Maximum amount of time (in hours) the task can run for.
`somaticVep.mergeVcfs_threads`|Int|4|Requested CPU threads.
`somaticVep.mergeVcfs_overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`somaticVep.mergeVcfs_jobMemory`|Int|24|Memory allocated to job (in GB).
`somaticVep.mergeVcfs_extraArgs`|String?|None|Additional arguments to be passed directly to the command.
`somaticVep.mergeVcfs_modules`|String|"gatk/4.1.7.0"|Required environment modules.
`somaticVep.mergeMafs_timeout`|Int|24|Maximum amount of time (in hours) the task can run for.
`somaticVep.mergeMafs_threads`|Int|4|Requested CPU threads.
`somaticVep.mergeMafs_jobMemory`|Int|24|Memory allocated to job (in GB).
`somaticVep.mergeMafs_modules`|String|"tabix/0.2.6"|Required environment modules
`somaticVep.vcf2maf_timeout`|Int|48|Hours before task timeout
`somaticVep.vcf2maf_threads`|Int|4|Requested CPU threads
`somaticVep.vcf2maf_jobMemory`|Int|32|Memory allocated for this job (GB)
`somaticVep.vcf2maf_bufferSize`|Int|200|The buffer size
`somaticVep.vcf2maf_minHomVaf`|Float|0.7|The minimum vaf for homozygous calls
`somaticVep.vcf2maf_vepStats`|Boolean|true|If vepStats is true, remove flag '--no_stats' from vep. If vepStats is false, running vep with flag '--no_stats'
`somaticVep.vcf2maf_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`somaticVep.tumorOnlyAlign_timeout`|Int|6|Hours before task timeout
`somaticVep.tumorOnlyAlign_threads`|Int|4|Requested CPU threads
`somaticVep.tumorOnlyAlign_jobMemory`|Int|32|Memory allocated for this job (GB)
`somaticVep.tumorOnlyAlign_modules`|String|"bcftools/1.9 tabix/0.2.6"|Required environment modules
`somaticVep.tumorOnlyAlign_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`somaticVep.vep_timeout`|Int|16|Hours before task timeout
`somaticVep.vep_threads`|Int|4|Requested CPU threads
`somaticVep.vep_jobMemory`|Int|32|Memory allocated for this job (GB)
`somaticVep.vep_vepStats`|Boolean|true|If vepStats is true, remove flag '--no_stats' from vep. If vepStats is false, running vep with flag '--no_stats'
`somaticVep.vep_addParam`|String?|None|Additional vep parameters
`somaticVep.vep_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`somaticVep.subsetVcf_timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`somaticVep.subsetVcf_threads`|Int|4|Requested CPU threads.
`somaticVep.subsetVcf_jobMemory`|Int|32|Memory allocated to job (in GB).
`somaticVep.subsetVcf_modules`|String|"bcftools/1.9"|Required environment modules
`somaticVep.subsetVcf_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`somaticVep.chromosomeArray_timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`somaticVep.chromosomeArray_threads`|Int|4|Requested CPU threads.
`somaticVep.chromosomeArray_jobMemory`|Int|1|Memory allocated to job (in GB).
`somaticVep.getSampleNames_timeout`|Int|1|Hours before task timeout
`somaticVep.getSampleNames_threads`|Int|4|Requested CPU threads
`somaticVep.getSampleNames_jobMemory`|Int|1|Memory allocated for this job (GB)
`somaticVep.targetBedTask_timeout`|Int|6|Hours before task timeout
`somaticVep.targetBedTask_threads`|Int|4|Requested CPU threads
`somaticVep.targetBedTask_jobMemory`|Int|32|Memory allocated for this job (GB)
`somaticVep.targetBedTask_modules`|String|"bedtools/2.27 tabix/0.2.6"|Required environment modules
`somaticVep.targetBedTask_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`somaticVep.targetBed`|String?|None|Target bed file
`somaticFilterMaf.mafNormalFile`|File?|None|input file for normal sample
`somaticFilterMaf.freqList`|String|"$MAF_FILTERING_ROOT/TGL.frequency.20210609.annot.txt"|frequency list used in maf annotation
`somaticFilterMaf.modules`|String|"python/3.9 pandas/1.4.2 maf-filtering/2024-07-10"|module for running preprocessing
`somaticFilterMaf.jobMemory`|Int|8|memory allocated to preprocessing, in GB
`somaticFilterMaf.timeout`|Int|1|timeout in hours
`somaticFilterMaf.threads`|Int|1|number of cpu threads to be used


### Outputs

Output | Type | Description | Labels
---|---|---|---
`tumorVcf`|File|vep vcf for tumor sample|vidarr_label: tumorVcf
`tumorVcfIndex`|File|vep vcf index for tumor sample|vidarr_label: tumorVcfIndex
`normalVcf`|File?|vep vcf for normal sample|vidarr_label: normalVcf
`normalVcfIndex`|File?|vep vcf index for normal sample|vidarr_label: normalVcfIndex
`somaticVcf`|File?|vep vcf for somatic samples|vidarr_label: somaticVcf
`somaticVcfIndex`|File?|vep vcf index for somatic samples|vidarr_label: somaticVcfIndex
`tumorMaf`|File?|maf file of tumor, before filtering|vidarr_label: tumorMaf
`normalMaf`|File?|maf file of normal, before filtering|vidarr_label: normalMaf
`somaticMaf`|File?|Unfiltered maf file generated from mutect2 run in somatic mode, with matched tumor and normal|vidarr_label: somaticMaf
`tumorFilteredMaf`|File?|tumour Maf with normal annotation and filtering|vidarr_label: filterredMaf
`somaticFilteredMaf`|File?|maf file after filtering for somaticMaf|vidarr_label: somaticFilterredMaf


## Commands
This section lists command(s) run by mutect2Consensus workflow

* Running mutect2Consensus

```
  python3<<CODE
  import subprocess
  import sys
  inputStrings = []
  v = "~{sep=' ' inputVcfs}"
  vcfFiles = v.split()
  w = "~{sep=' ' workflows}"
  workflowIds = w.split()
  priority = "~{priority}"
  
  if len(vcfFiles) != len(workflowIds):
      print("The arrays with input files and their respective workflow names are not of equal size!")
  else:
      for f in range(0, len(vcfFiles)):
          inputStrings.append("--variant:" + workflowIds[f] + " " + vcfFiles[f])

  javaMemory = ~{jobMemory} - 6 
  gatkCommand  = "$JAVA_ROOT/bin/java -Xmx" + str(javaMemory) + "G -jar $GATK_ROOT/GenomeAnalysisTK.jar "
  gatkCommand += "-T CombineVariants "
  gatkCommand += " ".join(inputStrings)
  gatkCommand += " -R ~{referenceFasta} "
  gatkCommand += "-o ~{outputPrefix}_combined.vcf.gz "
  gatkCommand += "-genotypeMergeOptions PRIORITIZE "
  gatkCommand += "-priority " + priority
  gatkCommand += " 2>&1"

  result_output = subprocess.run(gatkCommand, shell=True)
  sys.exit(result_output.returncode)
  CODE
```
```
  bcftools annotate -a ~{uniqueVcf} \
 -c FMT/AD,FMT/DP ~{mergedVcf} -Oz \
 -o "~{outputPrefix}.merged.vcf.gz"

 tabix -p vcf "~{outputPrefix}.merged.vcf.gz"
```
```
    python3<<CODE
    ## Adapted from https://github.com/oicr-gsi/djerba/blob/GCGI-806_v1.0.0-dev/src/lib/djerba/plugins/tar/snv_indel/plugin.py
    ## this code will filter a maf file, generated from tumor-only mutect2 calls to identify likely germline calls generated from a mutect2 calls from the matched normal
    import pandas as pd
    maf_file_path = "~{mafFile}"
    maf_normal_path = "~{mafNormalFile}"
    freq_list_path = "~{freqList}"
    output_path_prefix = "~{outputPrefix}"
    clean_columns = ["t_depth", "t_ref_count", "t_alt_count", "n_depth", "n_ref_count", "n_alt_count", "gnomAD_AF"]

    if maf_normal_path:
      df_bc = pd.read_csv(maf_normal_path,
                      sep = "\t",
                      on_bad_lines="error",
                      compression='gzip',
                      skiprows=[0])

      # Clean up df_bc if normal maf available
      for column in clean_columns:
        # Convert to numeric, setting errors='coerce' to turn non-numeric values into NaN
        df_bc[column] = pd.to_numeric(df_bc[column], errors='coerce')
        # Replace NaN with 0
        df_bc[column] = df_bc[column].fillna(0)

    df_pl = pd.read_csv(maf_file_path,
                    sep = "\t",
                    on_bad_lines="error",
                    compression='gzip',
                    skiprows=[0])
    
    # Clean up df_pl
    for column in clean_columns:
      # Convert to numeric, setting errors='coerce' to turn non-numeric values into NaN
      df_pl[column] = pd.to_numeric(df_pl[column], errors='coerce')
      # Replace NaN with 0
      df_pl[column] = df_pl[column].fillna(0)

    df_freq = pd.read_csv(freq_list_path,
                  sep = "\t")


    for row in df_pl.iterrows():
      chromosome = row[1]['Chromosome']
      start_position = row[1]['Start_Position']
      reference_allele = row[1]['Reference_Allele']
      allele = row[1]['Allele']

      # If there is normal input, annotate rows with information from the matched normal and from the frequency table
      if maf_normal_path:
        # Lookup the entry in the BC and annotate the tumour maf with
        #   n_depth, n_ref_count, n_alt_count

        row_lookup = df_bc[
                    (df_bc['Chromosome'] == chromosome) & 
                    (df_bc['Start_Position'] == start_position) &
                    (df_bc['Reference_Allele'] == reference_allele) &
                    (df_bc['Allele'] == allele)]


        # If there's only one entry, take its normal values
        if len(row_lookup) == 1:
            df_pl.at[row[0], "n_depth"] = row_lookup['n_depth'].item()
            df_pl.at[row[0], "n_ref_count"] = row_lookup['n_ref_count'].item()
            df_pl.at[row[0], "n_alt_count"] = row_lookup['n_alt_count'].item()
      
        # If the entry isn't in the table, 
        # or if there is more than one value and so you can't choose which normal values to take, 
        # set them as 0
        else:
            df_pl.at[row[0], "n_depth"] = 0
            df_pl.at[row[0], "n_ref_count"] = 0
            df_pl.at[row[0], "n_alt_count"] = 0
            
      # Lookup the entry in the frequency table and annotate the tumour maf with Freq
    
      row_lookup = df_freq[(df_freq['Start_Position'] == row[1]['Start_Position']) &
                          (df_freq['Reference_Allele'] == row[1]['Reference_Allele']) &
                          ((df_freq['Tumor_Seq_Allele'] == row[1]['Tumor_Seq_Allele1']) |
                          (df_freq['Tumor_Seq_Allele'] == row[1]['Tumor_Seq_Allele2']))]

      if len(row_lookup) > 0:
          df_pl.at[row[0], 'Freq'] = row_lookup['Freq'].item()
      else:
          df_pl.at[row[0], 'Freq'] = 0

    # Filter the maf to remove rows based on various criteria
    for row in df_pl.iterrows():
        frequency = row[1]['Freq']
        gnomAD_AF = row[1]['gnomAD_AF']
        n_alt_count = row[1]['n_alt_count']
        if  frequency > 0.1 or n_alt_count > 4 or gnomAD_AF > 0.001:
            df_pl = df_pl.drop(row[0])   

    df_pl.to_csv(output_path_prefix + '_filtered_maf.gz', sep = "\t", compression='gzip', index=False)
    CODE
  ```

 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
