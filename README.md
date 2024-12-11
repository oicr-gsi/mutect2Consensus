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
`normalInputGroup`|InputGroup|partitioned bam files from umiConsensus outputs for normal sample
`intervalFile`|String|interval file to subset variant calls
`inputIntervalsToParalellizeBy`|String|intervals for parallelization
`tumorName`|String|Name of the tumor sample
`normalName`|String|name of the normal sample
`reference`|String|reference version
`gatk`|String|gatk version to be used
`filterMafFile`|Boolean|whether filter the maf file
`combineVariants.workflows`|Array[String]|array of ids of producer workflows
`somaticCombineVariants.workflows`|Array[String]|array of ids of producer workflows


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---


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
`variantEffectPredictor.mergeVcfs_timeout`|Int|24|Maximum amount of time (in hours) the task can run for.
`variantEffectPredictor.mergeVcfs_threads`|Int|4|Requested CPU threads.
`variantEffectPredictor.mergeVcfs_overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`variantEffectPredictor.mergeVcfs_jobMemory`|Int|24|Memory allocated to job (in GB).
`variantEffectPredictor.mergeVcfs_extraArgs`|String?|None|Additional arguments to be passed directly to the command.
`variantEffectPredictor.mergeVcfs_modules`|String|"gatk/4.1.7.0"|Required environment modules.
`variantEffectPredictor.mergeMafs_timeout`|Int|24|Maximum amount of time (in hours) the task can run for.
`variantEffectPredictor.mergeMafs_threads`|Int|4|Requested CPU threads.
`variantEffectPredictor.mergeMafs_jobMemory`|Int|24|Memory allocated to job (in GB).
`variantEffectPredictor.mergeMafs_modules`|String|"tabix/0.2.6"|Required environment modules
`variantEffectPredictor.vcf2maf_timeout`|Int|48|Hours before task timeout
`variantEffectPredictor.vcf2maf_threads`|Int|4|Requested CPU threads
`variantEffectPredictor.vcf2maf_jobMemory`|Int|32|Memory allocated for this job (GB)
`variantEffectPredictor.vcf2maf_bufferSize`|Int|200|The buffer size
`variantEffectPredictor.vcf2maf_minHomVaf`|Float|0.7|The minimum vaf for homozygous calls
`variantEffectPredictor.vcf2maf_vepStats`|Boolean|true|If vepStats is true, remove flag '--no_stats' from vep. If vepStats is false, running vep with flag '--no_stats'
`variantEffectPredictor.vcf2maf_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`variantEffectPredictor.tumorOnlyAlign_timeout`|Int|6|Hours before task timeout
`variantEffectPredictor.tumorOnlyAlign_threads`|Int|4|Requested CPU threads
`variantEffectPredictor.tumorOnlyAlign_jobMemory`|Int|32|Memory allocated for this job (GB)
`variantEffectPredictor.tumorOnlyAlign_modules`|String|"bcftools/1.9 tabix/0.2.6"|Required environment modules
`variantEffectPredictor.tumorOnlyAlign_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`variantEffectPredictor.vep_timeout`|Int|16|Hours before task timeout
`variantEffectPredictor.vep_threads`|Int|4|Requested CPU threads
`variantEffectPredictor.vep_jobMemory`|Int|32|Memory allocated for this job (GB)
`variantEffectPredictor.vep_vepStats`|Boolean|true|If vepStats is true, remove flag '--no_stats' from vep. If vepStats is false, running vep with flag '--no_stats'
`variantEffectPredictor.vep_addParam`|String?|None|Additional vep parameters
`variantEffectPredictor.vep_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`variantEffectPredictor.subsetVcf_timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`variantEffectPredictor.subsetVcf_threads`|Int|4|Requested CPU threads.
`variantEffectPredictor.subsetVcf_jobMemory`|Int|32|Memory allocated to job (in GB).
`variantEffectPredictor.subsetVcf_modules`|String|"bcftools/1.9"|Required environment modules
`variantEffectPredictor.subsetVcf_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`variantEffectPredictor.chromosomeArray_timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`variantEffectPredictor.chromosomeArray_threads`|Int|4|Requested CPU threads.
`variantEffectPredictor.chromosomeArray_jobMemory`|Int|1|Memory allocated to job (in GB).
`variantEffectPredictor.getSampleNames_timeout`|Int|1|Hours before task timeout
`variantEffectPredictor.getSampleNames_threads`|Int|4|Requested CPU threads
`variantEffectPredictor.getSampleNames_jobMemory`|Int|1|Memory allocated for this job (GB)
`variantEffectPredictor.targetBedTask_timeout`|Int|6|Hours before task timeout
`variantEffectPredictor.targetBedTask_threads`|Int|4|Requested CPU threads
`variantEffectPredictor.targetBedTask_jobMemory`|Int|32|Memory allocated for this job (GB)
`variantEffectPredictor.targetBedTask_modules`|String|"bedtools/2.27 tabix/0.2.6"|Required environment modules
`variantEffectPredictor.targetBedTask_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`variantEffectPredictor.targetBed`|String?|None|Target bed file
`variantEffectPredictor.normalName`|String?|None|Name of the normal sample
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

Output | Type | Description
---|---|---
`tumorVcf`|File|{'description': 'vep vcf for tumor sample', 'vidarr_label': 'tumorVcf'}
`tumorVcfIndex`|File|{'description': 'vep vcf index for tumor sample', 'vidarr_label': 'tumorVcfIndex'}
`normalVcf`|File|{'description': 'vep vcf for normal sample', 'vidarr_label': 'normalVcf'}
`normalVcfIndex`|File|{'description': 'vep vcf index for normal sample', 'vidarr_label': 'normalVcfIndex'}
`somaticVcf`|File|{'description': 'vep vcf for somatic samples', 'vidarr_label': 'somaticVcf'}
`somaticVcfIndex`|File|{'description': 'vep vcf index for somatic samples', 'vidarr_label': 'somaticVcfIndex'}
`tumorMaf`|File?|{'description': 'maf file of tumor, before filtering', 'vidarr_label': 'tumorMaf'}
`normalMaf`|File?|{'description': 'maf file of normal, before filtering', 'vidarr_label': 'normalMaf'}
`somaticMaf`|File?|{'description': 'Unfiltered maf file generated from mutect2 run in somatic mode, with matched tumor and normal', 'vidarr_label': 'somaticMaf'}
`tumorFilteredMaf`|File?|{'description': 'tumour Maf with normal annotation and filtering', 'vidarr_label': 'filterredMaf'}
`somaticFilteredMaf`|File?|{'description': 'maf file after filtering for somaticMaf', 'vidarr_label': 'somaticFilterredMaf'}


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
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
