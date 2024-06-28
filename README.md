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
`combineVariants.workflows`|Array[String]|array of ids of producer workflows
`matchedMutect2.outputFileNamePrefix`|String|prefix of output file
`matchedCombineVariants.workflows`|Array[String]|array of ids of producer workflows
`matchedVep.tumorName`|String|Name of the tumor sample
`matchedVep.reference`|String|reference genome for input sample


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
`getFileName.jobMemory`|Int|4|memory allocated to preprocessing, in GB
`getFileName.timeout`|Int|1|timeout in hours
`getFileName.threads`|Int|1|number of cpu threads to be used
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
`variantEffectPredictor.normalName`|String?|None|Name of the normal sample
`filterMaf.freqList`|String|"$MAF_FILTERING_ROOT/TGL.frequency.20210609.annot.txt"|frequency list used in maf annotation
`filterMaf.genesToKeep`|String|"$MAF_FILTERING_ROOT/genes_to_keep.txt"|gene list in maf filtering
`filterMaf.modules`|String|"python/3.9 pandas/1.4.2 maf-filtering/2023-10-06"|module for running preprocessing
`filterMaf.jobMemory`|Int|8|memory allocated to preprocessing, in GB
`filterMaf.timeout`|Int|1|timeout in hours
`filterMaf.threads`|Int|1|number of cpu threads to be used
`matchedMutect2.filter_timeout`|Int|12|Hours before task timeout
`matchedMutect2.filter_memory`|Int|16|Memory allocated for job
`matchedMutect2.filter_filterExtraArgs`|String?|None|placehoulder for extra arguments
`matchedMutect2.mergeStats_timeout`|Int|5|Hours before task timeout
`matchedMutect2.mergeStats_memory`|Int|4|Memory allocated for job
`matchedMutect2.mergeVCFs_timeout`|Int|12|Hours before task timeout
`matchedMutect2.mergeVCFs_memory`|Int|4|Memory allocated for job
`matchedMutect2.runMutect2_timeout`|Int|24|Maximum amount of time (in hours) the task can run for.
`matchedMutect2.runMutect2_memory`|Int|32|Memory allocated to job (in GB).
`matchedMutect2.runMutect2_threads`|Int|4|Requested CPU threads
`matchedMutect2.runMutect2_mutect2ExtraArgs`|String?|None|placehoulder for extra arguments
`matchedMutect2.runMutect2_mutectTag`|String|"mutect2"|version tag for mutect
`matchedMutect2.splitStringToArray_timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`matchedMutect2.splitStringToArray_memory`|Int|1|Memory allocated to job (in GB)
`matchedMutect2.splitStringToArray_lineSeparator`|String|","|Interval group separator - these are the intervals to split by.
`matchedMutect2.pon`|File?|None|panel of normal
`matchedMutect2.ponIdx`|File?|None|index of pon
`matchedMutect2.gnomad`|File?|None|Genome Aggregation Database
`matchedMutect2.gnomadIdx`|File?|None|Index of gnomad
`matched_getFileName.jobMemory`|Int|4|memory allocated to preprocessing, in GB
`matched_getFileName.timeout`|Int|1|timeout in hours
`matched_getFileName.threads`|Int|1|number of cpu threads to be used
`matchedCombineVariants.jobMemory`|Int|24|memory allocated to preprocessing, in GB
`matchedCombineVariants.timeout`|Int|20|timeout in hours
`matchedCombineVariants.threads`|Int|8|number of cpu threads to be used
`matchedAnnotation.modules`|String|"samtools/1.9 bcftools/1.9 htslib/1.9 tabix/1.9"|module for running preprocessing
`matchedAnnotation.jobMemory`|Int|24|memory allocated to preprocessing, in GB
`matchedAnnotation.timeout`|Int|20|timeout in hours
`matchedAnnotation.threads`|Int|8|number of cpu threads to be used
`matchedVep.mergeVcfs_timeout`|Int|24|Maximum amount of time (in hours) the task can run for.
`matchedVep.mergeVcfs_threads`|Int|4|Requested CPU threads.
`matchedVep.mergeVcfs_overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`matchedVep.mergeVcfs_jobMemory`|Int|24|Memory allocated to job (in GB).
`matchedVep.mergeVcfs_extraArgs`|String?|None|Additional arguments to be passed directly to the command.
`matchedVep.mergeVcfs_modules`|String|"gatk/4.1.7.0"|Required environment modules.
`matchedVep.mergeMafs_timeout`|Int|24|Maximum amount of time (in hours) the task can run for.
`matchedVep.mergeMafs_threads`|Int|4|Requested CPU threads.
`matchedVep.mergeMafs_jobMemory`|Int|24|Memory allocated to job (in GB).
`matchedVep.mergeMafs_modules`|String|"tabix/0.2.6"|Required environment modules
`matchedVep.vcf2maf_timeout`|Int|48|Hours before task timeout
`matchedVep.vcf2maf_threads`|Int|4|Requested CPU threads
`matchedVep.vcf2maf_jobMemory`|Int|32|Memory allocated for this job (GB)
`matchedVep.vcf2maf_bufferSize`|Int|200|The buffer size
`matchedVep.vcf2maf_minHomVaf`|Float|0.7|The minimum vaf for homozygous calls
`matchedVep.vcf2maf_vepStats`|Boolean|true|If vepStats is true, remove flag '--no_stats' from vep. If vepStats is false, running vep with flag '--no_stats'
`matchedVep.vcf2maf_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`matchedVep.tumorOnlyAlign_timeout`|Int|6|Hours before task timeout
`matchedVep.tumorOnlyAlign_threads`|Int|4|Requested CPU threads
`matchedVep.tumorOnlyAlign_jobMemory`|Int|32|Memory allocated for this job (GB)
`matchedVep.tumorOnlyAlign_modules`|String|"bcftools/1.9 tabix/0.2.6"|Required environment modules
`matchedVep.tumorOnlyAlign_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`matchedVep.vep_timeout`|Int|16|Hours before task timeout
`matchedVep.vep_threads`|Int|4|Requested CPU threads
`matchedVep.vep_jobMemory`|Int|32|Memory allocated for this job (GB)
`matchedVep.vep_vepStats`|Boolean|true|If vepStats is true, remove flag '--no_stats' from vep. If vepStats is false, running vep with flag '--no_stats'
`matchedVep.vep_addParam`|String?|None|Additional vep parameters
`matchedVep.vep_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`matchedVep.subsetVcf_timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`matchedVep.subsetVcf_threads`|Int|4|Requested CPU threads.
`matchedVep.subsetVcf_jobMemory`|Int|32|Memory allocated to job (in GB).
`matchedVep.subsetVcf_modules`|String|"bcftools/1.9"|Required environment modules
`matchedVep.subsetVcf_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`matchedVep.chromosomeArray_timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`matchedVep.chromosomeArray_threads`|Int|4|Requested CPU threads.
`matchedVep.chromosomeArray_jobMemory`|Int|1|Memory allocated to job (in GB).
`matchedVep.getSampleNames_timeout`|Int|1|Hours before task timeout
`matchedVep.getSampleNames_threads`|Int|4|Requested CPU threads
`matchedVep.getSampleNames_jobMemory`|Int|1|Memory allocated for this job (GB)
`matchedVep.targetBedTask_timeout`|Int|6|Hours before task timeout
`matchedVep.targetBedTask_threads`|Int|4|Requested CPU threads
`matchedVep.targetBedTask_jobMemory`|Int|32|Memory allocated for this job (GB)
`matchedVep.targetBedTask_modules`|String|"bedtools/2.27 tabix/0.2.6"|Required environment modules
`matchedVep.targetBedTask_basename`|String|basename("~{vcfFile}",".vcf.gz")|Base name
`matchedVep.targetBed`|String?|None|Target bed file
`matchedFilterMaf.mafNormalFile`|File?|None|input file for normal sample
`matchedFilterMaf.freqList`|String|"$MAF_FILTERING_ROOT/TGL.frequency.20210609.annot.txt"|frequency list used in maf annotation
`matchedFilterMaf.genesToKeep`|String|"$MAF_FILTERING_ROOT/genes_to_keep.txt"|gene list in maf filtering
`matchedFilterMaf.modules`|String|"python/3.9 pandas/1.4.2 maf-filtering/2023-10-06"|module for running preprocessing
`matchedFilterMaf.jobMemory`|Int|8|memory allocated to preprocessing, in GB
`matchedFilterMaf.timeout`|Int|1|timeout in hours
`matchedFilterMaf.threads`|Int|1|number of cpu threads to be used


### Outputs

Output | Type | Description
---|---|---
`tumorVepVcf`|File|vep vcf for tumor sample
`tumorVepVcfIndex`|File|vep vcf index for tumor sample
`normalVepVcf`|File|vep vcf for normal sample
`normalVepVcfIndex`|File|vep vcf index for normal sample
`matchedVepVcf`|File|vep vcf for matched samples
`matchedVepVcfIndex`|File|vep vcf index for matched samples
`filteredMaf`|File?|maf file after filtering
`matchedfilteredMaf`|File?|maf file after filtering for matched maf(maf file of matched tumor/normal version)


## Commands
 This section lists command(s) run by mutect2Consensus workflow
 
 * Running mutect2Consensus
 
 ```
  
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
  ```
  
    bcftools annotate -a ~{uniqueVcf} \
   -c FMT/AD,FMT/DP ~{mergedVcf} -Oz \
   -o "~{outputPrefix}.merged.vcf.gz"
  
   tabix -p vcf "~{outputPrefix}.merged.vcf.gz"
  ```
   ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
