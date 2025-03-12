version 1.0

import "imports/pull_mutect2.wdl" as mutect2
import "imports/pull_variantEffectPredictor.wdl" as vep


struct BamAndBamIndex {
  File bam
  File bamIndex
}

struct InputGroup {
  String outputFileNamePrefix
  BamAndBamIndex dcsScBamAndIndex
  BamAndBamIndex sscsScBamAndIndex
  BamAndBamIndex allUniqueBamAndIndex
}

struct GenomeResources {
  String inputRefFasta
  String combineVariants_modules
}

workflow mutect2Consensus {
  input {
    InputGroup tumorInputGroup
    InputGroup? normalInputGroup
    String intervalFile
    String inputIntervalsToParalellizeBy
    String tumorName
    String? normalName
    String reference
    String gatk
    Boolean filterMafFile
  }

  parameter_meta {
    tumorInputGroup: "partitioned bam files from umiConsensus outputs for tumor sample"
    normalInputGroup: "partitioned bam files from umiConsensus outputs for normal sample"
    intervalFile: "interval file to subset variant calls"
    inputIntervalsToParalellizeBy: "intervals for parallelization"
    tumorName: "Name of the tumor sample"
    normalName: "name of the normal sample"
    reference: "reference version"
    gatk: "gatk version to be used"
    filterMafFile: "whether filter the maf file"
  }

 Map[String,GenomeResources] resources = {
    "hg19": {
      "inputRefFasta": "$HG19_ROOT/hg19_random.fa",
      "combineVariants_modules": "gatk/3.6-0 tabix/0.2.6 hg19/p13"
      },
    "hg38": {
      "inputRefFasta": "$HG38_ROOT/hg38_random.fa",
      "combineVariants_modules": "gatk/3.6-0 tabix/0.2.6 hg38/p12"
      }
  }

  Array[BamAndBamIndex]tumor_partitionedBams = [tumorInputGroup.dcsScBamAndIndex, tumorInputGroup.sscsScBamAndIndex, tumorInputGroup.allUniqueBamAndIndex]
  scatter ( bamAndIndex in tumor_partitionedBams ) {
    call mutect2.mutect2 {
      input:
        tumorBam = bamAndIndex.bam,
        tumorBai = bamAndIndex.bamIndex,
        intervalFile = intervalFile,
        intervalsToParallelizeBy = inputIntervalsToParalellizeBy,
        reference = reference,
        gatk = gatk,
        outputFileNamePrefix = tumorInputGroup.outputFileNamePrefix
    }
  }
  Array[File] mutect2_tumor_VcfFiles = mutect2.filteredVcfFile
  Array[File] mutect2_tumor_VcfIndexes = mutect2.filteredVcfIndex

  call combineVariants {
    input: 
      inputVcfs = [mutect2_tumor_VcfFiles[0], mutect2_tumor_VcfFiles[1]],
      inputIndexes = [mutect2_tumor_VcfIndexes[0],mutect2_tumor_VcfIndexes[1]],
      priority = "mutect2-dcsSc,mutect2-sscsSc",
      outputPrefix = tumorInputGroup.outputFileNamePrefix,
      referenceFasta = resources[reference].inputRefFasta,
      modules = resources[reference].combineVariants_modules
  }

  call annotation {
    input: 
      uniqueVcf = mutect2_tumor_VcfFiles[2],
      uniqueVcfIndex = mutect2_tumor_VcfIndexes[2],
      mergedVcf = combineVariants.combinedVcf,
      mergedVcfIndex = combineVariants.combinedIndex,
      outputPrefix = tumorInputGroup.outputFileNamePrefix
  }

  call vep.variantEffectPredictor as tumorVep {
    input: 
      vcfFile = annotation.annotatedCombinedVcf,
      vcfIndex = annotation.annotatedCombinedIndex,
      toMAF = true,
      onlyTumor = true,
      tumorOnlyAlign_updateTagValue = true,
      vcf2maf_retainInfoProvided = true,
      tumorName = tumorInputGroup.outputFileNamePrefix,
      reference = reference
  }

  File? tumor_Maf = tumorVep.outputMaf

  if (defined(normalInputGroup)) {
    InputGroup normal = select_first([normalInputGroup])
    
    call mutect2.mutect2 as normalMutect2_dcs {
      input:
        tumorBam = normal.dcsScBamAndIndex.bam,
        tumorBai = normal.dcsScBamAndIndex.bamIndex,
        intervalFile = intervalFile,
        intervalsToParallelizeBy = inputIntervalsToParalellizeBy,
        reference = reference,
        gatk = gatk,
        outputFileNamePrefix = normalName + "_dcs"
    }
    
    call mutect2.mutect2 as normalMutect2_sscs {
      input:
        tumorBam = normal.sscsScBamAndIndex.bam,
        tumorBai = normal.sscsScBamAndIndex.bamIndex,
        intervalFile = intervalFile,
        intervalsToParallelizeBy = inputIntervalsToParalellizeBy,
        reference = reference,
        gatk = gatk,
        outputFileNamePrefix = normalName + "_sscs"
    }
    
    call mutect2.mutect2 as normalMutect2_allUnique {
      input:
        tumorBam = normal.allUniqueBamAndIndex.bam,
        tumorBai = normal.allUniqueBamAndIndex.bamIndex,
        intervalFile = intervalFile,
        intervalsToParallelizeBy = inputIntervalsToParalellizeBy,
        reference = reference,
        gatk = gatk,
        outputFileNamePrefix = normalName + "_allUnique"
    }

    call combineVariants as normalCombineVariants{
      input: 
        inputVcfs = [normalMutect2_dcs.filteredVcfFile, normalMutect2_sscs.filteredVcfFile],
        inputIndexes = [normalMutect2_dcs.filteredVcfFile, normalMutect2_sscs.filteredVcfIndex],
        priority = "mutect2-dcsSc,mutect2-sscsSc",
        outputPrefix = normal.outputFileNamePrefix,
        referenceFasta = resources[reference].inputRefFasta,
        modules = resources[reference].combineVariants_modules
    }

    call annotation as normalAnnotation {
      input: 
        uniqueVcf = normalMutect2_allUnique.filteredVcfFile,
        uniqueVcfIndex = normalMutect2_allUnique.filteredVcfIndex,
        mergedVcf = normalCombineVariants.combinedVcf,
        mergedVcfIndex = normalCombineVariants.combinedIndex,
        outputPrefix = normal.outputFileNamePrefix
    }

    call vep.variantEffectPredictor as normalVep {
      input: 
        vcfFile = normalAnnotation.annotatedCombinedVcf,
        vcfIndex = normalAnnotation.annotatedCombinedIndex,
        toMAF = true,
        onlyTumor = true,
        tumorOnlyAlign_updateTagValue = true,
        vcf2maf_retainInfoProvided = true,
        tumorName = tumorInputGroup.outputFileNamePrefix,
        reference = reference
    }
    File? normal_Maf = normalVep.outputMaf
  }

  if (filterMafFile && defined(tumor_Maf) && defined(normal_Maf)) {
    call filterMaf {
      input:
      mafFile = tumor_Maf,
      mafNormalFile = normal_Maf,
      outputPrefix = tumorName
    }
  }

  if (defined(normalInputGroup)) {
    InputGroup normalInput = select_first([normalInputGroup])
    
    call mutect2.mutect2 as somaticMutect2_dcs {
      input:
        tumorBam = tumorInputGroup.dcsScBamAndIndex.bam,
        tumorBai = tumorInputGroup.dcsScBamAndIndex.bamIndex,
        normalBam = normalInput.dcsScBamAndIndex.bam,
        normalBai = normalInput.dcsScBamAndIndex.bamIndex,
        intervalFile = intervalFile,
        intervalsToParallelizeBy = inputIntervalsToParalellizeBy,
        reference = reference,
        gatk = gatk,
        outputFileNamePrefix = tumorName + "_dcs_somatic"
    }
    
    call mutect2.mutect2 as somaticMutect2_sscs {
      input:
        tumorBam = tumorInputGroup.sscsScBamAndIndex.bam,
        tumorBai = tumorInputGroup.sscsScBamAndIndex.bamIndex,
        normalBam = normalInput.sscsScBamAndIndex.bam,
        normalBai = normalInput.sscsScBamAndIndex.bamIndex,
        intervalFile = intervalFile,
        intervalsToParallelizeBy = inputIntervalsToParalellizeBy,
        reference = reference,
        gatk = gatk,
        outputFileNamePrefix = tumorName + "_sscs_somatic"
    }
    
    call mutect2.mutect2 as somaticMutect2_allUnique {
      input:
        tumorBam = tumorInputGroup.allUniqueBamAndIndex.bam,
        tumorBai = tumorInputGroup.allUniqueBamAndIndex.bamIndex,
        normalBam = normalInput.allUniqueBamAndIndex.bam,
        normalBai = normalInput.allUniqueBamAndIndex.bamIndex,
        intervalFile = intervalFile,
        intervalsToParallelizeBy = inputIntervalsToParalellizeBy,
        reference = reference,
        gatk = gatk,
        outputFileNamePrefix = tumorName + "_allUnique_somatic"
    }
    
    call combineVariants as somaticCombineVariants {
      input: 
        inputVcfs = [somaticMutect2_dcs.filteredVcfFile, somaticMutect2_sscs.filteredVcfFile],
        inputIndexes = [somaticMutect2_dcs.filteredVcfIndex, somaticMutect2_sscs.filteredVcfIndex],
        priority = "mutect2-dcsSc,mutect2-sscsSc",
        outputPrefix = tumorName + "_somatic",
        referenceFasta = resources[reference].inputRefFasta,
        modules = resources[reference].combineVariants_modules
    }

    call annotation as somaticAnnotation {
      input: 
        uniqueVcf = somaticMutect2_allUnique.filteredVcfFile,
        uniqueVcfIndex = somaticMutect2_allUnique.filteredVcfIndex,
        mergedVcf = somaticCombineVariants.combinedVcf,
        mergedVcfIndex = somaticCombineVariants.combinedIndex,
        outputPrefix = tumorName + "_somatic"
    }

    call vep.variantEffectPredictor as somaticVep{
      input: 
        vcfFile = somaticAnnotation.annotatedCombinedVcf,
        vcfIndex = somaticAnnotation.annotatedCombinedIndex,
        toMAF = true,
        onlyTumor = false,
        tumorName = tumorName,
        normalName = normalName,
        tumorOnlyAlign_updateTagValue = true,
        vcf2maf_retainInfoProvided = true,
        reference = reference
    }

    if (filterMafFile) {
      call filterMaf as somaticFilterMaf {
          input:
          mafFile = somaticVep.outputMaf,
          outputPrefix = tumorName + "_somatic"
        }
    }
  }

  meta {
    author: "Alexander Fortuna, Rishi Shah and Gavin Peng"
    email: "alexander.fortuna@oicr.on.ca, rshah@oicr.on.ca, and gpeng@oicr.on.ca"
    description: "The Mutect2Consensus workflow will process umiConsensus outputs for the tumour data through mutect2 in tumour only mode to call variants then use information from the matched normal to identify likely germline variants."
    dependencies: [
     {
      name: "gatk/3.6-0",
      url: "https://gatk.broadinstitute.org"
     },
     {
      name: "python/3.9",
      url: "https://www.python.org/downloads/"
     },
     {
      name: "vep/105.0",
      url: "https://useast.ensembl.org/info/docs/tools/vep/"
     },
     {
      name: "gatk/4.1.6.0",
      url: "https://gatk.broadinstitute.org/"
     },
     {
      name: "tabix/0.2.6",
      url: "https://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2/download"
     },
     {
      name: "vcf2maf/1.6",
      url: "https://github.com/mskcc/vcf2maf"
     },
     {
      name: "pandas/1.4.2",
      url: "https://pandas.pydata.org/"
     }
    ]
    output_meta: {
      tumorVcf: {
          description: "vep vcf for tumor sample",
          vidarr_label: "tumorVcf"
      },
      tumorVcfIndex: {
          description: "vep vcf index for tumor sample",
          vidarr_label: "tumorVcfIndex"
      },
      normalVcf: {
          description: "vep vcf for normal sample",
          vidarr_label: "normalVcf"
      },
      normalVcfIndex: {
          description: "vep vcf index for normal sample",
          vidarr_label: "normalVcfIndex"
      },
      somaticVcf: {
          description: "vep vcf for somatic samples",
          vidarr_label: "somaticVcf"
      },
      somaticVcfIndex: {
          description: "vep vcf index for somatic samples",
          vidarr_label: "somaticVcfIndex"
      },
      tumorMaf: {
          description: "maf file of tumor, before filtering",
          vidarr_label: "tumorMaf"
      },
      normalMaf: {
          description: "maf file of normal, before filtering",
          vidarr_label: "normalMaf"
      },
      tumorFilteredMaf: {
          description: "tumour Maf with normal annotation and filtering",
          vidarr_label: "filterredMaf"
      },
      somaticMaf: {
          description: "Unfiltered maf file generated from mutect2 run in somatic mode, with matched tumor and normal",
          vidarr_label: "somaticMaf"
      },
      somaticFilteredMaf: {
          description: "maf file after filtering for somaticMaf",
          vidarr_label: "somaticFilterredMaf"
      }
    }
  }

  output {
    File tumorVcf = tumorVep.outputVcf
    File tumorVcfIndex = tumorVep.outputTbi
    File? normalVcf = normalVep.outputVcf
    File? normalVcfIndex = normalVep.outputTbi
    File? somaticVcf = somaticVep.outputVcf
    File? somaticVcfIndex = somaticVep.outputTbi
    File? tumorMaf = tumorVep.outputMaf
    File? normalMaf = normalVep.outputMaf
    File? somaticMaf = somaticVep.outputMaf
    File? tumorFilteredMaf = filterMaf.filteredMaf
    File? somaticFilteredMaf = somaticFilterMaf.filteredMaf
  }
}

task combineVariants {
input {
 Array[File] inputVcfs
 Array[File] inputIndexes
 Array[String] workflows
 String referenceFasta
 String outputPrefix 
 String modules
 String priority
 Int jobMemory = 24
 Int timeout = 20
 Int threads = 8
}

parameter_meta {
 inputVcfs: "array of input vcf files"
 inputIndexes: "array of tabix indexes for vcf files"
 workflows: "array of ids of producer workflows"
 referenceFasta: "path to the reference FASTA file"
 outputPrefix: "prefix for output file"
 modules: "modules for running preprocessing"
 priority: "Comma-separated list defining priority of workflows when combining variants"
 jobMemory: "memory allocated to preprocessing, in GB"
 timeout: "timeout in hours"
 threads: "number of cpu threads to be used"
}

command <<<
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
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  cpu:     "~{threads}"
  timeout: "~{timeout}"
}

output {
  File combinedVcf = "~{outputPrefix}_combined.vcf.gz"
  File combinedIndex = "~{outputPrefix}_combined.vcf.gz.tbi"
}
}

task annotation {
input {
 File uniqueVcf 
 File uniqueVcfIndex
 File mergedVcf
 File mergedVcfIndex
 String outputPrefix
 String modules = "samtools/1.9 bcftools/1.9 htslib/1.9 tabix/1.9"
 Int jobMemory = 24
 Int timeout = 20
 Int threads = 8
}

parameter_meta {
 uniqueVcf: "input unique vcf files"
 uniqueVcfIndex: "input unique tabix indexes for vcf files"
 mergedVcf: "input merged vcf"
 mergedVcfIndex: "input merged vcf index"
 outputPrefix: "prefix for output file"
 modules: "module for running preprocessing"
 jobMemory: "memory allocated to preprocessing, in GB"
 timeout: "timeout in hours"
 threads: "number of cpu threads to be used"
}

command <<<
  bcftools annotate -a ~{uniqueVcf} \
 -c FMT/AD,FMT/DP ~{mergedVcf} -Oz \
 -o "~{outputPrefix}.merged.vcf.gz"

 tabix -p vcf "~{outputPrefix}.merged.vcf.gz"
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  cpu:     "~{threads}"
  timeout: "~{timeout}"
}

output {
  File annotatedCombinedVcf = "~{outputPrefix}.merged.vcf.gz"
  File annotatedCombinedIndex = "~{outputPrefix}.merged.vcf.gz.tbi"
}
}

task filterMaf {
  input {
    File? mafFile
    File? mafNormalFile
    String freqList ="$MAF_FILTERING_ROOT/TGL.frequency.20210609.annot.txt"
    String outputPrefix 
    String modules = "python/3.9 pandas/1.4.2 maf-filtering/2024-07-10"
    Int jobMemory = 8
    Int timeout = 1
    Int threads = 1
  }

  parameter_meta {
    mafFile: "input maf file for tumor sample"
    mafNormalFile: "input file for normal sample"
    freqList: "frequency list used in maf annotation"
    outputPrefix: "prefix for output file"
    modules: "module for running preprocessing"
    jobMemory: "memory allocated to preprocessing, in GB"
    timeout: "timeout in hours"
    threads: "number of cpu threads to be used"
  }


  command <<<
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
  >>>

  runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  cpu:     "~{threads}"
  timeout: "~{timeout}"
  }

  output {
    File filteredMaf = "~{outputPrefix}_filtered_maf.gz"
  }
}
