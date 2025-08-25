nextflow.enable.dsl=2

// PARAMETERS
params.expression_h5ad = null
params.prefix = "pidc_output"

// PROCESS
process PIDC {
  tag "$prefix"
  publishDir "results", mode: 'copy'

  input:
  path expression_h5ad
  val prefix
  path convert_script
  path run_script
  path pidc_julia_script

  output:
  path "${prefix}_rankedEdges.csv"

  script:
  """
  echo "Using Python from: \$(which python)"
  python3 -c "import scanpy; print('scanpy is available')"

  python3 ${convert_script} ${expression_h5ad} expression.tsv
  bash ${run_script} expression.tsv ${prefix}_rankedEdges.csv
  """
}

// WORKFLOW
workflow {
  Channel
    .fromPath(params.expression_h5ad)
    .ifEmpty { error "File not found: ${params.expression_h5ad}" }
    .set { expression_h5ad_ch }

  convert_script_ch = Channel.value(file("scripts/h5ad_to_tsv.py"))
  run_script_ch     = Channel.value(file("modules/grn/pidc/run_pidc.sh"))
  julia_script_ch   = Channel.value(file("modules/grn/pidc/runPIDC.jl")) 

  PIDC(
    expression_h5ad_ch,
    params.prefix,
    convert_script_ch,
    run_script_ch,
    julia_script_ch
  )
}
