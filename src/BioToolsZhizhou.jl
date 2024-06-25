module BioToolsZhizhou

# Write your package code here.
using CSV, DataFrames
using FASTX, BioSequences

# export samtools_sort_and_index, summarize_multi_bams, parse_gff, write_fa_from_table

include("fastx_io.jl")
include("gene_expression_normalization.jl")
include("gff_parser.jl")
include("cd_hit_parser.jl")
include("summarize_alignment_samtools.jl")
include("wig_parser.jl")

end
