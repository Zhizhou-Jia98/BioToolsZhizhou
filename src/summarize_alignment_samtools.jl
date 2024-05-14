#=
The ways to use samtools to summarize alignments, numbers are for testing
ok    "num_total"              => "raw total sequences:" 49802 (samtools view -c -F 0x900)
ok    "num_mapped"             => "reads mapped:" 49541 (samtools view -c -F 0x904) or (num_total - num_unmapped)
ok    "num_unmapped"           => "reads unmapped:" 261 (samtools view -c -f 0x4)
ok    "num_total_pair"         => "reads paired:" / 2 = 24901 (samtools view -c -f 0x1 -F 0x900) / 2
ok    "num_unique"             => 48480 (samtools view -c -d 'NH:1' -F 0x904)
ok    "num_multi"              => 1061 ("num_mapped" - "num_unique")
ok    "num_mapped_pair"        => 24640 (samtools view -c -F 0x904 -f 2) / 2
ok    "num_unique_pair"        => 24112 (samtools view -c -d 'NH:1' -F 0x904 -f 3) / 2
ok    "num_multi_pair"         => 528 ("num_mapped_pair" - "num_unique_pair")

ok    "perc_mapped"            => 0.994759 ("num_mapped" / "num_total")
ok    "perc_uniq"              => 0.978583 ("num_multi" / "num_mapped")
ok    "perc_multi"             => 0.0214166 ("num_multi" / "num_mapped")

ok    "NC_004603.1-num_mapped" => 45154
ok    "NC_004603.1-num_unique" => 44095
ok    "NC_004603.1-num_multi"  => 1059
ok    "NC_004605.1-num_mapped" => 4387
ok    "NC_004605.1-num_unique" => 4385
ok    "NC_004605.1-num_multi"  => 2
=#

function samtools_sort_and_index(sam_file::String, num_threads::Integer=1, keep_sam::Bool=true)
    bam_file = replace(sam_file, ".sam" => "_sort.bam")
    println("Sorting $sam_file and output to $bam_file")
    @time run(pipeline(`samtools view -@ $num_threads -b $sam_file`, `samtools sort -@ $num_threads -o $bam_file`))
    println("Indexing $bam_file")
    @time run(`samtools index -@ $num_threads $bam_file`)
    if !keep_sam
        rm(sam_file)
    end
    return bam_file
end

function samtools_sort_and_index(sam_files::Vector{String}, num_threads::Integer=1, keep_sam::Bool=true)
    bam_files = String[]
    for sam_file in sam_files
        push!(bam_files, samtools_sort_and_index(sam_file, num_threads, keep_sam))
    end
    return bam_files
end

function samtools_count(type::String, bam_file::String, num_threads::Integer=1; chromosome::String="all")
    so = IOBuffer()

    if occursin("multi", type)
        if type == "multi" # "mapped" - "unique"
            println("To get the number of multimapping alignments:")
            println("\tnum_mapped - num_unique")
        elseif type == "multi_pair" # "mapped_pair" - "unique_pair"
            println("To get the number of multimapping and paired alignments:")
            println("\tnum_mapped_pair - num_unique_pair")
        end
        return nothing

    elseif occursin("pair", type)
        if type == "total_pair" # Need to divide by 2 after samtools view -c
            filters = "-f 0x3 -F 0x900"
        elseif type == "mapped_pair"
            filters = "-F 0x904 -f 0x3"
        elseif type == "unique_pair"
            filters = "-d NH:1 -F 0x904 -f 0x3"
        end

        if chromosome == "all"
            run(pipeline(`bash -c "samtools view -@ $num_threads -c $filters $bam_file"`, stdout=so))
            return Int64(parse(Int64, strip(String(take!(so)))) / 2)
        else
            run(pipeline(`bash -c "samtools view -@ $num_threads -c $filters $bam_file $chromosome"`, stdout=so))
            return Int64(parse(Int64, strip(String(take!(so)))) / 2)
        end

    else
        if type == "total"
            filters = "-F 0x900"
        elseif type == "mapped"
            filters = "-F 0x904"      
        elseif type == "unmapped"
            filters = "-f 0x4"
        elseif type == "unique"
            filters = "-d NH:1 -F 0x904"
        end
        if chromosome == "all"
            run(pipeline(`bash -c "samtools view -@ $num_threads -c $filters $bam_file"`, stdout=so))
            return parse(Int64, strip(String(take!(so))))
        else
            run(pipeline(`bash -c "samtools view -@ $num_threads -c $filters $bam_file $chromosome"`, stdout=so))
            return parse(Int64, strip(String(take!(so))))
        end
    end    
end

mutable struct AlignmentSummary
    library_name::String
    chromosomes::Vector{String}
    num_total::Dict{String, Integer}
    num_mapped::Dict{String, Integer}
    num_unmapped::Dict{String, Integer}
    num_unique::Dict{String, Integer}
    num_total_pair::Dict{String, Integer}
    num_mapped_pair::Dict{String, Integer}
    num_unique_pair::Dict{String, Integer}
    num_multi::Dict{String, Integer} # ("num_mapped" - "num_unique")
    num_multi_pair::Dict{String, Integer} # ("num_mapped_pair" - "num_unique_pair")
    perc_mapped::Dict{String, AbstractFloat}
    perc_uniq::Dict{String, AbstractFloat}
    perc_multi::Dict{String, AbstractFloat}
end

function AlignmentSummary(bam_file::String, lib_name::String, chromosomes::Vector{String}, num_threads::Int64)
    
    num_total::Dict{String, Integer} = Dict("all" => samtools_count("total", bam_file, num_threads))
    num_mapped::Dict{String, Integer} = Dict("all" => samtools_count("mapped", bam_file, num_threads))
    num_unmapped::Dict{String, Integer} = Dict("all" => samtools_count("unmapped", bam_file, num_threads))
    num_unique::Dict{String, Integer} = Dict("all" => samtools_count("unique", bam_file, num_threads))
    num_total_pair::Dict{String, Integer} = Dict("all" => samtools_count("total_pair", bam_file, num_threads))
    num_mapped_pair::Dict{String, Integer} = Dict("all" => samtools_count("mapped_pair", bam_file, num_threads))
    num_unique_pair::Dict{String, Integer} = Dict("all" => samtools_count("unique_pair", bam_file, num_threads))
    num_multi::Dict{String, Integer} = Dict("all" => num_mapped["all"] - num_unique["all"])
    num_multi_pair::Dict{String, Integer} = Dict("all" => num_mapped_pair["all"] - num_unique_pair["all"])
    perc_mapped::Dict{String, AbstractFloat} = Dict("all" => num_mapped["all"] / num_total["all"])
    perc_uniq::Dict{String, AbstractFloat} = Dict("all" => num_unique["all"] / num_mapped["all"])
    perc_multi::Dict{String, AbstractFloat} = Dict("all" => num_multi["all"] / num_mapped["all"])

    for chr in chromosomes
        num_total[chr] = samtools_count("total", bam_file, num_threads; chromosome=chr)
        num_mapped[chr] = samtools_count("mapped", bam_file, num_threads; chromosome=chr)
        num_unmapped[chr] = samtools_count("unmapped", bam_file, num_threads; chromosome=chr)
        num_unique[chr] = samtools_count("unique", bam_file, num_threads; chromosome=chr)
        num_total_pair[chr] = samtools_count("total_pair", bam_file, num_threads; chromosome=chr)
        num_mapped_pair[chr] = samtools_count("mapped_pair", bam_file, num_threads; chromosome=chr)
        num_unique_pair[chr] = samtools_count("unique_pair", bam_file, num_threads; chromosome=chr)
        num_multi[chr] = num_mapped[chr] - num_unique[chr]
        num_multi_pair[chr] = num_mapped_pair[chr] - num_unique_pair[chr]
        perc_mapped[chr] = num_mapped[chr] / num_total[chr]
        perc_uniq[chr] = num_unique[chr] / num_mapped[chr]
        perc_multi[chr] = num_multi[chr] / num_mapped[chr]
    end

    return AlignmentSummary(
        lib_name, 
        chromosomes, 
        num_total,
        num_mapped,
        num_unmapped,
        num_unique,
        num_total_pair,
        num_mapped_pair,
        num_unique_pair,
        num_multi,
        num_multi_pair,
        perc_mapped,
        perc_uniq,
        perc_multi
        )
end

function summary_table(summary_dtype::DataType, chromosomes::AbstractVector{String})
    summaries = String[]
    for (i, field) in enumerate(fieldnames(summary_dtype))
        if summary_dtype.types[i] <: Dict
            push!(summaries, [string(field) * "_" * x for x in ["all"; chromosomes]]...)
        end
    end
    return DataFrame("summaries" => summaries)
end

function summarize_multi_bams(bam_files::AbstractVector{String}, chromosomes::AbstractVector{String}, num_threads::Integer=1, keep_orginial_sam::Bool=true)
    align_sum_tbl = summary_table(AlignmentSummary, chromosomes)

    for bam_file in bam_files
        lib_name = replace(last(splitdir(bam_file)), ".bam" => "", "_sort" => "")
        if isfile(bam_file * ".bai")
            println("Index for BAM $bam_file already exists 😄")
        else
            @time "Indexing $bam_file" run(`bash -c "samtools index -@ $num_threads $bam_file"`)
        end
        @time "Summarize $lib_name" align_sum = AlignmentSummary(bam_file, lib_name, chromosomes, num_threads)
        println()
        align_sum_tbl[:, align_sum.library_name] = Vector{Number}(undef, size(align_sum_tbl)[1])
        for field in fieldnames(typeof(align_sum))
            if typeof(getfield(align_sum, field)) <: Dict
                for chr in ["all"; chromosomes]
                    i = findfirst(align_sum_tbl.:summaries .== string(field) * "_" * chr)
                    align_sum_tbl[i, align_sum.library_name] = getproperty(align_sum, field)[chr]
                end
            end
        end
    end

    return align_sum_tbl
end

# # Test
# @time run_samtools_counting("total", "alignment_RIMD2210633NCBI/DFL0041_25000_sort.bam", 2) # 49802, correct
# @time run_samtools_counting("mapped", "alignment_RIMD2210633NCBI/DFL0041_25000_sort.bam", 2) # 49541, correct
# @time run_samtools_counting("unmapped", "alignment_RIMD2210633NCBI/DFL0041_25000_sort.bam", 2) # 261, correct
# @time run_samtools_counting("unique", "alignment_RIMD2210633NCBI/DFL0041_25000_sort.bam", 2) # 48480, correct
# @time run_samtools_counting("total_pair", "alignment_RIMD2210633NCBI/DFL0041_25000_sort.bam", 2) # 24901, correct
# @time run_samtools_counting("mapped_pair", "alignment_RIMD2210633NCBI/DFL0041_25000_sort.bam", 2) # 24640, correct
# @time run_samtools_counting("unique_pair", "alignment_RIMD2210633NCBI/DFL0041_25000_sort.bam", 2) # 24112, correct
# @time run_samtools_counting("multi", "alignment_RIMD2210633NCBI/DFL0041_25000_sort.bam", 2) # correct
# @time run_samtools_counting("multi_pair", "alignment_RIMD2210633NCBI/DFL0041_25000_sort.bam", 2) # 24112, correct

# @time run_samtools_counting("mapped", "alignment_RIMD2210633NCBI/DFL0041_25000_sort.bam", 2; chromosome="NC_004603.1") # 45154, correct
# @time run_samtools_counting("unique", "alignment_RIMD2210633NCBI/DFL0041_25000_sort.bam", 2; chromosome="NC_004603.1") # 44095, correct
# @time run_samtools_counting("multi", "alignment_RIMD2210633NCBI/DFL0041_25000_sort.bam", 2; chromosome="NC_004603.1") # 45154 - 44095 = 1059, correct

# @time run_samtools_counting("mapped", "alignment_RIMD2210633NCBI/DFL0041_25000_sort.bam", 2; chromosome="NC_004605.1") # 4387, correct
# @time run_samtools_counting("unique", "alignment_RIMD2210633NCBI/DFL0041_25000_sort.bam", 2; chromosome="NC_004605.1") # 4385, correct
# @time run_samtools_counting("multi", "alignment_RIMD2210633NCBI/DFL0041_25000_sort.bam", 2; chromosome="NC_004605.1") # 4387 - 4385 = 2, correct
