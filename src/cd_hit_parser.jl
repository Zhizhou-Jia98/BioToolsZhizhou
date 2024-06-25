# using DataFrames


function parse_cd_hit_to_table(cd_hit_file)
    cdhit_table = DataFrame(
        :cluster_id => Int64[],
        :is_representative=> Bool[],
        :seq_id=> Int64[],
        :seq_length=> Int64[],
        :seq_name=> AbstractString[],
        :seq_chrom=> AbstractString[],
        :seq_description=> AbstractString[],
        :align_start_self=> Union{Missing, Int64}[],
        :align_end_self=> Union{Missing, Int64}[],
        :align_start_representative=> Union{Missing, Int64}[],
        :align_end_representative=> Union{Missing, Int64}[],
        :align_strand=> Union{Missing, AbstractString}[],
        :align_identity=> Union{Missing, Float64}[],
    )
    open(cd_hit_file) do f
        cluster_id = 0
        while !eof(f)
            line = readline(f)
            if startswith(line, '>')
                cluster_id = parse(Int64, last(split(line)))
            else
                # Example of a non-representative line
                # ["0", "430nt,", ">DFS0007_02733:DFS0007_chr1:2873081-2873510:-...", "at", "1:430:1:430/+/99.77%"]
                is_representative = endswith(line, "*") ? true : false
                parsed_line = split(line)
                seq_id = parse(Int64, parsed_line[1])
                seq_length = parse(Int64, chop(parsed_line[2], tail=3))
                seq_name = chop(first(split(parsed_line[3], ":")), head=1)
                seq_chrom = split(parsed_line[3], ":")[2]
                seq_description = chop(parsed_line[3], head=1, tail=3)
                align_start_self = is_representative ? missing : parse(Int64, split(first(split(parsed_line[5], "/")), ":")[1])
                align_end_self = is_representative ? missing : parse(Int64, split(first(split(parsed_line[5], "/")), ":")[2])
                align_start_representative = is_representative ? missing : parse(Int64, split(first(split(parsed_line[5], "/")), ":")[3])
                align_end_representative = is_representative ? missing : parse(Int64, split(first(split(parsed_line[5], "/")), ":")[4])
                align_strand = is_representative ? missing : split(parsed_line[5], "/")[2]
                align_identity = is_representative ? missing : parse(Float64, chop(split(parsed_line[5], "/")[end], tail=1))
                push!(cdhit_table, [cluster_id, is_representative, seq_id, seq_length, seq_name, seq_chrom, seq_description, align_start_self, align_end_self, align_start_representative, align_end_representative, align_strand, align_identity])
            end
        end
    end

    cdhit_table.:cluster_size = ones(Int64, size(cdhit_table)[1])
    println(names(cdhit_table))
    for cluster_id in unique(cdhit_table.:cluster_id)
        cluster_size = sum(cdhit_table.:cluster_id .== cluster_id)
        cdhit_table[findall(cdhit_table.:cluster_id .== cluster_id), :cluster_size] .= cluster_size
    end
    return cdhit_table
end

# Test
# @time cd_hit_table = parse_cdhit_to_table("/Users/zhizhoujia/Desktop/srna_cdhit/est_ncrna_c85_s85.clstr")
# CSV.write("/Users/zhizhoujia/Desktop/test.csv", cdhit_tbl)

struct CdHitSequence
    id_in_cluster::Int64
    length::Int64
    name::AbstractString
    chromosome::AbstractString
    description::AbstractString
    alignment::AbstractString
    is_representative::Bool
end

struct CdHitCluster
    id::Int64
    size::Int64
    sequences::Vector{CdHitSequence}
    # CdHitCluster(id, size, sequences) = new(id, size, deepcopy(sequences)) # Slower
end

function parse_cd_hit(cd_hit_file)
    cd_hit_clusters = CdHitCluster[]
    open(cd_hit_file) do f
        cluster_id = -1
        cluster_size = 0
        cd_hit_sequences = CdHitSequence[]
        while !eof(f)
            line = readline(f)
            if startswith(line, '>')
                cluster_size = length(cd_hit_sequences)
                cd_hit_cluster = CdHitCluster(cluster_id, cluster_size, cd_hit_sequences)
                push!(cd_hit_clusters, cd_hit_cluster)
                cd_hit_sequences = CdHitSequence[] # Either create a new vector, or deepcopy and empty the vector (slower)!!!!! ALWAYS remember!!!!!
                # empty!(cd_hit_sequences)
                cluster_id = parse(Int64, last(split(line)))
            else
                # Example of a non-representative line
                # ["0", "430nt,", ">DFS0007_02733:DFS0007_chr1:2873081-2873510:-...", "at", "1:430:1:430/+/99.77%"]
                parsed_line = split(line)
                seq_id = parse(Int64, parsed_line[1])
                seq_length = parse(Int64, chop(parsed_line[2], tail=3))
                seq_name = chop(split(parsed_line[3], ":")[1], head=1, tail=0)
                seq_chrom = split(parsed_line[3], ":")[2]
                seq_description = chop(parsed_line[3], head=1, tail=3)
                align_info = last(parsed_line)
                is_representative = align_info == "*" ? true : false
                cd_hit_sequence = CdHitSequence(seq_id, seq_length, seq_name, seq_chrom, seq_description, align_info, is_representative)
                push!(cd_hit_sequences, cd_hit_sequence)
            end
        end
        push!(cd_hit_clusters, CdHitCluster(cluster_id, cluster_size, cd_hit_sequences))
    end
    return cd_hit_clusters[2:end]
end

# Test
# @time cd_hit_clusters = parse_cd_hit("/Users/zhizhoujia/Research/IPS/daniel_group/Vibrio_parahaemolyticus/Vp_reference_sequences_and_annotation/pangenome/vp_ecogroup/trycycler_assemblies/srna_cdhit/est_ncrna_c85_s85.clstr");
