# using CSV, DataFrames, XLSX

# gff file has 3 parts
# Part 1: directives (starts with ##) & comments (starts with #)
# Part 2: annotation -> can be parsed into a DataFrame
# Part 3: sequences -> can be parsed into a Dict

## Read Part 1 of gff file
mutable struct GffHeader{T<:Vector{String}}
    directives::T
    comments::T
end

function get_gff_header(gff_file::String)
    gff_header = GffHeader(String[], String[])
    open(gff_file) do f
        while !eof(f)
            line = strip(readline(f))
            # println(line)
            if occursin(r"^##", line)
                push!(gff_header.directives, line)
            elseif occursin(r"^#[a-zA-Z0-9]", line)
                push!(gff_header.comments, line)
            end
        end
    end
    return gff_header
end

function parse_gff(gff_file::String)
    gff = CSV.read(gff_file, DataFrame; 
                   comment="#", delim="\t", header=false,
                   pool=false, stringtype=String)
    if any(startswith.(gff.:Column1, ">"))
        gff = gff[1:first(findall(startswith.(gff.:Column1, ">"))) - 1, :]
    end
    rename!(gff, ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])
    for i in 1:size(gff)[1]
        attributes = Dict(split.(split(gff.:attributes[i], ";"), "="))
        for (key, value) in attributes
            if key in names(gff)
                gff[i, key] = value
            else
                gff[!, key] = Vector{Union{Missing, String}}(missing, size(gff)[1])
                gff[i, key] = value
            end
        end
    end
    return gff
end


function parse_gff(gff_raw::AbstractDataFrame)
    gff = copy(gff_raw)
    rename!(gff, ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])
    for i in 1:size(gff)[1]
        attributes = Dict(split.(split(gff.:attributes[i], ";"), "="))
        for (key, value) in attributes
            if key in names(gff)
                gff[i, key] = value
            else
                gff[!, key] = Vector{Union{Missing, String}}(missing, size(gff)[1])
                gff[i, key] = value
            end
        end
    end
    return gff
end

function gff_meta_info(gff_file::String)
    metainfo = open(gff_file) do f
        comment = Vector{String}()
        while !eof(f)
            line = readline(f)
            if startswith(line, "#")
                push!(comment, line)
            else
                break
            end
        end
        comment
    end
    return metainfo
end

function merge_attributes(gff::AbstractDataFrame)
    attributes = gff[:, 10:end]

    for i in axes(attributes, 1), j in axes(attributes, 2)
        if !ismissing(attributes[i, j])
            # Check if ";" is present in attributes, if so, replace it with "%3B"
            # Check and replace reserved characters in "attributes" column:
            #   ";" => "%3B"
            #   "=" => "%3D"
            #   "&" => "%26"
            #   "," => "%2C"
            if occursin(r"[;,=&]", attributes[i, j])
                # @info "Replacing character \";\" with \"%3B\" in row $i of attribute \"$(names(attributes)[j])\""
                attributes[i, j] = replace(attributes[i, j], ";" => "%3B ", "=" => "%3D", "&" => " %26 ", "," => "%2C ")
            end
            attributes[i, j] = string(names(attributes)[j], "=", attributes[i,j])
        end
    end
    return select(attributes, AsTable(:) => ByRow(x -> join(collect(skipmissing(x)), ";")))[:, 1]
end

function gff_fasta(gff_file::String)
    open(gff_file) do f
        lines = readlines(f)
        i = first(findall(occursin.("##FASTA", lines)))
        fa = join(lines[i+1:end], "\n")
        rdr = FASTAReader(IOBuffer(fa))
        records = collect(rdr)
        close(rdr)
        return records
    end
end