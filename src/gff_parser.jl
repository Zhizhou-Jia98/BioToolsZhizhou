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

function get_gff_meta_info(gff_file::String)
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

function parse_gff_attribute(attr::AbstractString)
    attr_split = split.(strip.(split(attr, ";")), "=")
    return Dict(attr_split)
end

function parse_gtf_attribute(attr::AbstractString)
    attr_temp = rstrip.(rstrip.(attr), ';')
    attr_split = split.(strip.(split(attr_temp, "\"; ")), " ")
    for i in eachindex(attr_split)
        if length(attr_split[i]) > 2
            name = first(attr_split[i])
            content = join(attr_split[i][2:end], " ")
            attr_split[i] = [name, content]
        end
        attr_split[i][2] = strip(attr_split[i][2], '\"')
    end
    return Dict(attr_split)
end

function parse_gff(gff_file::String)
    gff = CSV.read(gff_file, DataFrame;
                   comment="#", delim="\t", header=false,
                   pool=false, stringtype=String)
    filetype = last(split(gff_file, "."))
    if any(startswith.(gff.:Column1, ">"))
        gff = gff[1:first(findall(startswith.(gff.:Column1, ">"))) - 1, :]
    end
    rename!(gff, ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])
    # if endswith(gff_file, "gff")
    #     attr_delim = "\t"
    # elseif endswith(gff_file, "gtf")
    #     gff.:attributes = rstrip.(rstrip.(gff.:attributes), ';')
    #     attr_delim = " "
    # end
    for i in 1:size(gff)[1]
        # println(split.(strip.(split(gff.:attributes[i], ";")), attr_delim))
        # attributes = Dict()
        # attributes = try
        #     attributes_split = split.(strip.(split(gff.:attributes[i], ";")), attr_delim)
        #     Dict(split.(strip.(split(gff.:attributes[i], ";")), attr_delim))
        # catch
        #     println(split.(strip.(split(gff.:attributes[i], ";")), attr_delim))
        # end
        if filetype == "gff"
            attributes = parse_gff_attribute(gff.attributes[i])
        elseif filetype == "gtf"
            attributes = parse_gtf_attribute(gff.attributes[i])
        end
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


function merge_attributes(gff::AbstractDataFrame)
    if names(gff)[1:9] == ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
        attributes = gff[:, 10:end]
    else
        attributes = gff[:, :]
    end
    for j in axes(attributes, 2)
        attributes[!, j] = map(x -> ismissing(x) ? missing : string(x), attributes[!, j])
    end
    for i in axes(attributes, 1), j in axes(attributes, 2)
        if !ismissing(attributes[i, j])
            # Check if ";" is present in attributes, if so, replace it with "%3B"
            # Check and replace reserved characters in "attributes" column:
            #   ";" => "%3B"
            #   "=" => "%3D"
            #   "&" => "%26"
            #   "," => "%2C"
            #   "\n" => ""
            if isa(attributes[i, j], AbstractString) && occursin(r"[;,=&]", attributes[i, j])
                # @info "Replacing character \";\" with \"%3B\" in row $i of attribute \"$(names(attributes)[j])\""
                attributes[i, j] = replace(attributes[i, j], ";" => "%3B ", "=" => "%3D", "&" => " %26 ", "," => "%2C ", "\n" => "")
                attributes[i, j] = replace(attributes[i, j], ";" => "%3B ", "=" => "%3D", "&" => " %26 ", "," => "%2C ", "\n" => "")
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

# function write_gff()
#     
# end

