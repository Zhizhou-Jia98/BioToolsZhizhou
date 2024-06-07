# using CSV, DataFrames


function wig_to_df(wig_file::AbstractString)
    wig_df = DataFrame(:chrom => AbstractString[], :pos => Integer[], :value => Float64[])
    open(wig_file) do f
        chrom = ""
        first_line = readline(f)
        if !occursin("track type=wiggle_0", first_line)
            @info "Maybe the file is wrongly formatted? Please check 😈"
        end
        while !eof(f)
            line = readline(f)
            row_content = split(line, " ")
            if !isnothing(tryparse(Int64, row_content[1]))
                push!(wig_df, [chrom, parse(Int64, row_content[1]), parse(Float64, row_content[2])])
            elseif row_content[1] == "variableStep"
                # @info occursin.("chrom=", row_content)
                # @info row_content[findfirst(occursin.("chrom=", row_content))]
                chrom = last(split(row_content[findfirst(occursin.("chrom=", row_content))], "="))
            elseif row_content[1] == "fixedStep"
                @warn "Don't like fixedStep!!!😉"
            end
        end
    end
    return wig_df
end
