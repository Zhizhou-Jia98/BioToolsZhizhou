function run_blast(program, query_file, subject_file, out_file, output_cols, args...; kwargs...)
    if isfile(subject_file)
        cm = "$(program) -query $(query_file) -subject $(subject_file) -out $(out_file)"
    else
        cm = "$(program) -query $(query_file) -db $(subject_file) -out $(out_file)"
    end
    if !isempty(args)
        for arg in args
            cm *= " "
            cm *= "-$(arg)"
        end
    end
    for (key, value) in kwargs
        cm *= " "
        cm *= "-$(key) $(value)"
    end
    cm *= " -outfmt \"6 $(join(output_cols, ' '))\""
    run(`bash -c $cm`)
end
