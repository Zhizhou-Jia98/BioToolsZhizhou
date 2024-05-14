function write_fa_from_table(id, start, stop, strand, chromosome, genome_file, out_file::AbstractString; full_id::Bool=false)
    FASTX.FASTA.Reader(open(genome_file); index = "$genome_file.fai") do f_genome
        FASTX.FASTA.Writer(open(out_file, "w")) do f_write
            for i in eachindex(id)
                if full_id
                    seq_id = id[i] * ":" * chromosome[i] * ":" * start[i] * "-" * stop[i]
                else
                    seq_id = id[i]
                end

                if strand[i] == "+"
                    seq = sequence(f_genome[chromosome[i]])[start[i]:stop[i]]
                elseif strand[i] == "-"
                    seq = reverse_complement(sequence(LongDNA{4}, f_genome[chromosome[i]])[start[i]:stop[i]])
                end

                write(f_write, FASTA.Record(seq_id, seq))
            end
        end
    end
end