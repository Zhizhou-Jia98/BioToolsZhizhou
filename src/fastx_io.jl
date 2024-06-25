function write_multiple_fa(fa_descriptions::AbstractVector,
                           starts::AbstractVector,
                           ends::AbstractVector,
                           strands::AbstractVector,
                           chromosomes::AbstractVector,
                           genome_fa::AbstractString,
                           out_file::AbstractString;
                           full_description::Bool=false)
    ispath("$genome_fa.fai") ? nothing : faidx(genome_fa)
    FASTX.FASTA.Reader(open(genome_fa); index = "$genome_fa.fai") do f_genome
        FASTX.FASTA.Writer(open(out_file, "w")) do f_write
            for i in eachindex(fa_descriptions)
                if full_description
                    seq_id = string(fa_descriptions[i], " ", chromosomes[i], ":", starts[i], "-", ends[i], ":", strands[i])
                else
                    seq_id = fa_descriptions[i]
                end

                if strands[i] == "+"
                    seq = sequence(f_genome[chromosomes[i]])[starts[i]:ends[i]]
                elseif strands[i] == "-"
                    seq = reverse_complement(sequence(LongDNA{4}, f_genome[chromosomes[i]])[starts[i]:ends[i]])
                end

                write(f_write, FASTA.Record(seq_id, seq))
            end
        end
    end
end
