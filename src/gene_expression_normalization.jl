
function featurecounts_normalization(fc_tbl_file::AbstractString, output_file::AbstractString)
	fc_tbl = CSV.read(fc_tbl_file, DataFrame; delim = "\t", comment = "#")

	lib_names = names(fc_tbl)[7:end]
	lib_names = replace.(last.(splitdir.(lib_names)), ".bam" => "_raw")
	rename!(fc_tbl, vcat(names(fc_tbl)[1:6], lib_names))

	for lib_name in lib_names
		lib_name_tpm = replace(lib_name, "_raw" => "_tpm")
		fc_tbl[!, lib_name_tpm] = calculate_tpm(fc_tbl.:Length, fc_tbl[!, lib_name])
	end
	for lib_name in lib_names
		lib_name_rpkm = replace(lib_name, "_raw" => "_rpkm")
		fc_tbl[!, lib_name_rpkm] = calculate_rpkm(fc_tbl.:Length, fc_tbl[!, lib_name])
	end

	CSV.write(output_file, fc_tbl; delim="\t")
	return fc_tbl
end


function calculate_rpkm(gene_length::AbstractVector, expression::AbstractVector)
	permil = sum(expression) / 1_000_000
	rpm = expression ./ permil
	return rpm ./ (gene_length ./ 1_000)
end

function calculate_tpm(gene_length::AbstractVector, expression::AbstractVector)
	rpk = expression ./ (gene_length ./ 1_000)
	permil = sum(rpk) / 1_000_000
	return rpk ./ permil
end
