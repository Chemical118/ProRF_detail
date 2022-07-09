using XLSX, DataFrames, Printf, FASTX

seq = "IEKFKLLAEKVEEIVAKNARAEIDYSDAPDEFRDPLMDTLMTDPVRLPSGTVMDRSIILRHLLNSPTDPFNRQMLTESMLEPVPELKEQIQAWMREKQSSDH"
ref = collect(seq)

excel_data = DataFrame(XLSX.readtable("Data/ube_data.xlsx", "Sheet1", infer_eltypes=true)...)[:, [:seqID, :nscor_log2_ratio]]

data_vector = Vector{Tuple{String, Float64}}()
seq_vector = Vector{String}()
for (mutstr, val) in map(Tuple, eachrow(excel_data))
    if '*' ∉ mutstr && !issubset("NA", mutstr)
        ref_seq = deepcopy(ref)
        mutstr_vector = Vector{String}()
        for (snum, mut) in zip(map(x -> split(x, ','), split(mutstr, '-'))...)
            num = parse(Int, snum) + 1
            push!(mutstr_vector, ref_seq[num] * string(num) * mut)
            ref_seq[num] = only(mut)
        end

        push!(data_vector, (join(mutstr_vector, ':'), val))
        push!(seq_vector, join(ref_seq))
    end
end

@printf "%d %d\n" length(seq_vector) length(data_vector)

open(FASTA.Writer, "Data/eUBedata.fasta") do io
    for (seq, id) in zip(seq_vector,  map(x -> x[1], data_vector))
        write(io, FASTA.Record(id, seq))
    end
end

XLSX.openxlsx("Data/eUBedata.xlsx", mode="w") do xf
    sheet = xf[1]
    XLSX.rename!(sheet, "Sheet1")
    sheet["A1"] = "aaMutations"
    sheet["B1"] = "Enrichment_score"
    sheet["A2", dim=1] = map(x -> x[1], data_vector)
    sheet["B2", dim=1] = map(x -> x[2], data_vector);
end