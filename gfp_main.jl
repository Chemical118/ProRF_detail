using XLSX, DataFrames, Printf, FASTX, BioSequences

seq = [FASTA.sequence(record) for record in open(FASTA.Reader, "Data/avref.fasta")][1]
ref = collect("M" * string(translate(seq)[1:end-1])) # end codon

excel_data = DataFrame(XLSX.readtable("Data/Gdata.xlsx", "Sheet1", infer_eltypes=true)...)[:, [:aaMutations, :medianBrightness]]

data_vector = Vector{Tuple{String, Float64}}()
seq_vector = Vector{String}()
XLSX.openxlsx("Data/eGdata.xlsx", mode="w") do xf

sheet = xf[1]
XLSX.rename!(sheet, "Sheet1")
sheet["A1"] = "aaMutations"
sheet["B1"] = "medianBrightness"

for (mutstr, val) in map(Tuple, eachrow(excel_data))
    if '*' ∉ mutstr
        ref_seq = deepcopy(ref)
        mutstr_vector = Vector{String}()

        for mutar in split(mutstr, ":")
            mut = mutar[end]
            num = parse(Int, mutar[3:end-1]) + 2
            if only(mutar[2]) != ref[num]
                error("Data!")
            end
            ref_seq[num] = mut
            push!(mutstr_vector, mutar[2] * string(num) * mutar[end])
        end

        push!(data_vector, (join(mutstr_vector, ':'), val))
        push!(seq_vector, join(ref_seq))
    end
end

sheet["A2", dim=1] = map(x -> x[1], data_vector)
sheet["B2", dim=1] = map(x -> x[2], data_vector)
end # xlsx end

open(FASTA.Writer, "Data/eGdata.fasta") do io
    for (seq, id) in zip(seq_vector,  map(x -> x[1], data_vector))
        write(io, FASTA.Record(id, seq))
    end
end