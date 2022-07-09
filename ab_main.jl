using XLSX, DataFrames, Printf, FASTX

seq = "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA"
ref = collect(seq)

data_vector = Vector{Tuple{String, Float64, Float64}}()
seq_vector = Vector{String}()

# 1 mutation
excel_data = DataFrame(XLSX.readtable("Data/amb_data.xlsx", "1mut", infer_eltypes=true)...)
for (ml1, wt1, m1, v, s) in map(Tuple, eachrow(excel_data))
    if ref[ml1] != only(wt1)
        error("Data Error!")
    end
    ref_seq = deepcopy(ref)
    ref_seq[ml1] = only(m1)
    push!(seq_vector, join(ref_seq))
    push!(data_vector, (String(@sprintf "%c%d%c" ref[ml1] ml1 m1), v, s))
end

excel_data = DataFrame(XLSX.readtable("Data/amb_data.xlsx", "2mut", infer_eltypes=true)...)
for (ml2, m2, ml1, m1, wt1, wt2, v, s) in map(Tuple, eachrow(excel_data))
    if ref[ml1] != only(wt1) || ref[ml2] != only(wt2)
        error("Data Error!")
    end
    ref_seq = deepcopy(ref)
    ref_seq[ml1] = only(m1)
    ref_seq[ml2] = only(m2)
    push!(seq_vector, join(ref_seq))
    push!(data_vector, (String(@sprintf "%c%d%c:%c%d%c" ref[ml1] ml1 m1 ref[ml2] ml2 m2), v, s))
end

@printf "%d %d\n" length(seq_vector) length(data_vector)

open(FASTA.Writer, "Data/eABdata.fasta") do io
    for (seq, id) in zip(seq_vector,  map(x -> x[1], data_vector))
        write(io, FASTA.Record(id, seq))
    end
end

XLSX.openxlsx("Data/eABdata.xlsx", mode="w") do xf
    sheet = xf[1]
    XLSX.rename!(sheet, "Sheet1")
    sheet["A1"] = "aaMutations"
    sheet["B1"] = "score"
    sheet["C1"] = "dscore"
    sheet["A2", dim=1] = map(x -> x[1], data_vector)
    sheet["B2", dim=1] = map(x -> x[2], data_vector)
    sheet["C2", dim=1] = map(x -> x[3], data_vector)
end