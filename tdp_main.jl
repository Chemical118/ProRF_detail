using XLSX, DataFrames, Printf, FASTX

seq = "GNSRGGGAGLGNNQGSNMGGGMNFGAFSINPAMMAAAQAALQSSWGMMGMLASQQNQSGPSGNNQNQGNMQREPNQAFGSGNNS"
ref = collect(seq)

excel_loc = "Data/tdp_data.xlsx"
data_vector = Vector{Tuple{String, Float64, Float64}}()
seq_vector = Vector{String}()

# 1 mutation
excel_data = DataFrame(XLSX.readtable(excel_loc, "1mut", infer_eltypes=true)...)
for (ml1, wt1, m1, s, v, mstr1) in map(Tuple, eachrow(excel_data))
    ml1 = parse(Int, mstr1[2:end-1]) - 289
    if ref[ml1] != only(wt1)
        error("Data!")
    end
    ref_seq = deepcopy(ref)
    ref_seq[ml1] = only(m1)
    push!(seq_vector, join(ref_seq))
    push!(data_vector, (mstr1, v, s))
end

# 2 mutation
excel_data = DataFrame(XLSX.readtable(excel_loc, "2mut", infer_eltypes=true)...)
for (ml1, ml2, wt1, wt2, m1, m2, s, v, mstr1, mstr2) in map(Tuple, eachrow(excel_data))
    ml1 = parse(Int, mstr1[2:end-1]) - 289
    ml2 = parse(Int, mstr2[2:end-1]) - 289
    if ref[ml1] != only(wt1) || ref[ml2] != only(wt2)
        error("Data!")
    end
    ref_seq = deepcopy(ref)
    ref_seq[ml1] = only(m1)
    ref_seq[ml2] = only(m2)
    push!(seq_vector, join(ref_seq))
    push!(data_vector, (mstr1 * ":" * mstr2, v, s))
end

@printf "%d %d\n" length(seq_vector) length(data_vector)

open(FASTA.Writer, "Data/eTDPdata.fasta") do io
    for (seq, id) in zip(seq_vector,  map(x -> x[1], data_vector))
        write(io, FASTA.Record(id, seq))
    end
end

XLSX.openxlsx("Data/eTDPdata.xlsx", mode="w") do xf
    sheet = xf[1]
    XLSX.rename!(sheet, "Sheet1")
    sheet["A1"] = "aaMutations"
    sheet["B1"] = "toxicity"
    sheet["C1"] = "dtoxicity"
    sheet["A2", dim=1] = map(x -> x[1], data_vector)
    sheet["B2", dim=1] = map(x -> x[2], data_vector)
    sheet["C2", dim=1] = map(x -> x[3], data_vector)
end