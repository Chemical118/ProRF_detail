using XLSX, DataFrames, Printf, FASTX

seq = "MQYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE"
ref = collect(seq)

excel_data = DataFrame(XLSX.readtable("Data/gb1_raw.xlsx", "double_mutants", infer_eltypes=true)...)

WinputC, WselectC = 1759616, 3041819

data_vector = Vector{Tuple{String, Float64}}()
seq_vector = Vector{String}()
for (w1, l1, m1, w2, l2, m2, inputC, selectC) in map(Tuple, eachrow(excel_data))
    if m1[end] != '*' && m2[end] != '*'
        ref_seq = deepcopy(ref)
        ref_seq[l1] = m1[end]
        ref_seq[l2] = m2[end]
        push!(seq_vector, join(ref_seq))
        push!(data_vector, (String(@sprintf "%c%d%c:%c%d%c" ref[l1] l1 m1[end] ref[l2] l2 m2[end]), log2((WinputC + 1) / (WselectC + 1)) - log2((inputC + 1) / (selectC + 1))))
    end
end

@printf "%d %d" length(seq_vector) length(data_vector)

open(FASTA.Writer, "Data/eGBdata.fasta") do io
    for (seq, id) in zip(seq_vector,  map(x -> x[1], data_vector))
        write(io, FASTA.Record(id, seq))
    end
end

XLSX.openxlsx("Data/eGBdata.xlsx", mode="w") do xf
    sheet = xf[1]
    XLSX.rename!(sheet, "Sheet1")
    sheet["A1"] = "aaMutations"
    sheet["B1"] = "Enrichment_score"
    sheet["A2", dim=1] = map(x -> x[1], data_vector)
    sheet["B2", dim=1] = map(x -> x[2], data_vector)
end