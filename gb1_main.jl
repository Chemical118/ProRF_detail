using XLSX, DataFrames, Printf, FASTX

seq = "MQYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE"
ref = collect(seq);

WinputC, WselectC = 1759616, 3041819

excel_data = DataFrame(XLSX.readtable("Data/gb1_raw.xlsx", "single_mutants", infer_eltypes=true)...)

data_vector = Vector{Tuple{String, Float64}}()
seq_vector = Vector{String}()
for (w1, l1, m1, inputC, selectC) in map(Tuple, eachrow(excel_data))
    if only(m1) != '*'
        if only(w1) != ref[l1]
            error("Data?")
        end
        ref_seq = deepcopy(ref)
        ref_seq[l1] = only(m1)
        push!(seq_vector, join(ref_seq))
        push!(data_vector, (String(@sprintf "%c%d%c" ref[l1] l1 only(m1)), (selectC / inputC) / (WselectC / WinputC)))
    end
end

excel_data = DataFrame(XLSX.readtable("Data/gb1_raw.xlsx", "double_mutants", infer_eltypes=true)...)

for (w1, l1, m1, w2, l2, m2, inputC, selectC) in map(Tuple, eachrow(excel_data))
    if only(m1) != '*' && only(m2) != '*'
        if only(w1) != ref[l1] || only(w2) != ref[l2]
            error("Data?")
        end
        ref_seq = deepcopy(ref)
        ref_seq[l1] = only(m1)
        ref_seq[l2] = only(m2)
        push!(seq_vector, join(ref_seq))
        push!(data_vector, (String(@sprintf "%c%d%c:%c%d%c" ref[l1] l1 only(m1) ref[l2] l2 only(m2)), (selectC / inputC) / (WselectC / WinputC)))
    end
end

@printf "%d %d\n" length(seq_vector) length(data_vector)

fit_vector = map(x -> x[2], data_vector)
min_val = minimum(filter(x -> x > 0, fit_vector))
log_fit_vector = log.(fit_vector .+ min_val);

open(FASTA.Writer, "Data/eGBdata.fasta") do io
    for (seq, id) in zip(seq_vector,  map(x -> x[1], data_vector))
        write(io, FASTA.Record(id, seq))
    end
end

XLSX.openxlsx("Data/eGBdata.xlsx", mode="w") do xf
    sheet = xf[1]
    XLSX.rename!(sheet, "Sheet1")
    sheet["A1"] = "aaMutations"
    sheet["B1"] = "Fitness"
    sheet["C1"] = "log Fitness"
    sheet["A2", dim=1] = map(x -> x[1], data_vector)
    sheet["B2", dim=1] = fit_vector
    sheet["C2", dim=1] = log_fit_vector
end