using XLSX, DataFrames, Printf, FASTX

excel_data = DataFrame(XLSX.readtable("Data/gGdata.xlsx", "Sheet1", infer_eltypes=true)...)

data_vector = Vector{Tuple{String, Float64}}()
seq_vector = Vector{String}()
for (seq, mutstr, val) in map(Tuple, eachrow(excel_data))
    push!(data_vector, (mutstr, val))
    push!(seq_vector, seq)
end

@printf "%d %d\n" length(seq_vector) length(data_vector)

open(FASTA.Writer, "Data/egGBedata.fasta") do io
    for (seq, id) in zip(seq_vector,  map(x -> x[1], data_vector))
        write(io, FASTA.Record(id, seq))
    end
end

XLSX.openxlsx("Data/egGBedata.xlsx", mode="w") do xf
    sheet = xf[1]
    XLSX.rename!(sheet, "Sheet1")
    sheet["A1"] = "aaMutations"
    sheet["B1"] = "ddG_mean (kcal/mol)"
    sheet["A2", dim=1] = map(x -> x[1], data_vector)
    sheet["B2", dim=1] = map(x -> x[2], data_vector)
end
