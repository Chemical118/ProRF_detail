using XLSX, DataFrames, Printf, FASTX

# https://www.ncbi.nlm.nih.gov/protein/417441
seq = [FASTA.sequence(record) for record in open(FASTA.Reader, "Data/pab1ref.fasta")][1]
ref = map(x -> only(string(x)), collect(seq))

excel_data = DataFrame(XLSX.readtable("Data/pab1data.xlsx", "All_Epistasis", infer_eltypes=true))

data_vector = Vector{Tuple{String, Float64}}()
seq_vector = Vector{String}()
for (m12, v12, m1, v1, m2, v2) in map(Tuple, eachrow(excel_data))
    if m1[end] != '*'
        ml1 = parse(Int, split(m1, '-')[1])
    
        ref_seq = deepcopy(ref)
        ref_seq[ml1] = m1[end]
        push!(seq_vector, join(ref_seq))
        push!(data_vector, (String(@sprintf "%c%d%c" ref[ml1] ml1 m1[end]), v1))
    
        if m2[end] != '*'
            ml2 = parse(Int, split(m2, '-')[1])
    
            ref_seq = deepcopy(ref)
            ref_seq[ml2] = m2[end]
            push!(seq_vector, join(ref_seq))
            push!(data_vector, (String(@sprintf "%c%d%c" ref[ml2] ml2 m2[end]), v2))
    
            ref_seq = deepcopy(ref)
            ref_seq[ml1] = m1[end]
            ref_seq[ml2] = m2[end]
            push!(seq_vector, join(ref_seq))
            push!(data_vector, (String(@sprintf "%c%d%c:%c%d%c" ref[ml1] ml1 m1[end] ref[ml2] ml2 m2[end]), v12))
        end
    elseif m2[end] != '*'
        ml2 = parse(Int, split(m2, '-')[1])
    
        ref_seq = deepcopy(ref)
        ref_seq[ml2] = m2[end]
        push!(seq_vector, join(ref_seq))
        push!(data_vector, (String(@sprintf "%c%d%c" ref[ml2] ml2 m2[end]), v2))
    end
end

@printf "%d %d" length(seq_vector) length(data_vector)

a = length(filter(x -> '*' ∉ x, eachcol(excel_data)[1]))
b = length(filter(x -> '*' ∉ x, eachcol(excel_data)[3]))
c = length(filter(x -> '*' ∉ x, eachcol(excel_data)[5]))
println(a + b + c)

open(FASTA.Writer, "Data/ePdata.fasta") do io
    for (seq, id) in zip(seq_vector,  map(x -> x[1], data_vector))
        write(io, FASTA.Record(id, seq))
    end
end

XLSX.openxlsx("Data/ePdata.xlsx", mode="w") do xf
    sheet = xf[1]
    XLSX.rename!(sheet, "Sheet1")
    sheet["A1"] = "aaMutations"
    sheet["B1"] = "Enrichment score"
    sheet["C1"] = "Log Enrichment score"
    sheet["A2", dim=1] = map(x -> x[1], data_vector)
    sheet["B2", dim=1] = map(x -> x[2], data_vector)
    sheet["C2", dim=1] = log.(map(x -> x[2], data_vector))
end