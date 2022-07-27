using XLSX, DataFrames

excel_data = DataFrame(XLSX.readtable(raw"Data\GB1p\data.xlsx", "Sheet1", infer_eltypes=true))

XLSX.openxlsx(raw"Data\GB1p\data.xlsx", mode="rw") do xf
    sheet = xf[1]
    min_val = minimum(filter(x -> x > 0, excel_data[!, 5]))
    sheet["F2", dim=1] = log.(excel_data[!, 5] .+ min_val)
end
