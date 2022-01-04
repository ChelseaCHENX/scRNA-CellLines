using NetworkInference
cd("/bgfs/alee/chelsea/projects/10X/CellLine/codes")

CellTypes = ["MCF7", "T47D WT", "T47D KO", "MM134", "SUM44", "BCK4",  "MCF10A", "HEK293"]
for cell in CellTypes
    infer_network("../data/Counts//1000g.$cell.logNormCts.tsv", PIDCNetworkInference(), out_file_path="../data/Network/1000g.$cell.net")
end