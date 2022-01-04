from utils import *
import dask
from dask.diagnostics import ProgressBar

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

import pyscenic 
from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell

def calcTFs(expr, tf_names, db, prefix,
            motif_path='../data/pySCENIC/ref/motifs-v9-nr.hgnc-m0.001-o0.0.tbl', 
            out_path='../data/pySCENIC', 
            ppn=8):
    """Computes motifs, regulons and trancriptional factor activation using pySCENIC.

    Arguments
    ---------
    expr: `pandas DataFrame` 
        cell X gene raw counts; FPKM; not TPM as coexpression will be calculated
    tf_names: `list` (`str`)
        curated human transcriptional factor downloaded from github: pySCENIC/ref/hs_hgnc_curated_tfs.txt
    db: `list` (`FeatherRankingDatabase()`)
        feather files, ranking genome [FeatherRankingDatabase(name="hg38__refseq-r80__10kb_up_and_down_tss")]
    prefix: `str` (default: `None`)
        Specify name to save files (eg, cell line names)

    Returns
    -------
    Do not return but write files (the calc takes too long...)
    """

    # Inference of co-expression modules
    adjacencies = grnboost2(expr, tf_names=tf_names, verbose=True)
    modules = list(modules_from_adjacencies(adjacencies, expr))

    # Calculate a list of enriched motifs and the corresponding target genes for all modules.
    with ProgressBar():
        df = prune2df(db, modules, motif_path, num_workers=ppn)

    # Create regulons from this table of enriched motifs.
    regulons = df2regulons(df)

    # Save the enriched motifs and the discovered regulons to disk.
    with open('{}/{}_motifs.csv'.format(out_path, prefix), "wb") as f:
        pickle.dump(regulons, f)

    auc_mtx = aucell(expr, regulons, num_workers=ppn)
    tfs = [tf.strip('(+)') for tf in auc_mtx.columns]
    auc_mtx.to_csv('{}/{}_auc_mtx.csv'.format(out_path, prefix))

    print('finished calculation for %s'%(prefix))

if __name__ == "__main__":

    wkdir = '/bgfs/alee/chelsea/projects/10X/AL1/codes'
    os.chdir(wkdir)

    db = [RankingDatabase(fname='../data/pySCENIC/ref/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather', name='hg38__refseq-r80__10kb_up_and_down_tss.mc9nr')]
    tf_names = load_tf_names('../data/pySCENIC/ref/hs_hgnc_curated_tfs.txt')

    CellTypes = ['MCF7','T47D WT','T47D KO',  'MM134','SUM44','BCK4', 'MCF10A','HEK293'] 

    for cell in CellTypes:
        
        tmp = adata_raw[adata_raw.obs['CellType'] == cell]
        RawSplicedCts = pd.DataFrame(tmp.layers['spliced'].todense(), index=tmp.obs.index, columns=tmp.var.index) # cell X gene
        print(cell, RawSplicedCts.shape)
        
        if not isfile('{}/{}_auc_mtx.csv'.format('../data/pySCENIC', cell)):      
            calcTFs(RawSplicedCts, tf_names, db, cell)
            print('calculating TF for %s'%(cell))

# srun -n 1 --cpus-per-task 8 --mem=10g -t 24:00:00 --pty zsh
# python calcTF.py > ../logs/calcTF.log 2>&1 &
# tail -n 5 ../logs/calcTF.log