from utils import *

def read_gmt(fpath):
    gdict = {}
    
    with open(fpath, 'r') as f:
        lines = f.readlines()
        for line in lines:
            items = line.strip().rstrip().split('\t')
            gset = items[0].lower()
            genes = items[2:]
            gdict[gset] = genes
            
    return gdict

def write_gmt(gdict, fpath):
    
    with open(fpath, 'w') as f:
        for key, val in gdict.items():
            val_string = '\t'.join(val)
            f.write('%s\thttp\t%s\n'%(key, val_string))
            
    print('finished writing %d gene sets to %s'%(len(gdict), fpath))

if __name__ == "__main__":

    wkdir = '/bgfs/alee/chelsea/projects/10X/AL1/codes'
    os.chdir(wkdir)

    # write gmt
    keyword = 'junction' ### 
    gdict = read_gmt('../db/msigdb.v7.0.symbols.gmt')
    gmt_path = '../data/GSVA/%s.gmt'%(keyword)
    sel_sets = [key for key in gdict.keys() if keyword in key.lower()]
    sel_dict = {}
    for key in sel_sets:
        sel_dict[key] = gdict[key]
    write_gmt(sel_dict, gmt_path)

    # junction, all samples (3000 genes X 8 cell lines)

    for keyword in ['h.all.v6.2.symbols', 'junction']:
        for sample in ['All', 'T47D_WT_KO']:
            for method in ['gsva', 'ssgsea']:

                gmt_path = '../data/GSVA/%s.gmt'%keyword
                expr_path = '../data/Counts/%s.logNormSplicedCts.csv'%sample
                gsvaPrefix = '../data/GSVA/{}_{}'.format(sample, keyword)

                if not isfile('{}.{}.csv'.format(gsvaPrefix, method)):
                    cmd = '/ihome/crc/install/gcc-8.2.0/r/3.6.0/bin/Rscript ./GSVA.R %s %s %s %s'%(expr_path, gmt_path, gsvaPrefix, method)
                    code, stdout, stderr = run(cmd)
                    logg.info('{}, {}, {}'.format(code, stdout, stderr))

# srun -n 1 --cpus-per-task 8 --mem=10g -t 24:00:00 --pty zsh
# python calcGSVA.py > ../logs/calcGSVA.log 2>&1 &