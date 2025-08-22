## sh step3.pyscenic_from_loom.sh > pyscenic.log 2>&1
## sh step3.pyscenic_from_loom.sh -i ./02.result/mouse.loom > process.log 2>&1
## nohup sh step3.pyscenic_from_loom.sh > pyscenic.log 2>&1 &

#default value
# input_loom=${1}
# output_path=${2}
input_loom='/data/work/02.result/others/ST_NC/02.result/05.scenic/anno_0818_5w_1wHVG.loom'
output_path='/data/work/02.result/others/ST_NC/02.result/05.scenic'
n_workers=20
#help function
function usage() {
echo -e "OPTIONS:\n-i|--input_loom:\t input loom file"
echo -e "-n|--n_workers:\t working core number"
echo -e "-h|--help:\t Usage information"
exit 1
}
#get value
while getopts :i:n:h opt
do
    case "$opt" in
        i) input_loom="$OPTARG" ;;
        n) n_workers="$OPTARG" ;;
        h) usage ;;
        :) echo "This option -$OPTARG requires an argument."
           exit 1 ;;
        ?) echo "-$OPTARG is not an option"
           exit 2 ;;
    esac
done
# 需要更改路径
tfs=/data/work/02.result/others/ST_NC/02.result/05.scenic/database/allTFs_hg38.txt
feather=/data/work/02.result/others/ST_NC/02.result/05.scenic/database/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
tbl=/data/work/02.result/others/ST_NC/02.result/05.scenic/database/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl

pyscenic=pyscenic

# grn
$pyscenic grn \
--num_workers $n_workers \
--output ${output_path}/grn.tsv \
--method grnboost2 \
$input_loom  $tfs

# cistarget
$pyscenic ctx \
${output_path}/grn.tsv $feather \
--annotations_fname $tbl \
--expression_mtx_fname $input_loom \
--mode "dask_multiprocessing" \
--output ${output_path}/ctx.csv \
--num_workers $n_workers   \
--mask_dropouts

# AUCell
$pyscenic aucell \
$input_loom \
${output_path}/ctx.csv \
--output ${output_path}/aucell.loom \
--num_workers $n_workers
