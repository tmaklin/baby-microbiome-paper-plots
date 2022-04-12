niters=10000
~/Projects/software/fastspar/src/fastspar -y --otu_table $1 --iterations 1000 -x 200 --correlation $2""_correlation.tsv --covariance $2""_covariance.tsv --threads 4
mkdir tmp_bootstrap
~/Projects/software/fastspar/src/fastspar_bootstrap --otu_table $1  --number $niters --prefix tmp_bootstrap/$2 --threads 4
mkdir tmp_correlation
parallel -j 4 ~/Projects/software/fastspar/src/fastspar -y --otu_table {} --correlation tmp_correlation/cor_{/} --covariance tmp_correlation/cov_{/} -i 5 ::: tmp_bootstrap/*
~/Projects/software/fastspar/src/fastspar_pvalues --otu_table $1 --correlation $2""_correlation.tsv --prefix tmp_correlation/cor_ --permutations $niters --outfile $2""_pvalues.tsv --threads 4
rm -rf tmp_bootstrap
rm -rf tmp_correlation
