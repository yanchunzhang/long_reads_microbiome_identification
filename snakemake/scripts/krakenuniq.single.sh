
prefix=$3
classified_out="$prefix.classified_by_krakenuniq"
input=$2
DB=$1
thread=$4
krakenuniq --preload --db $DB --threads $thread --classified-out $classified_out --report $prefix.krakenuniq.report $input > $prefix.krakenuniq

