for d in */; do
    n=${d%/}
    if [[ "$n" =~ ^[0-9.eE+-]+$ ]] && awk "BEGIN {exit !($n > 2e7)}"; then
        rm -rf -- "$d"
    fi
done