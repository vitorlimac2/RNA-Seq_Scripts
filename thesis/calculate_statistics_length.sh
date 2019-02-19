# MAPPED READ INFO
REP=$1

echo "Comprimento de reads"
cat "$REP" | awk '$2 < 50' | wc -l
cat "$REP" | awk '$2 >= 50 && $2 < 100' | wc -l
cat "$REP" | awk '$2 >= 100 && $2 < 150' | wc -l
cat "$REP" | awk '$2 >= 150 && $2 < 200' | wc -l
cat "$REP" | awk '$2 >= 200 && $2 < 250' | wc -l
cat "$REP" | awk ' $2 >= 250' | wc -l

echo "Mapeamentos unicos"
cat "$REP" | awk '$3== 1 && $2 < 50' | wc -l
cat "$REP" | awk '$3== 1 && $2 >= 50 && $2 < 100' | wc -l
cat "$REP" | awk '$3== 1 && $2 >= 100 && $2 < 150' | wc -l
cat "$REP" | awk '$3== 1 && $2 >= 150 && $2 < 200' | wc -l
cat "$REP" | awk '$3== 1 && $2 >= 200 && $2 < 250' | wc -l
cat "$REP" | awk '$3== 1 && $2 >= 250' | wc -l


echo "Mapeamentos totais"
cat "$REP" | awk '$3!= 0 && $2 < 50' | wc -l
cat "$REP" | awk '$3!= 0 && $2 >= 50 && $2 < 100' | wc -l
cat "$REP" | awk '$3!= 0 && $2 >= 100 && $2 < 150' | wc -l
cat "$REP" | awk '$3!= 0 && $2 >= 150 && $2 < 200' | wc -l
cat "$REP" | awk '$3!= 0 && $2 >= 200 && $2 < 250' | wc -l
cat "$REP" | awk '$3!= 0 && $2 >= 250' | wc -l

