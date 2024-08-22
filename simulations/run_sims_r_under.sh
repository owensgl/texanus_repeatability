rep=$1

for recomb in 1e-3 5e-4 1e-4 5e-5 1e-5 5e-6 1e-6 5e-7 1e-7
do
	slim -d selection=0.5 -d recom=$recomb BC1_nonWF_underdominant.v0.2.slim | grep "^chr\|^01" > underdom/output.$rep.$recomb.txt
	python3 ../cvtk/run_cvtk.py underdom/output.$rep.$recomb.txt > underdom/conv.re.$rep.$recomb.txt
done
