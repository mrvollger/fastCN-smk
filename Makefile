# change LDFLAGS to suit your needs, for example below is for the eichler lab
# LDFLAGS=-Wl,-rpath=/net/gs/vol3/software/modules-sw/gcc/8.2.0/Linux/CentOS7/x86_64/lib64

QuicK-mer2/quicKmer2: fastCN/GC_control_gen
	git clone https://github.com/KiddLab/QuicK-mer2.git \
		&& cd QuicK-mer2 \
		&& gcc -LDFLAGS=${LDFLAGS} QuicKmer.c -O3 -g -pthread -std=c99 -lm -o quicKmer2 \
		&& cd .. \
		&& ln -s ../QuicK-mer2/quicKmer2 bin/




fastCN/GC_control_gen:
	git clone https://github.com/KiddLab/fastCN.git \
		&& cd fastCN \
		&& sed -i -e 's/Use_strict 1/Use_strict 0/g' GC_control_gen.cc \
		&& g++ ${LDFLAGS} -o GC_control_gen GC_control_gen.cc \
		&& g++ ${LDFLAGS} -o SAM_GC_correction SAM_GC_correction.cc \
		&& cd .. \
		&& ln -s fastCN bin \
		&& cargo install rustybam --root .




.PHONY: clean

clean:
	rm -rf fastCN && rm -rf QuicK-mer2
