# installation

mkdir $HOME/bin
export PATH=$HOME/bin:$PATH

# install samtools
wget https://github.com/samtools/samtools/releases/download/1.6/samtools-1.6.tar.bz2
tar jxvf samtools-1.6.tar.bz2
cd samtools-1.6
make
cd ..
cp samtools-1.6/samtools $HOME/bin

# install hisat2
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.0.1-beta-OSX_x86_64.zip
unzip hisat2-2.0.1-beta-OSX_x86_64.zip
cp hisat2-2.0.1-beta/hisat2* hisat2-2.0.1-beta/*.py $HOME/bin

# download hisat index
mkdir indexes
cd indexes
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38_tran.tar.gz

# install bowtie
# manually downloaded bowtie2-2.2.9-macos-x86_64
unzip bowtie2-2.2.9-macos-x86_64.zip
cp bowtie2-2.2.9/bowtie* $HOME/bin
