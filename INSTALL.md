# install dependencies 

On a fresh Ubuntu 14.04 installation, install homebrew and python3.7

    sh -c "$(curl -fsSL https://raw.githubusercontent.com/Linuxbrew/install/master/install.sh)"
    
    export PATH=/home/linuxbrew/.linuxbrew/bin:$PATH
    
    python3 --version
    # Python 3.5.2
    # python3.5.2 is broken: https://stackoverflow.com/a/56010650/178297. 
    # Need to upgrade to >=3.5.3
    # pyenv is a convenient way to do this: https://github.com/pyenv/pyenv
    sudo apt-get install -y make build-essential libssl-dev zlib1g-dev libbz2-dev \
    libreadline-dev libsqlite3-dev wget curl llvm libncurses5-dev libncursesw5-dev \
    xz-utils tk-dev libffi-dev liblzma-dev python-openssl git
    brew install pyenv
    eval "$(pyenv init -)"
    pyenv install 3.5.3
    pyenv local 3.5.3
    python --version
    # Python 3.5.3
    pip install --upgrade pip

    # merge DNAnexus' PYTHONPATH with this from PYTHON3
    export PYTHONPATH=/home/linuxbrew/.linuxbrew/lib/python3.7/site-packages:/usr/share/dnanexus/lib/python2.7/site-packages


Then install the system dependencies: freebayes (>=1.2.0), htslib (tabix+bgzip), gsort.        
    brew tap brewsci/bio
    brew install freebayes
    brew install htslib
    
    curl -s https://api.github.com/repos/brentp/gsort/releases/latest \
      | grep browser_download_url \
      | grep -i $(uname) \
      | cut -d '"' -f 4 \
      | wget -O gsort -qi -
    chmod +x gsort
    export PATH=.:$PATH

# install mity

Either install mity globally:

    export PYTHONPATH=/usr/share/dnanexus/lib/python2.7/site-packages
    export PYTHONPATH=/usr/local/lib/python3.5/dist-packages:/usr/lib/python3/dist-packages:/usr/share/dnanexus/lib/python2.7/site-packages
    
    # fix a python version incompatibility bug in futures
    sudo perl -pi -e 's|raise exception_type, self._exception, self._traceback|raise Exception(self._exception).with_traceback(self._traceback)|' /usr/share/dnanexus/lib/python2.7/site-packages/concurrent/futures/_base.py
    
    VERSION=0.0.1b7
    pip3 install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple mity==$VERSION
    
Or install mity using a virtualenv

    sudo apt-get install python3-venv
    unset PYTHONPATH
    VERSION=0.0.1a9
    python3 -m venv .
    source bin/activate
    ./bin/pip install wheel
    ./bin/pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple mity==$VERSION

# test
(These URLs valid until 26/7/2020)
    wget https://dl.dnanex.us/F/D/XJfjx2X139ZkzY7b29QQKBppzfj9p5V794Bfqf4G/A1.dedup.realigned.recalibrated.chrMT.bam
    wget https://dl.dnanex.us/F/D/qyV40Qgfj6Jgy3zZfJ07vkgXqZvJ6Fb2kXb24fyv/A1.dedup.realigned.recalibrated.chrMT.bam.bai
    wget https://dl.dnanex.us/F/D/pVG7PjZy4qKBB6ZKbkkF0X6kB0kxf7ZzjpK7fXjY/hs37d5.fasta-index.tar.gz
    tar -xzvf hs37d5.fasta-index.tar.gz; mv genome.dict hs37d5.dict; mv genome.fa hs37d5.fa; mv genome.fa.fai hs37d5.fa.fai

    
    
