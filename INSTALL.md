# install dependencies 

On a fresh Ubuntu 14.04 installation, install homebrew and python3.7

    sh -c "$(curl -fsSL https://raw.githubusercontent.com/Linuxbrew/install/master/install.sh)"
    export PATH=/home/linuxbrew/.linuxbrew/bin:$PATH
    export PYTHONPATH=/home/linuxbrew/.linuxbrew/lib/python3.7/site-packages
    python3 --version
    # Python 3.7.4

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

    VERSION=0.0.1a4
    pip3 install wheel
    pip3 install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple mity==$VERSION
    
Or install mity using a virtualenv

    sudo apt-get install python3-venv
    unset PYTHONPATH
    VERSION=0.0.1a4
    python3 -m venv .
    source bin/activate
    ./bin/pip install wheel
    ./bin/pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple mity==$VERSION
