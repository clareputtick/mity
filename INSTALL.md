# install dependencies 

    brew tap brewsci/bio
    brew install freebayes
    
    brew install htslib
    
    curl -s https://api.github.com/repos/brentp/gsort/releases/latest \
      | grep browser_download_url \
      | grep -i $(uname) \
      | cut -d '"' -f 4 \
      | wget -O gsort -qi -
    chmod +x gsort

# install mity

    VERSION=0.0.1a4
    python -m venv .
    source bin/activate
    bin/pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple mity==$VERSION