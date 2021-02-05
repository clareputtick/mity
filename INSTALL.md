# Docker
The simplest way to run `mity` is via `docker`:
```
docker run drmjc/mity -h
```

# pip
If you have `freebayes` >=1.2 and Brent Pederson's `gsort` installed, then `pip` should work well
```
pip3 install mitywgs
```

# manual installation 
If you would prefer to install `mity` on a fresh Ubuntu installation, the following should work.
We have tested this on a fresh Ubuntu 14.04 image; We use `pyenv` to install python 3.7.4, though there
are a number of alternatives. YMMV.

## install dependencies 
### install `homebrew` and `python3.7.4`
* `pyenv` is a convenient way to manage multiple python distros: https://github.com/pyenv/pyenv
```
sh -c "$(curl -fsSL https://raw.githubusercontent.com/Linuxbrew/install/master/install.sh)"
export PATH=/home/linuxbrew/.linuxbrew/bin:$PATH

sudo apt-get install -y make build-essential libssl-dev zlib1g-dev libbz2-dev \
  libreadline-dev libsqlite3-dev wget curl llvm libncurses5-dev libncursesw5-dev \
  xz-utils tk-dev libffi-dev liblzma-dev python-openssl git
brew install pyenv
eval "$(pyenv init -)"
pyenv install 3.7.4
pyenv local 3.7.4
python --version
# Python 3.7.4
pip3 install --upgrade pip
export PATH=$PATH:.local/bin:$HOME/.pyenv/versions/3.7.4/bin
# if running on a DNANexus cloud instance, then merge DNAnexus' PYTHONPATH with this from PYTHON3
export PYTHONPATH=/home/linuxbrew/.linuxbrew/lib/python3.7/site-packages:/usr/share/dnanexus/lib/python2.7/site-packages
```

### install dependencies: `freebayes` (>=1.2.0), `htslib` (tabix+bgzip), `gsort`, 'tabix'.
```
brew tap brewsci/bio
brew install freebayes
brew install htslib
sudo apt-get install -y tabix

curl -s https://api.github.com/repos/brentp/gsort/releases/latest \
  | grep browser_download_url \
  | grep -i $(uname) \
  | cut -d '"' -f 4 \
  | wget -O gsort -qi -
chmod +x gsort
export PATH=.:$PATH
```

### install `mity`

#### Either install `mity` globally:
```
# for most users
export PYTHONPATH=/usr/local/lib/python3.7/dist-packages:/usr/lib/python3/dist-packages
# for those using a DNANexus cloud instance
export PYTHONPATH=/usr/local/lib/python3.7/dist-packages:/usr/lib/python3/dist-packages:/usr/share/dnanexus/lib/python2.7/site-packages

# fix a python version incompatibility bug in futures
sudo perl -pi -e 's|raise exception_type, self._exception, self._traceback|raise Exception(self._exception).with_traceback(self._traceback)|' /usr/share/dnanexus/lib/python2.7/site-packages/concurrent/futures/_base.py

pip3 install mitywgs
```

#### Or install `mity` using a `virtualenv`
```
sudo apt-get install python3-venv
unset PYTHONPATH
python3 -m venv .
source bin/activate
./bin/pip install wheel
./bin/pip install mitywgs
```

# test mity
* These URLs valid until 26/7/2020
```
wget https://dl.dnanex.us/F/D/XJfjx2X139ZkzY7b29QQKBppzfj9p5V794Bfqf4G/A1.dedup.realigned.recalibrated.chrMT.bam
wget https://dl.dnanex.us/F/D/qyV40Qgfj6Jgy3zZfJ07vkgXqZvJ6Fb2kXb24fyv/A1.dedup.realigned.recalibrated.chrMT.bam.bai
mity call --normalise A1.dedup.realigned.recalibrated.chrMT.bam
mity report A1.dedup.realigned.recalibrated.chrMT.mity.vcf.gz
```

## test using docker
```
docker run drmjc/mity:0.2.0 -h
wget https://dl.dnanex.us/F/D/XJfjx2X139ZkzY7b29QQKBppzfj9p5V794Bfqf4G/A1.dedup.realigned.recalibrated.chrMT.bam
wget https://dl.dnanex.us/F/D/qyV40Qgfj6Jgy3zZfJ07vkgXqZvJ6Fb2kXb24fyv/A1.dedup.realigned.recalibrated.chrMT.bam.bai
docker run --rm -it -v $(pwd):/home drmjc/mity:0.2.0 call --prefix A1 A1.dedup.realigned.recalibrated.chrMT.bam
docker run --rm -it -v $(pwd):/home drmjc/mity:0.2.0 normalise --outfile A1.mity.norm.vcf.gz A1.mity.vcf.gz
docker run --rm -it -v $(pwd):/home drmjc/mity:0.2.0 report A1.mity.norm.vcf.gz
```
