SET VS90COMNTOOLS=%VS120COMNTOOLS%
python setup.py build_ext -i
pyinstaller vipster.spec