. ./env/bin/activate
. ~/.bash_functions
CURR_VERSION=$(awk '/version="/{print $1}' setup.py | tr -d '[version=",]')
sed -i 's|'"$CURR_VERSION"'|'"$1"'|g' setup.py
python setup.py bdist_wheel
twine upload --skip-existing dist/*

commit $PWD
