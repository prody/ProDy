# rebuild extension modules inplace
rm -rf lib/prody/*/*.so
rm -rf lib/prody/*/*.pyd
python setup.py build_ext --inplace
