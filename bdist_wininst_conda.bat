del /S /Q /F build\*.*
python setup.py bdist_wininst
activate py35
python setup.py bdist_wininst
activate py27
python setup.py bdist_wininst