
read -p "Enter tutorial title (e.g. ENM Analysis using ProDy): " TITLE
read -p "Enter a short title (e.g. enm analysis): " SHORT
SHORT=${SHORT,,}
TUTORIAL=${SHORT// /_}
REFERENCE=${SHORT// /-}


mkdir -p tutorials/$TUTORIAL
mkdir -p tutorials/$TUTORIAL/$TUTORIAL_files
ln -fs ../../_static tutorials/$TUTORIAL
cp -f tutorials/template/conf.py tutorials/$TUTORIAL
cp -f tutorials/template/Makefile tutorials/$TUTORIAL
sed 's/TITLE/'"$TITLE"'/g' tutorials/template/index.rst > tutorials/$TUTORIAL/index.rst 
sed -i 's/REFERENCE/'"$REFERENCE"'/g' tutorials/$TUTORIAL/index.rst
sed 's/TUTORIAL/'"$TUTORIAL"'/g' tutorials/template/intro.rst > tutorials/$TUTORIAL/intro.rst
echo
echo "Tutorial folders and files are prepared, see tutorials/$TUTORIAL"
