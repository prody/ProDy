
read -p "Enter tutorial title (e.g. ENM Analysis using ProDy): " TITLE
while [ -z "$TITLE" ]
do
	read -p "Enter tutorial title: " TITLE
done
read -p "Enter a short title (e.g. ENM Analysis): " SHORT
while [ -z "$SHORT" ]
do
	read -p "Enter a short title: " SHORT
done
read -p "Enter author name (seperate multiple names with comma): " AUTHOR
while [ -z "$AUTHOR" ]
do
	read -p "Enter author name: " AUTHOR
done
SHORT=${SHORT,,}
TUTORIAL=${SHORT// /_}
REFERENCE=${SHORT// /-}


mkdir -p tutorials/$TUTORIAL
mkdir -p tutorials/$TUTORIAL/$TUTORIAL\_files
ln -fs ../../_pkginv tutorials/$TUTORIAL
ln -fs ../../_static tutorials/$TUTORIAL
ln -fs ../../_templates tutorials/$TUTORIAL
ln -fs ../../funding.rst tutorials/$TUTORIAL/acknowledgments.rst
ln -fs ../template/Makefile tutorials/$TUTORIAL
cp -f tutorials/template/Makefile tutorials/$TUTORIAL
sed 's/AUTHOR/'"$AUTHOR"'/g' tutorials/template/conf.py > tutorials/$TUTORIAL/conf.py
sed -i 's/TITLE/'"$TITLE"'/g' tutorials/$TUTORIAL/conf.py
sed 's/TITLE/'"$TITLE"'/g' tutorials/template/index.rst > tutorials/$TUTORIAL/index.rst
sed -i 's/REFERENCE/'"$REFERENCE"'/g' tutorials/$TUTORIAL/index.rst
sed 's/TUTORIAL/'"$TUTORIAL"'/g' tutorials/template/intro.rst > tutorials/$TUTORIAL/intro.rst
echo
echo "Tutorial folders and files are prepared, see tutorials/$TUTORIAL"
