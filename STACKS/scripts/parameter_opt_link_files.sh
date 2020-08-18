### Make directory a variable
dir=$(pwd)
dir2="/03_clone_filter_out/"

### Create links
while IFS= read -r line
do
ln -s $dir$dir2$line* ./05_parameter_opt/clone_filter_test/
done < hs.txt
