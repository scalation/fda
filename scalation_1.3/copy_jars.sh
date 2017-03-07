# copy_jars.sh - copy jar files from subprojects to lib directories

echo // cp scalation_mathstat/target/scala-2.12/scalation_mathstat_2.12-1.3.jar scalation_modeling/lib
echo // cp scalation_mathstat/target/scala-2.12/scalation_mathstat_2.12-1.3.jar scalation_models/lib
echo // cp scalation_modeling/target/scala-2.12/scalation_modeling_2.12-1.3.jar scalation_models/lib

cp scalation_mathstat/target/scala-2.12/scalation_mathstat_2.12-1.3.jar scalation_modeling/lib
cp scalation_mathstat/target/scala-2.12/scalation_mathstat_2.12-1.3.jar scalation_models/lib
cp scalation_modeling/target/scala-2.12/scalation_modeling_2.12-1.3.jar scalation_models/lib

