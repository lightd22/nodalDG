#!/bin/sh
echo 'What are maxPolyDegree Values?'
read -a maxPolyVals
echo 'What are startRes Values?'
read -a startResVals

echo "maxPolyDegree Values: ${maxPolyVals[*]}"
echo "startRes Values: ${startResVals[*]}"

ind=${!maxPolyVals[*]}
for i in $ind;
do
  echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  echo "  Starting run using maxPolyDegree = ${maxPolyVals[$i]}"
  echo "  startRes = ${startResVals[$i]}"
  echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  newMaxPoly="maxPolyDegree = ${maxPolyVals[$i]}"
  newStartRes="startRes = ${startResVals[$i]}"

  sed -i '' 's/maxPolyDegree = ./'"$newMaxPoly"'/' inputs.nl
  sed -i '' 's/startRes = [0-9][0-9]/'"$newStartRes"'/' inputs.nl

  make 2d_test
done
