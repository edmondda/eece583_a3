#! /bin/tcsh -f


echo
echo ""Enhanced  "Kernighan-Lin  based  bi-partitioning algorithm for hypergraphs"
echo

# create out folder, if not existing
if (!(-d out)) then
  mkdir out
endif


# for all infiles
foreach inFile ($1/*.txt)
  set designName=`echo $inFile|cut -d. -f1|rev|cut -d/ -f1|rev`
  echo " Run design: $designName ..."
  ./assignment3 "$designName.txt" 
  #./assignment3 "$designName.txt"  -g -r f
end

unset designName
