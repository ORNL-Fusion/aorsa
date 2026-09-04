cd test/DIIID-whistler
resultA=$(bash -c '(./test.sh); exit $?' 2>&1)
exitA=$?
echo codes $resultA $exitA
cd -

pwd
cd test/DIIID-helicon
resultA=$(bash -c '(./test.sh); exit $?' 2>&1)
exitA=$?
echo codes $resultA $exitA
cd -


