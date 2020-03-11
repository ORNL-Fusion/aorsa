add_test(Helicon "./test.sh")
set_tests_properties(Helicon PROPERTIES  WORKING_DIRECTORY "test/DIIID-helicon")
add_test(Whistler "./test.sh")
set_tests_properties(Whistler PROPERTIES  WORKING_DIRECTORY "test/DIIID-whistler")
