library(testthat)

test_path <- if (dir.exists("tests/testthat")) "tests/testthat" else "testthat"
test_dir(test_path, load_package = "source")
