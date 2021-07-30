test_that(".NeverDrop() works", {
  expect_error(.NeverDrop(FALSE, letters))
  expect_equal(character(0), expect_warning(.NeverDrop('NOT THERE', letters)))
  expect_equal(character(0), .NeverDrop(labels = letters))
  expect_equal(letters[2:3], sort(.NeverDrop(c('b', 'c'), labels = letters)))
  expect_equal(letters[2:3], sort(.NeverDrop(3:2, labels = letters)))
  expect_equal(letters[26], expect_warning(.NeverDrop(26:27, labels = letters)))
})
