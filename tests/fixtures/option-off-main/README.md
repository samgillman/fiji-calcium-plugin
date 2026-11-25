# Files main writes with the option off and no ring files

What the plugin of release 2026.11.24.1 with pull request 35 (main at b224305) saved beside the tables for the steps of the
option-off part of `tests/test_review_marks.groovy`: an analysis, two manual marks, then a second analysis that carries them
over. The time of each mark is blanked as `TIME`. `messages.txt` holds the finish messages and the Methods sentences in the
order the test reads them.

The test runs the same steps with this build and compares. Every file must be the same, apart from the new last column
`QC_flags` in the call checks tables, empty without a flag. The review files are the same byte for byte: with manual
marks only they have no `Marked_by` column.
