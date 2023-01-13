'Usage:
  msstats.R --input=<DIR> --output=<DIR> --anno=<FILE> [--log]
  msstats.R (-h | --help)
  msstats.R --version

Options:
  --input=<DIR>    MaxQuant Input Folder 
  --output=<DIR>   MSstats Output Folder
  --anno=<FILE>    Experiment Designe Table
  --log            Keep Log Files
  -h --help        Show this screen.
  --version        Show version.
' -> doc

library(docopt)
arguments <- docopt(doc, version = paste('MSstats', packageVersion('MSstats'), '\n'))
print(arguments)
arguments$input
arguments$output
arguments$anno
arguments$log

